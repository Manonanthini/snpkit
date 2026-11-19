"""Batch annotation / liftover of a text file of variants against the local db."""

import logging
import sqlite3

import pandas as pd

from .common import norm_chrom, revcomp
from .genes import GeneIndex
from .lift import Lifter

log = logging.getLogger("snpkit")

ALIASES = {
    "rsid": ["rsid", "snp", "id", "rs", "variant_id", "markername", "snpid", "rsids"],
    "chr": ["chr", "chrom", "chromosome", "#chrom", "chr_name"],
    "pos": ["pos", "bp", "position", "base_pair_location", "chr_position", "bp_pos"],
    "ref": ["ref", "reference", "ref_allele", "other_allele", "non_effect_allele"],
    "alt": ["alt", "alt_allele", "effect_allele", "a1"],
}


def read_variants(path, colmap=None):
    """Read whitespace/tab text file; return (df_raw, df_std) with std cols
    rsid, chr, pos, ref, alt (missing ones are None).  Single-column files may
    hold rsIDs or chr:pos[:ref:alt] strings."""
    raw = pd.read_csv(path, sep=r"\s+", dtype=str, comment=None, engine="python")
    raw.columns = [c.lstrip("#") if c.startswith("#CHROM") else c for c in raw.columns]
    low = {c.lower(): c for c in raw.columns}
    colmap = colmap or {}
    std = pd.DataFrame(index=raw.index)
    found = False
    for key, names in ALIASES.items():
        src = colmap.get(key) or next((low[n] for n in names if n in low), None)
        std[key] = raw[src] if src else None
        found |= src is not None
    if not found:
        # headerless single column: first token was parsed as header, re-read
        raw = pd.read_csv(path, sep=r"\s+", dtype=str, header=None, engine="python")
        col = raw.iloc[:, 0]
        std = pd.DataFrame(index=raw.index, columns=list(ALIASES), dtype=object)
        is_cp = col.str.contains(":", regex=False)
        parts = col.where(is_cp).str.split(":", expand=True)
        std.loc[~is_cp, "rsid"] = col[~is_cp]
        if is_cp.any():
            std.loc[is_cp, "chr"] = parts[0]
            std.loc[is_cp, "pos"] = parts[1]
            if parts.shape[1] > 3:
                std.loc[is_cp, "ref"] = parts[2]
                std.loc[is_cp, "alt"] = parts[3]
        raw.columns = ["input"]
    std["chr"] = std["chr"].map(lambda c: norm_chrom(c) if isinstance(c, str) else None)
    std["pos"] = pd.to_numeric(std["pos"], errors="coerce").astype("Int64")
    for k in ("ref", "alt"):
        std[k] = std[k].map(lambda a: a.upper() if isinstance(a, str) else None)
    return raw, std


def _allele_status(qref, qalt, r, a):
    if not qref or not qalt or not r or not a:
        return "no_alleles"
    if (qref, qalt) == (r, a):
        return "exact"
    if (qref, qalt) == (a, r):
        return "swapped"
    fr, fa = revcomp(qref), revcomp(qalt)
    if (fr, fa) == (r, a):
        return "strand_flip"
    if (fr, fa) == (a, r):
        return "strand_flip_swapped"
    return "mismatch"


_RANK = {"exact": 0, "swapped": 1, "strand_flip": 2, "strand_flip_swapped": 3,
         "no_alleles": 4, "mismatch": 5}

OUT_COLS = ["rsid", "chr", "pos37", "pos38", "ref37", "alt37", "ref38", "alt38", "gene"]


def annotate(db, path, build, colmap=None, chain=None, chunk=500_000):
    """Annotate variants given in `build` (37 or 38). Returns DataFrame."""
    assert build in (37, 38)
    other = 75 - build
    raw, std = read_variants(path, colmap)
    std = std.copy()
    std["qi"] = range(len(std))

    con = sqlite3.connect(str(db))
    con.execute("DROP TABLE IF EXISTS temp.q")
    con.execute(
        "CREATE TEMP TABLE q(qi INTEGER, rsid TEXT, chr TEXT, pos INTEGER, ref TEXT, alt TEXT)"
    )
    con.executemany(
        "INSERT INTO q VALUES(?,?,?,?,?,?)",
        [
            (r.qi, r.rsid if isinstance(r.rsid, str) else None, r.chr,
             None if pd.isna(r.pos) else int(r.pos), r.ref, r.alt)
            for r in std.itertuples()
        ],
    )
    con.execute("CREATE INDEX temp.q_r ON q(rsid)")
    con.execute("CREATE INDEX temp.q_p ON q(chr,pos)")

    sel = ("SELECT q.qi, q.ref AS qref, q.alt AS qalt, v.rsid, v.chr, v.pos37, v.pos38, "
           "v.ref37, v.alt37, v.ref38, v.alt38, v.gene FROM q JOIN variants v ON ")
    by_rsid = pd.read_sql_query(sel + "v.rsid=q.rsid", con)
    by_rsid["source"] = "db_rsid"
    hit = set(by_rsid.qi)
    by_pos = pd.read_sql_query(
        sel + f"v.chr=q.chr AND v.pos{build}=q.pos WHERE q.chr IS NOT NULL", con
    )
    by_pos = by_pos[~by_pos.qi.isin(hit)]
    by_pos["source"] = "db_pos"
    m = pd.concat([by_rsid, by_pos], ignore_index=True)

    rb, ab = f"ref{build}", f"alt{build}"
    m["allele_status"] = [
        _allele_status(a, b, c, d) for a, b, c, d in zip(m.qref, m.qalt, m[rb], m[ab])
    ]
    # Position-based hits with alleles supplied: drop rows that disagree entirely
    has_q = m.qref.notna() & m.qalt.notna()
    bad = (m.source == "db_pos") & has_q & (m.allele_status == "mismatch")
    m = m[~bad]
    m["_rank"] = m.allele_status.map(_RANK)
    m = m.sort_values(["qi", "_rank"])
    # Queries with alleles: keep the single best-matching db row.
    # Queries without alleles: keep every ALT row (multi-allelic expansion).
    with_alleles = set(std.qi[std.ref.notna() & std.alt.notna()])
    best = pd.concat(
        [m[m.qi.isin(with_alleles)].drop_duplicates("qi"), m[~m.qi.isin(with_alleles)]]
    ).sort_values("qi", kind="stable")

    out = std[["qi"]].merge(best.drop(columns=["_rank", "qref", "qalt"]), on="qi", how="left")
    out["source"] = out["source"].fillna("none")
    out["allele_status"] = out["allele_status"].fillna("not_found")

    # Fallback: chain-file liftover for positions not in db
    miss = out[(out.source == "none")].index
    if len(miss) and chain:
        lifter = Lifter(chain)
        gi = GeneIndex(
            [(c, s, e, n) for c, s, e, n in
             con.execute("SELECT chr,start,end,name FROM genes")]
        )
        for i in miss:
            q = std.iloc[out.at[i, "qi"]]
            if q.chr is None or pd.isna(q.pos):
                continue
            hit = lifter.lift(q.chr, int(q.pos))
            if hit is None:
                out.at[i, "allele_status"] = "unmapped"
                continue
            tc, tp, strand = hit
            out.at[i, "chr"] = tc
            out.at[i, f"pos{build}"] = int(q.pos)
            out.at[i, f"pos{other}"] = tp
            g38 = tp if other == 38 else int(q.pos)
            out.at[i, "gene"] = gi.genes_at(tc, g38)
            out.at[i, "source"] = "chain" if strand == "+" else "chain_minus_strand"
            out.at[i, "allele_status"] = "chain_no_allele_check"
    con.close()

    for c in ("pos37", "pos38"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    ann = out[["qi"] + OUT_COLS + ["source", "allele_status"]].rename(
        columns={c: f"out_{c}" for c in OUT_COLS + ["source", "allele_status"]}
    )
    raw = raw.reset_index(drop=True)
    raw["qi"] = range(len(raw))
    res = raw.merge(ann, on="qi", how="left").drop(columns="qi")
    n = len(std)
    nf = int((out.source == "none").sum())
    log.info("annotated %s variants (%s not found)", f"{n:,}", f"{nf:,}")
    return res
