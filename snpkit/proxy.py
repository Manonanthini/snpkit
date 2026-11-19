"""LD proxy search.

PLINK (1.9) does the LD arithmetic on a reference panel (1000G / TOPMed / your
own data set) but only for the target SNPs (--ld-snp-list), so there is no
genome-wide pairwise table to precompute.  Results are cached in the db `ld`
table and annotated (b37/b38 position, alleles, gene) from the variants table.

Panel variant IDs in the .bim must be rsIDs.  The signed correlation R is
requested so the allele phase between target and proxy can be reported, which
is needed to carry an effect allele over to the proxy (e.g. for PRS).
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .common import chunked

log = logging.getLogger("snpkit")


def read_bim(prefix, wanted=None):
    """Return {rsid: (chr, pos, a1, a2)} from prefix.bim (optionally filtered)."""
    out = {}
    with open(f"{prefix}.bim") as fh:
        for ln in fh:
            c, rid, _, pos, a1, a2 = ln.split()[:6]
            if wanted is None or rid in wanted:
                out[rid] = (c, int(pos), a1, a2)
    return out


def parse_plink_ld(path):
    """Parse a PLINK 1.9 --r .ld file -> DataFrame(snp_a, snp_b, chr, pos_a, pos_b, R)."""
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df.rename(columns={"SNP_A": "snp_a", "SNP_B": "snp_b", "BP_A": "pos_a",
                            "BP_B": "pos_b", "CHR_A": "chr", "R": "R"})
    return df[["snp_a", "snp_b", "chr", "pos_a", "pos_b", "R"]]


def run_plink(plink, prefix, targets, r2, kb, workdir, threads=1):
    tfile = Path(workdir) / "targets.txt"
    tfile.write_text("\n".join(targets) + "\n")
    outp = Path(workdir) / "ld"
    cmd = [
        plink, "--bfile", prefix, "--r", "--ld-snp-list", str(tfile),
        "--ld-window-kb", str(kb), "--ld-window", "99999",
        "--ld-window-r2", str(r2),
        "--threads", str(threads), "--out", str(outp), "--silent",
    ]
    log.info("running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    ld = Path(f"{outp}.ld")
    if not ld.exists():  # plink writes no file when nothing passes
        return pd.DataFrame(columns=["snp_a", "snp_b", "chr", "pos_a", "pos_b", "R"])
    return parse_plink_ld(ld)


def _directional(ld, targets):
    """Plink lists each pair once; expand to (target -> proxy) rows."""
    t = set(targets)
    a = ld[ld.snp_a.isin(t)].rename(
        columns={"snp_a": "target", "snp_b": "proxy", "pos_a": "pos_target", "pos_b": "pos_proxy"})
    b = ld[ld.snp_b.isin(t)].rename(
        columns={"snp_b": "target", "snp_a": "proxy", "pos_b": "pos_target", "pos_a": "pos_proxy"})
    d = pd.concat([a, b], ignore_index=True).drop_duplicates(["target", "proxy"])
    return d


def find_proxies(db, targets, panel, r2=0.8, kb=500, available=None, top=5,
                 plink="plink", panel_name="panel", panel_build=38, threads=1):
    """Return DataFrame of proxies ranked by r2 for each target rsID.

    panel : bfile prefix; may contain '{chr}' for per-chromosome filesets.
    available : optional set of rsIDs that proxies must belong to (e.g. SNPs
                present in the target data set / summary statistics).
    """
    import sqlite3

    if shutil.which(plink) is None:
        raise RuntimeError(f"'{plink}' not found on PATH (set --plink)")
    targets = list(dict.fromkeys(targets))
    frames, bims = [], {}
    with tempfile.TemporaryDirectory() as wd:
        if "{chr}" in panel:
            con = sqlite3.connect(str(db))
            tchr = dict(con.execute(
                f"SELECT rsid,chr FROM variants WHERE rsid IN ({','.join('?'*len(targets))})",
                targets).fetchall()) if len(targets) < 30000 else {}
            con.close()
            groups = {}
            for t in targets:
                groups.setdefault(tchr.get(t), []).append(t)
            jobs = [(panel.format(chr=c), ts) for c, ts in groups.items() if c]
        else:
            jobs = [(panel, targets)]
        for prefix, ts in jobs:
            frames.append(run_plink(plink, prefix, ts, r2, kb, wd, threads))
            bims.update(read_bim(prefix))
    cols = ["snp_a", "snp_b", "chr", "pos_a", "pos_b", "R"]
    ld = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    d = _directional(ld, targets)
    d["r2"] = d.R ** 2
    d = d[d.r2 >= r2]
    if available is not None:
        d = d[d.proxy.isin(available)]

    # allele phase from the panel .bim: R>0 -> A1 pairs with A1
    def allele_info(row):
        ta = bims.get(row.target, (None, None, None, None))
        pa = bims.get(row.proxy, (None, None, None, None))
        same = row.R > 0
        return pd.Series({
            "target_a1": ta[2], "target_a2": ta[3],
            "proxy_a1": pa[2], "proxy_a2": pa[3],
            "proxy_allele_for_target_a1": pa[2] if same else pa[3],
            "proxy_allele_for_target_a2": pa[3] if same else pa[2],
        })
    if len(d):
        d = pd.concat([d.reset_index(drop=True),
                       d.apply(allele_info, axis=1).reset_index(drop=True)], axis=1)
    d["dist_bp"] = (d.pos_proxy - d.pos_target).abs()
    d = d.sort_values(["target", "r2", "dist_bp"], ascending=[True, False, True])
    d = d.groupby("target", sort=False).head(top).reset_index(drop=True)

    # annotate from db (both builds, genes) and cache
    con = sqlite3.connect(str(db))
    ann = []
    for ch in chunked(list(d.proxy.unique()), 900):
        ann.append(pd.read_sql_query(
            f"SELECT rsid AS proxy, pos37 AS proxy_pos37, pos38 AS proxy_pos38, "
            f"ref37 AS proxy_ref37, alt37 AS proxy_alt37, ref38 AS proxy_ref38, "
            f"alt38 AS proxy_alt38, gene AS proxy_gene FROM variants "
            f"WHERE rsid IN ({','.join('?'*len(ch))})", con, params=ch))
    if ann:
        a = pd.concat(ann).drop_duplicates("proxy")
        d = d.merge(a, on="proxy", how="left")
    con.executemany(
        "INSERT INTO ld VALUES(?,?,?,?,?,?,?,?,?)",
        [(r.target, r.proxy, str(r.chr), int(r.pos_target), int(r.pos_proxy),
          float(r.r2), None, panel_build, panel_name) for r in d.itertuples()])
    con.commit()
    con.close()

    # report targets with no proxy
    found = set(d.target)
    d.attrs["no_proxy"] = [t for t in targets if t not in found]
    return d
