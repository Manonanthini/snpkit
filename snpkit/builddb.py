"""One-time construction of the SQLite lookup database.

variants : rsid, chr, pos37, pos38, ref37, alt37, ref38, alt38, gene, lifted
genes    : chr, start, end, name, gene_id, biotype, strand   (GRCh38 coordinates)
ld       : cache of proxy searches
meta     : key/value

Both builds' coordinates sit side by side on one row (same idea as PRSedm's
position_hg19/position_hg38 columns), so a later b37<->b38 conversion of a known
variant is a single indexed SELECT instead of a chain remap.
"""

import logging
import sqlite3

from .common import chunked, norm_chrom, open_text, revcomp
from .genes import GeneIndex, parse_gtf_genes
from .lift import Lifter

log = logging.getLogger("snpkit")

SCHEMA = """
CREATE TABLE IF NOT EXISTS variants(
    rsid TEXT, chr TEXT, pos37 INTEGER, pos38 INTEGER,
    ref37 TEXT, alt37 TEXT, ref38 TEXT, alt38 TEXT,
    gene TEXT, lifted INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS genes(
    chr TEXT, start INTEGER, end INTEGER, name TEXT, gene_id TEXT,
    biotype TEXT, strand TEXT);
CREATE TABLE IF NOT EXISTS ld(
    target TEXT, proxy TEXT, chr TEXT, pos_target INTEGER, pos_proxy INTEGER,
    r2 REAL, dprime REAL, build INTEGER, panel TEXT);
CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT);
"""


def connect(path):
    con = sqlite3.connect(str(path))
    con.executescript(SCHEMA)
    return con


def load_rsid_list(path):
    with open(path) as fh:
        return {ln.split()[0] for ln in fh if ln.strip() and not ln.startswith("#")}


def _stage_dbsnp(con, vcf, build, rsids=None, snv_only=False, batch=200_000):
    """Stream a dbSNP VCF into staging table s{build}(rsid,chr,pos,ref,alt)."""
    t = f"s{build}"
    con.execute(f"DROP TABLE IF EXISTS {t}")
    con.execute(f"CREATE TABLE {t}(rsid TEXT, chr TEXT, pos INTEGER, ref TEXT, alt TEXT)")

    def rows():
        with open_text(vcf) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 5)
                rid = f[2]
                if rid == "." or (rsids is not None and rid not in rsids):
                    continue
                c = norm_chrom(f[0])
                if c is None:
                    continue
                ref = f[3]
                for alt in f[4].split(","):
                    if snv_only and (len(ref) != 1 or len(alt) != 1):
                        continue
                    yield (rid, c, int(f[1]), ref, alt)

    n = 0
    for b in chunked(rows(), batch):
        con.executemany(f"INSERT INTO {t} VALUES(?,?,?,?,?)", b)
        n += len(b)
        log.info("b%s: %s variants staged", build, f"{n:,}")
    con.execute(f"CREATE INDEX {t}_k ON {t}(rsid, alt)")
    con.commit()
    return n


def _merge_builds(con, have37, have38):
    con.execute("DELETE FROM variants")
    if have38 and have37:
        con.execute(
            """INSERT INTO variants(rsid,chr,pos37,pos38,ref37,alt37,ref38,alt38)
               SELECT a.rsid,a.chr,b.pos,a.pos,b.ref,b.alt,a.ref,a.alt
               FROM s38 a LEFT JOIN s37 b ON a.rsid=b.rsid AND a.alt=b.alt"""
        )
        con.execute(
            """INSERT INTO variants(rsid,chr,pos37,ref37,alt37)
               SELECT b.rsid,b.chr,b.pos,b.ref,b.alt FROM s37 b
               WHERE NOT EXISTS(SELECT 1 FROM s38 a WHERE a.rsid=b.rsid AND a.alt=b.alt)"""
        )
    elif have38:
        con.execute(
            """INSERT INTO variants(rsid,chr,pos38,ref38,alt38)
               SELECT rsid,chr,pos,ref,alt FROM s38"""
        )
    else:
        con.execute(
            """INSERT INTO variants(rsid,chr,pos37,ref37,alt37)
               SELECT rsid,chr,pos,ref,alt FROM s37"""
        )
    con.commit()


def _fill_missing_build(con, chain, src, dst):
    """Chain-lift rows that only have src-build coordinates into dst columns."""
    lifter = Lifter(chain)
    cur = con.execute(
        f"""SELECT rowid,chr,pos{src},ref{src},alt{src} FROM variants
            WHERE pos{src} IS NOT NULL AND pos{dst} IS NULL"""
    )
    upd = []
    for rid, c, p, ref, alt in cur.fetchall():
        hit = lifter.lift(c, p)
        if hit is None:
            continue
        tc, tp, strand = hit
        if tc != c:  # moved chromosome: don't trust
            continue
        if strand == "-":
            ref, alt = revcomp(ref), revcomp(alt)
        upd.append((tp, ref, alt, rid))
    con.executemany(
        f"UPDATE variants SET pos{dst}=?, ref{dst}=?, alt{dst}=?, lifted=1 WHERE rowid=?",
        upd,
    )
    con.commit()
    log.info("chain-lifted %s variants %s->%s", f"{len(upd):,}", src, dst)


def _load_genes_and_annotate(con, gtf, biotypes):
    con.execute("DELETE FROM genes")
    rows = list(parse_gtf_genes(gtf, biotypes))
    con.executemany("INSERT INTO genes VALUES(?,?,?,?,?,?,?)", rows)
    con.commit()
    gi = GeneIndex(rows)
    cur = con.execute("SELECT rowid,chr,pos38 FROM variants WHERE pos38 IS NOT NULL")
    for batch in chunked(cur, 200_000):
        con.executemany(
            "UPDATE variants SET gene=? WHERE rowid=?",
            [(gi.genes_at(c, p), rid) for rid, c, p in batch],
        )
        con.commit()
    log.info("annotated genes using %s gene models", f"{len(rows):,}")


def build_db(
    out,
    dbsnp38=None,
    dbsnp37=None,
    chain_38to37=None,
    chain_37to38=None,
    gtf=None,
    rsid_list=None,
    snv_only=False,
    biotypes=None,
):
    if not (dbsnp38 or dbsnp37):
        raise ValueError("need at least one dbSNP VCF (--dbsnp38 and/or --dbsnp37)")
    con = connect(out)
    rsids = load_rsid_list(rsid_list) if rsid_list else None
    if dbsnp38:
        _stage_dbsnp(con, dbsnp38, 38, rsids, snv_only)
    if dbsnp37:
        _stage_dbsnp(con, dbsnp37, 37, rsids, snv_only)
    _merge_builds(con, bool(dbsnp38), bool(dbsnp37))
    for t in ("s37", "s38"):
        con.execute(f"DROP TABLE IF EXISTS {t}")
    if chain_38to37:
        _fill_missing_build(con, chain_38to37, 38, 37)
    if chain_37to38:
        _fill_missing_build(con, chain_37to38, 37, 38)
    if gtf:
        _load_genes_and_annotate(con, gtf, set(biotypes) if biotypes else None)
    for sql in (
        "CREATE INDEX IF NOT EXISTS v_rsid ON variants(rsid)",
        "CREATE INDEX IF NOT EXISTS v_p37 ON variants(chr,pos37)",
        "CREATE INDEX IF NOT EXISTS v_p38 ON variants(chr,pos38)",
        "CREATE INDEX IF NOT EXISTS g_chr ON genes(chr,start)",
        "CREATE INDEX IF NOT EXISTS ld_t ON ld(target,build)",
    ):
        con.execute(sql)
    con.execute("INSERT OR REPLACE INTO meta VALUES('dbsnp','b156 (user-supplied VCFs)')")
    con.commit()
    n = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    con.execute("VACUUM")
    con.close()
    log.info("db ready: %s (%s variants)", out, f"{n:,}")
    return n
