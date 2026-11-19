"""Command line interface: snpkit build-db | annotate | liftover | proxy."""

import argparse
import logging
import sys

import pandas as pd


def _colmap(args):
    return {k: getattr(args, f"{k}_col") for k in ("rsid", "chr", "pos", "ref", "alt")
            if getattr(args, f"{k}_col", None)}


def _add_cols(p):
    for k in ("rsid", "chr", "pos", "ref", "alt"):
        p.add_argument(f"--{k}-col", help=f"input column name for {k} (auto-detected by default)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="snpkit", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build-db", help="build the SQLite db from dbSNP VCF(s), chain files, GTF")
    b.add_argument("--out", required=True)
    b.add_argument("--dbsnp38", help="dbSNP b156 GRCh38 VCF (GCF_000001405.40.gz)")
    b.add_argument("--dbsnp37", help="dbSNP b156 GRCh37 VCF (GCF_000001405.25.gz)")
    b.add_argument("--chain-38to37", help="hg38ToHg19.over.chain(.gz): fills missing b37 coords")
    b.add_argument("--chain-37to38", help="hg19ToHg38.over.chain(.gz): fills missing b38 coords")
    b.add_argument("--gtf", help="GENCODE/Ensembl GTF (GRCh38) for gene annotation")
    b.add_argument("--biotypes", nargs="*", help="keep only these gene types, e.g. protein_coding lncRNA")
    b.add_argument("--rsid-list", help="only load these rsIDs (keeps the db small)")
    b.add_argument("--snv-only", action="store_true")

    a = sub.add_parser("annotate", help="b37<->b38 coords, ref/alt, gene for a variant file")
    a.add_argument("input")
    a.add_argument("--db", required=True)
    a.add_argument("--build", type=int, choices=[37, 38], required=True, help="build of the input")
    a.add_argument("--chain", help="chain file (input build -> other) for variants missing from db")
    a.add_argument("--out", required=True)
    _add_cols(a)

    l = sub.add_parser("liftover", help="pure chain-file liftover of chr/pos (no db needed)")
    l.add_argument("input")
    l.add_argument("--chain", required=True)
    l.add_argument("--out", required=True)
    _add_cols(l)

    p = sub.add_parser("proxy", help="LD proxies for target rsIDs via PLINK on a reference panel")
    p.add_argument("targets", help="text file, first column = rsID")
    p.add_argument("--db", required=True)
    p.add_argument("--panel", required=True, help="PLINK bfile prefix; may contain {chr}")
    p.add_argument("--panel-name", default="panel")
    p.add_argument("--panel-build", type=int, default=38, choices=[37, 38])
    p.add_argument("--r2", type=float, default=0.8)
    p.add_argument("--kb", type=int, default=500)
    p.add_argument("--top", type=int, default=5)
    p.add_argument("--available", help="text file of rsIDs proxies must belong to")
    p.add_argument("--plink", default="plink")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out", required=True)

    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", stream=sys.stderr)

    if args.cmd == "build-db":
        from .builddb import build_db
        build_db(args.out, args.dbsnp38, args.dbsnp37, args.chain_38to37, args.chain_37to38,
                 args.gtf, args.rsid_list, args.snv_only, args.biotypes)
    elif args.cmd == "annotate":
        from .annotate import annotate
        annotate(args.db, args.input, args.build, _colmap(args), args.chain).to_csv(
            args.out, sep="\t", index=False)
    elif args.cmd == "liftover":
        from .annotate import read_variants
        from .lift import Lifter
        raw, std = read_variants(args.input, _colmap(args))
        lf = Lifter(args.chain)
        res = []
        for c, p in zip(std.chr, std.pos):
            h = lf.lift(c, int(p)) if c and not pd.isna(p) else None
            res.append(h if h else (None, None, None))
        raw = raw.reset_index(drop=True)
        raw["lift_chr"], raw["lift_pos"], raw["lift_strand"] = zip(*res) if res else ([], [], [])
        raw["lift_pos"] = pd.to_numeric(raw["lift_pos"], errors="coerce").astype("Int64")
        raw["lift_status"] = ["ok" if pd.notna(x) else "unmapped" for x in raw["lift_chr"]]
        raw.to_csv(args.out, sep="\t", index=False)
    elif args.cmd == "proxy":
        from .common import chunked  # noqa: F401
        from .proxy import find_proxies
        with open(args.targets) as fh:
            targets = [ln.split()[0] for ln in fh if ln.strip() and not ln.startswith("#")]
        avail = None
        if args.available:
            with open(args.available) as fh:
                avail = {ln.split()[0] for ln in fh if ln.strip()}
        d = find_proxies(args.db, targets, args.panel, args.r2, args.kb, avail, args.top,
                         args.plink, args.panel_name, args.panel_build, args.threads)
        d.to_csv(args.out, sep="\t", index=False)
        if d.attrs.get("no_proxy"):
            logging.info("no proxy found for %d targets: %s", len(d.attrs["no_proxy"]),
                         ",".join(d.attrs["no_proxy"][:20]))


if __name__ == "__main__":
    main()
