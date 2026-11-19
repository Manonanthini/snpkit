"""Shared helpers: contig naming, fast text opening, chunking."""

import gzip
import io
import shutil
import subprocess

# NCBI RefSeq accession -> bare chromosome name, for dbSNP VCFs
# (GRCh38 = NC_0000NN.1x, GRCh37 = NC_0000NN.1x with a different version).
_REFSEQ = {f"NC_{i:06d}": str(i) for i in range(1, 23)}
_REFSEQ.update({"NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT"})

VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def norm_chrom(c):
    """Return bare chromosome name ('1'..'22','X','Y','MT') or None if unusable."""
    c = str(c)
    if c.startswith("NC_"):
        c = _REFSEQ.get(c.split(".")[0])
        return c
    if c.lower().startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    c = c.upper() if c.upper() in ("X", "Y", "MT") else c
    return c if c in VALID_CHROMS else None


def ucsc_chrom(c):
    """Bare chromosome -> UCSC style ('chr1', 'chrM')."""
    return "chrM" if c == "MT" else f"chr{c}"


def open_text(path):
    """Open (possibly gzipped) text file for fast streaming reads.

    Uses an external pigz/gzip process when available (much faster than
    python's gzip module on dbSNP-sized files).
    """
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        for tool in ("pigz", "gzip"):
            exe = shutil.which(tool)
            if exe:
                proc = subprocess.Popen(
                    [exe, "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 20
                )
                return io.TextIOWrapper(proc.stdout, encoding="utf-8", newline="")
        return gzip.open(path, "rt")
    return open(path, "rt")


def chunked(it, n):
    buf = []
    for x in it:
        buf.append(x)
        if len(buf) >= n:
            yield buf
            buf = []
    if buf:
        yield buf


COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMP)[::-1]
