"""Gene models from a GTF and fast position->gene lookup (no per-SNP API calls)."""

import bisect
import re

import numpy as np

from .common import norm_chrom, open_text

_ATTR = re.compile(r'(\w+) "([^"]*)"')


def parse_gtf_genes(gtf_path, biotypes=None):
    """Yield (chr, start, end, gene_name, gene_id, biotype, strand) for 'gene' rows."""
    with open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            c = norm_chrom(f[0])
            if c is None:
                continue
            attrs = dict(_ATTR.findall(f[8]))
            bt = attrs.get("gene_type") or attrs.get("gene_biotype") or ""
            if biotypes and bt not in biotypes:
                continue
            gid = attrs.get("gene_id", "")
            name = attrs.get("gene_name", gid)
            yield (c, int(f[3]), int(f[4]), name, gid, bt, f[6])


class GeneIndex:
    """Per-chromosome sorted intervals; overlap query by bisect + prefix-max-end."""

    def __init__(self, rows):
        by_chr = {}
        for c, s, e, name, *_ in rows:
            by_chr.setdefault(c, []).append((s, e, name))
        self.idx = {}
        for c, lst in by_chr.items():
            lst.sort()
            starts = [x[0] for x in lst]
            ends = np.array([x[1] for x in lst])
            pmax = np.maximum.accumulate(ends)
            self.idx[c] = (starts, ends, pmax, [x[2] for x in lst])

    def genes_at(self, chrom, pos):
        """Names of genes overlapping a 1-based position ('' if intergenic)."""
        if chrom not in self.idx:
            return ""
        starts, ends, pmax, names = self.idx[chrom]
        i = bisect.bisect_right(starts, pos) - 1
        out = []
        while i >= 0 and pmax[i] >= pos:
            if ends[i] >= pos:
                out.append(names[i])
            i -= 1
        return ",".join(sorted(set(out)))
