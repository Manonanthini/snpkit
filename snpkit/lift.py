"""Chain-file liftover (UCSC chain, text) with bulk helpers.

Used for (a) filling the second build's coordinates when building the db and
(b) positions that are not in the db.  Coordinates are 1-based in/out.
"""

from pyliftover import LiftOver

from .common import norm_chrom, ucsc_chrom


class Lifter:
    """Wraps pyliftover; loads the chain file once, O(log n) per position."""

    def __init__(self, chain_path):
        self.lo = LiftOver(str(chain_path))

    def lift(self, chrom, pos):
        """Return (chrom, pos, strand) or None. chrom bare, pos 1-based."""
        c = norm_chrom(chrom)
        if c is None:
            return None
        hits = self.lo.convert_coordinate(ucsc_chrom(c), int(pos) - 1)
        if not hits:
            return None
        # Take the highest-scoring hit; ignore ambiguous multi-mapping ties
        hit = max(hits, key=lambda h: h[3])
        tc = norm_chrom(hit[0])
        if tc is None:
            return None
        return tc, hit[1] + 1, hit[2]

    def lift_many(self, rows):
        """rows: iterable of (chrom, pos). Returns list of (chrom,pos,strand)|None."""
        return [self.lift(c, p) for c, p in rows]
