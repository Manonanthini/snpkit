"""Generate tiny synthetic dbSNP VCFs, chain files, GTF, PLINK bim and a fake plink."""

import gzip
import os
import stat
import sys

d = sys.argv[1]
os.makedirs(d, exist_ok=True)


def w(name, text, gz=False):
    p = os.path.join(d, name)
    (gzip.open if gz else open)(p, "wt").write(text)
    return p


hdr = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
v38 = [
    ("NC_000001.11", 101000, "rs1", "A", "G"),
    ("NC_000001.11", 102000, "rsP1", "C", "T"),
    ("NC_000001.11", 103000, "rsP2", "G", "A"),
    ("NC_000001.11", 201000, "rs2", "C", "T"),
    ("NC_000001.11", 301000, "rs3", "G", "T,C"),
    ("NC_000001.11", 5000000, "rs4", "A", "C"),
]
v37 = [
    ("NC_000001.10", 100000, "rs1", "A", "G"),
    ("NC_000001.10", 101000, "rsP1", "C", "T"),
    ("NC_000001.10", 102000, "rsP2", "G", "A"),
    ("NC_000001.10", 200000, "rs2", "C", "T"),
    ("NC_000001.10", 300000, "rs3", "G", "T,C"),
    ("NC_000001.10", 400000, "rs5", "T", "C"),
]
fmt = lambda rows: "".join(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\tRS=1\n" for c, p, i, r, a in rows)
w("dbsnp38.vcf.gz", hdr + fmt(v38), gz=True)
w("dbsnp37.vcf", hdr + fmt(v37))


# chain: hg38 pos = hg19 pos + 1000 on chr1 (blocks cover 0..1,000,000 of hg19)
def chain(t, q, tstart, qstart):
    return (f"chain 1000 chr1 3000000 + {tstart} {tstart+1000000} chr1 3000000 + "
            f"{qstart} {qstart+1000000} 1\n1000000\n\n")


w("hg38ToHg19.chain.gz", chain(38, 37, 1000, 0), gz=True)
w("hg19ToHg38.chain", chain(37, 38, 0, 1000))

w("genes.gtf",
  '1\tHAVANA\tgene\t50000\t250000\t.\t+\t.\tgene_id "G0"; gene_type "lncRNA"; gene_name "LONGGENE";\n'
  '1\tHAVANA\tgene\t100500\t101500\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "GENEA";\n'
  '1\tHAVANA\tgene\t200800\t202000\t.\t-\t.\tgene_id "G2"; gene_type "protein_coding"; gene_name "GENEB";\n'
  '1\tHAVANA\texon\t100500\t101500\t.\t+\t.\tgene_id "G1"; gene_name "GENEA";\n')

# panel .bim (variant IDs = rsIDs): chr rsid cm pos a1 a2
w("panel.bim",
  "1\trs1\t0\t101000\tG\tA\n1\trsP1\t0\t102000\tT\tC\n1\trsP2\t0\t103000\tA\tG\n"
  "1\trs2\t0\t201000\tT\tC\n")

# fake plink: ignores real LD, writes a canned .ld file
fp = os.path.join(d, "fake_plink")
open(fp, "w").write(f"""#!{sys.executable}
import sys
a = sys.argv
out = a[a.index("--out") + 1]
H = "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R\\n"
rows = [
 "1 101000 rs1 1 102000 rsP1 0.95",
 "1 101000 rs1 1 103000 rsP2 -0.90",
 "1 201000 rs2 1 102000 rsP1 0.40",
]
open(out + ".ld", "w").write(H + "\\n".join(rows) + "\\n")
""")
os.chmod(fp, os.stat(fp).st_mode | stat.S_IEXEC)

w("in38.txt",
  "rsid chr pos ref alt\n"
  "rs1 1 101000 A G\n"          # exact
  "rs2 1 201000 T C\n"          # swapped
  ". 1 301000 G C\n"            # pos match, no rsid, multiallelic -> best alt
  "rsX 1 700000 A C\n"          # not in db -> chain
  "rsY 1 9000000 A C\n")        # not mappable
w("in37_rsids.txt", "rs1\nrs5\nrsNOPE\n")
w("targets.txt", "rs1\nrs2\n")
