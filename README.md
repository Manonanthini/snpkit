# snpkit

snpkit is a command-line tool for variant work across human genome builds. It operates on plain-text files and a local SQLite database, and provides:

- **Liftover** between GRCh37 (b37) and GRCh38 (b38) using chain files, with precomputed positions for both builds stored in the database.
- **Annotation** of variants with rsID, chromosome, positions in both builds, REF/ALT alleles and gene.
- **LD proxy search** using PLINK and a user-supplied reference panel.

All look-ups are local; no network requests are made at run time. The storage design (both builds' positions on a single database record) was inspired by [PRSedm](https://github.com/sethsh7/PRSedm); no PRSedm code is included.

## Installation

Requirements: Python 3.9 or later. PLINK 1.9 is required only for proxy search.

```bash
git clone <repository-url>
cd snpkit
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
snpkit --help
```

## Commands

| Command | Function |
| --- | --- |
| `snpkit build-db` | Build the look-up database from dbSNP VCF files, chain files and a GTF file (performed once) |
| `snpkit annotate` | Add b37/b38 positions, alleles and gene to a variant file |
| `snpkit liftover` | Convert chromosome and position with a chain file (no database required) |
| `snpkit proxy` | Identify LD proxies for target rsIDs |

## Quick start

```bash
# Build the database (restrict to the rsIDs of interest to limit its size)
snpkit build-db --out variants.db \
    --dbsnp38 GCF_000001405.40.gz --dbsnp37 GCF_000001405.25.gz \
    --chain-37to38 hg19ToHg38.over.chain.gz --chain-38to37 hg38ToHg19.over.chain.gz \
    --gtf gencode.annotation.gtf.gz --rsid-list rsids.txt

# Annotate variants given in build 37
snpkit annotate sumstats.txt --db variants.db --build 37 --out sumstats.annotated.tsv

# Identify LD proxies
snpkit proxy targets.txt --db variants.db --panel /path/to/panel_chr{chr} --out proxies.tsv
```

## Documentation

The full guide, including installation, verification with synthetic test data, reference-data sources, input and output formats, and troubleshooting, is in [docs/USER_GUIDE.md](docs/USER_GUIDE.md).

## Status

The software has been tested on small synthetic data sets and has not yet been validated against complete dbSNP b156 data or real PLINK output. A small subset should be tested and compared manually with dbSNP before large-scale use. dbSNP, UCSC chain files and GENCODE annotations are distributed under their own terms.
