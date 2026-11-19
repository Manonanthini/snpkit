# snpkit User Guide

## 1. Introduction

snpkit converts genetic variants between the GRCh37 (b37) and GRCh38 (b38) reference assemblies, annotates them with rsID, chromosome, positions in both builds, REF/ALT alleles and gene, and identifies linkage disequilibrium (LD) proxy variants. It reads and writes plain-text files and performs all look-ups against a single local SQLite database; no network or web-service requests are made at run time.

| Command | Function | Requirements |
| --- | --- | --- |
| `snpkit build-db` | Builds the look-up database (performed once) | dbSNP VCF file(s), chain file(s), GTF file |
| `snpkit annotate` | Adds b37 and b38 positions, alleles and gene to a variant file | The database |
| `snpkit liftover` | Converts chromosome and position only | Chain file |
| `snpkit proxy` | Identifies LD proxies for target rsIDs | The database, PLINK 1.9, a reference panel |

Both builds' coordinates are stored on the same database row, so converting a variant that is present in the database requires a single indexed look-up. A chain file is used only for variants absent from the database. This storage design was inspired by PRSedm, which keeps hg19 and hg38 positions side by side in a local SQLite database.

## 2. Requirements

| Item | Required for |
| --- | --- |
| Python 3.9 or later | All commands |
| pandas (2.0 or later), numpy, pyliftover | All commands (installed automatically by pip) |
| PLINK 1.9 | `snpkit proxy` only |
| pigz | Optional; accelerates reading of compressed VCF files (gzip is used if absent) |

## 3. Installation

1. Obtain the source code.

    ```bash
    git clone <repository-url>
    cd snpkit
    ```

2. Create and activate an isolated Python environment.

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```

3. Install the package in editable mode. Editable mode allows later updates to take effect without reinstallation.

    ```bash
    pip install -e .
    ```

4. Verify the installation. The output should list the commands `build-db`, `annotate`, `liftover` and `proxy`.

    ```bash
    snpkit --help
    ```

5. For proxy search only, confirm that PLINK 1.9 is available. If it is not on the executable path, supply its location with the `--plink` option.

    ```bash
    plink --version
    ```

If the message `snpkit: command not found` appears, the environment is not active or step 3 was not completed.

## 4. Verification with test data

The test data set is synthetic and requires no downloads. Run the commands from the repository root (the directory containing `pyproject.toml`).

```bash
python tests/make_toy.py toy_data
cd toy_data
snpkit build-db --out v.db --dbsnp38 dbsnp38.vcf.gz --dbsnp37 dbsnp37.vcf \
    --chain-37to38 hg19ToHg38.chain --chain-38to37 hg38ToHg19.chain.gz --gtf genes.gtf
snpkit annotate in38.txt --db v.db --build 38 --chain hg38ToHg19.chain.gz --out a.tsv
head a.tsv
```

The build log ends with `db ready: v.db (8 variants)`, and the annotation step reports `annotated 5 variants (1 not found)`. The expected results for the five input variants are as follows.

| Input | b37 / b38 position | Gene | Source | Allele status | Interpretation |
| --- | --- | --- | --- | --- | --- |
| `rs1` chr1:101000 A/G | 100000 / 101000 | GENEA, LONGGENE | `db_rsid` | `exact` | Matched by rsID; alleles agree |
| `rs2` chr1:201000 T/C | 200000 / 201000 | GENEB, LONGGENE | `db_rsid` | `swapped` | Input REF/ALT are reversed relative to dbSNP (C/T) |
| `.` chr1:301000 G/C | 300000 / 301000 | none | `db_pos` | `exact` | No rsID supplied; matched by position and rsID `rs3` recovered |
| `rsX` chr1:700000 A/C | 699000 / 700000 | none | `chain` | `chain_no_allele_check` | Absent from the database; position converted by chain file only |
| `rsY` chr1:9000000 A/C | none | none | `none` | `unmapped` | Absent from the database and outside the chain blocks |

The message "1 not found" in the log refers to `rsY`. Variant `rsX` was resolved by the chain-file fallback.

To test proxy search, run the following. The test uses a stand-in for PLINK that returns fixed LD values, so it verifies only the program logic. The expected result is two proxies for `rs1` (`rsP1`, r² = 0.9025; `rsP2`, r² = 0.81) and a log message that no proxy was found for `rs2`.

```bash
snpkit proxy targets.txt --db v.db --panel panel --plink ./fake_plink --r2 0.8 --out p.tsv
cat p.tsv
```

## 5. Reference data

Real analyses require dbSNP build 156, a chain file for each direction of conversion, a GRCh38 GTF file and, for proxy search, a PLINK reference panel. Existing institutional copies should be used where available.

| File | Source | Purpose |
| --- | --- | --- |
| dbSNP b156 GRCh38 VCF (`GCF_000001405.40.gz`) | [NCBI dbSNP releases](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/) | rsID, b38 position, REF/ALT |
| dbSNP b156 GRCh37 VCF (`GCF_000001405.25.gz`) | Same location | b37 position, REF/ALT |
| `hg19ToHg38.over.chain.gz` | [UCSC hg19 liftOver](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) | b37 to b38 conversion |
| `hg38ToHg19.over.chain.gz` | [UCSC hg38 liftOver](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) | b38 to b37 conversion |
| GENCODE GTF (GRCh38) | [GENCODE human](https://www.gencodegenes.org/human/) | Gene annotation |
| PLINK panel with rsIDs as variant identifiers | User-provided (for example 1000 Genomes) | LD proxy search |

The dbSNP files are several tens of gigabytes each. NCBI contig names (for example `NC_000001.11`) and both `chr1` and `1` styles are handled. The `.bim` file of the panel must use rsIDs as variant identifiers; identifiers of the form `chr:pos:ref:alt` must be converted to rsIDs beforehand, otherwise proxies cannot be matched to the database.

## 6. Building the database

The database is built once. A complete dbSNP database contains approximately one billion variants and would exceed 100 GB; the build should therefore be restricted to a list of required rsIDs. On shared computing systems the build should be submitted as a scheduled job rather than run interactively.

1. Prepare a file containing one rsID per line, for example from the first column of a summary-statistics file with a header row.

    ```bash
    cut -f1 sumstats.txt | tail -n +2 | sort -u > rsids.txt
    ```

2. Build the database.

    ```bash
    snpkit build-db --out variants.db \
        --dbsnp38 GCF_000001405.40.gz --dbsnp37 GCF_000001405.25.gz \
        --chain-37to38 hg19ToHg38.over.chain.gz --chain-38to37 hg38ToHg19.over.chain.gz \
        --gtf gencode.annotation.gtf.gz --biotypes protein_coding lncRNA \
        --rsid-list rsids.txt --snv-only
    ```

3. Confirm that the log ends with `db ready: variants.db (N variants)`.

| Option | Description |
| --- | --- |
| `--dbsnp38`, `--dbsnp37` | One or both dbSNP VCF files; if only one is given, the other build is filled by chain conversion |
| `--chain-37to38`, `--chain-38to37` | Fill coordinates missing from the other build's VCF (recorded as `lifted=1`) |
| `--gtf`, `--biotypes` | Gene annotation; `--biotypes` retains only the listed gene types. The GTF must use GRCh38 |
| `--rsid-list` | Load only the listed rsIDs (strongly recommended) |
| `--snv-only` | Exclude insertions and deletions |

Before use, a small number of variants (five to ten) should be compared manually against dbSNP or the UCSC Genome Browser to confirm positions and alleles.

## 7. Annotating variants

The `annotate` command matches variants by rsID and then by chromosome and position, and adds both builds' positions, alleles and genes. The `--build` option states the build of the input file; the output always contains both builds.

```bash
snpkit annotate sumstats.txt --db variants.db --build 37 \
    --chain hg19ToHg38.over.chain.gz --out sumstats.annotated.tsv
```

The `--chain` option is optional and is used only for variants not present in the database. Additional columns in the input file are preserved.

### 7.1 Input formats

Two input formats are accepted.

1. A tab- or space-separated file with a header row. Columns are identified by the following names (case-insensitive).

    | Field | Recognised header names |
    | --- | --- |
    | rsID | `rsid`, `snp`, `id`, `rs`, `variant_id`, `markername`, `snpid` |
    | Chromosome | `chr`, `chrom`, `chromosome`, `#CHROM` (`chr1` or `1`) |
    | Position | `pos`, `bp`, `position`, `base_pair_location`, `chr_position` |
    | REF / other allele | `ref`, `ref_allele`, `other_allele`, `non_effect_allele` |
    | ALT / effect allele | `alt`, `alt_allele`, `effect_allele`, `a1` |

    Where column names differ, they can be specified explicitly:

    ```bash
    snpkit annotate in.txt --db variants.db --build 38 --out out.tsv \
        --rsid-col SNP --chr-col CHROM --pos-col BP --ref-col A2 --alt-col A1
    ```

2. A single-column file without a header, containing rsIDs (for example `rs1`) or strings of the form `chr:pos` or `chr:pos:ref:alt`.

When alleles are supplied, the best-matching database record is selected. When they are not, all ALT alleles of a multi-allelic site are returned.

### 7.2 Output columns

The following columns are appended after the input columns.

| Column | Description |
| --- | --- |
| `out_rsid`, `out_chr` | rsID and chromosome from the database |
| `out_pos37`, `out_pos38` | Position in each build |
| `out_ref37`, `out_alt37`, `out_ref38`, `out_alt38` | Alleles in each build |
| `out_gene` | Overlapping gene(s); empty indicates an intergenic position |
| `out_source` | `db_rsid`, `db_pos`, `chain`, `chain_minus_strand` or `none` |
| `out_allele_status` | `exact`, `swapped`, `strand_flip`, `strand_flip_swapped`, `mismatch`, `no_alleles`, `not_found`, `unmapped` or `chain_no_allele_check` |

The column `out_allele_status` must be inspected before alleles are used in downstream analyses. A status of `swapped` indicates that the input effect allele corresponds to the database REF allele, so the direction of effect must be reversed. Every `mismatch` should be reviewed. Records with source `chain` have a converted position only and contain no rsID, alleles or gene.

## 8. Liftover of positions

The `liftover` command converts chromosome and position using a chain file and does not require the database. It does not verify alleles or rsIDs; the `annotate` command should be used when these are required.

```bash
snpkit liftover in38.txt --chain hg38ToHg19.chain.gz --out lifted.tsv
```

The output retains the input columns and adds `lift_chr`, `lift_pos`, `lift_strand` and `lift_status` (`ok` or `unmapped`). With the test data, the first four records are converted successfully (for example, `rs1` at position 101000 becomes 100000, reflecting the 1,000 bp offset of the test chain file), and `rsY` at position 9,000,000 is reported as `unmapped` because it lies outside the chain blocks.

## 9. Proxy search

The `proxy` command identifies variants in LD with the target rsIDs by running PLINK on a reference panel, and annotates each proxy from the database. LD is calculated only for the target variants; no genome-wide LD table is precomputed.

```bash
snpkit proxy targets.txt --db variants.db \
    --panel /path/to/panel_chr{chr} --panel-name PANEL --panel-build 37 \
    --plink plink --r2 0.8 --kb 500 --top 5 \
    --available available_rsids.txt --out proxies.tsv
```

The file `targets.txt` contains one rsID per line (first column). The ancestry of the reference panel should match that of the study population.

| Option | Description |
| --- | --- |
| `--panel` | PLINK file prefix without `.bed`; `{chr}` may be used for per-chromosome file sets |
| `--panel-build` | Build of the panel positions (37 or 38), recorded with the results |
| `--r2` | Minimum r² (default 0.8) |
| `--kb` | Window around each target, in kilobases (default 500) |
| `--top` | Number of proxies retained per target (default 5), ranked by r² and then by distance |
| `--available` | File of rsIDs to which proxies must belong, for example variants present in the target data set |
| `--threads` | Number of PLINK threads |

Output columns:

| Column | Description |
| --- | --- |
| `target`, `proxy` | The two rsIDs |
| `R`, `r2`, `dist_bp` | Signed correlation, its square, and distance in base pairs |
| `proxy_pos37`, `proxy_pos38`, `proxy_ref*`, `proxy_alt*`, `proxy_gene` | Proxy positions, alleles and gene from the database |
| `proxy_allele_for_target_a1`, `proxy_allele_for_target_a2` | The proxy allele that is co-inherited with each target allele |

The last two columns allow an effect allele to be transferred from a target to its proxy. A positive R indicates that allele a1 of the target pairs with allele a1 of the proxy; a negative R indicates that it pairs with allele a2.

In the test example, `rs1` (alleles G/A) has two proxies: `rsP1` with R = 0.95 (r² = 0.9025), where target allele G corresponds to proxy allele T, and `rsP2` with R = −0.90 (r² = 0.81), where target allele G corresponds to proxy allele G. No proxy is found for `rs2` at r² ≥ 0.8, and this is reported in the log.

Proxies with A/T or C/G alleles are not removed; strand ambiguity must be assessed by the user.

## 10. Troubleshooting

| Problem | Cause and resolution |
| --- | --- |
| `snpkit: command not found` | The environment is not active or installation was not completed. Activate the environment and run `pip install -e .` in the repository directory |
| `FileNotFoundError` for a test input file | The test files have not been generated, or the command was run from the wrong directory. Run `python tests/make_toy.py toy_data` and change into `toy_data` |
| `'plink' not found on PATH` | Install PLINK 1.9 or supply its full path with `--plink` |
| All records report `not_found` | The `--build` value does not match the input, or the rsIDs were excluded by `--rsid-list` when the database was built |
| Many `swapped` or `mismatch` records | Allele coding differs from dbSNP; inspect `out_allele_status` before using alleles |
| No proxies for any target | Panel variant identifiers are not rsIDs, chromosome naming differs, or `--r2` is too strict |
| Database build is slow or the database is very large | Use `--rsid-list` and `--snv-only`, install pigz, and run the build as a scheduled job |

## 11. Updating

Because the installation is editable, updating requires only that the latest source be retrieved from within the repository directory.

```bash
git pull
```

## 12. Validation status and limitations

The software has been tested on small synthetic data sets, in a speed test of 200,000 look-ups against a database of one million variants (approximately 3 seconds), and with the verification procedure in Section 4. It has not yet been validated against complete dbSNP b156 data or real PLINK output. A small subset should be tested, and a number of variants compared manually with dbSNP, before large-scale use.

Data files such as databases, VCF files and genotype or summary-statistic files should not be committed to the source repository. dbSNP, UCSC chain files and GENCODE annotations are distributed under their own terms.

## 13. Acknowledgements

The use of a local SQLite database with positions for both genome builds stored on a single record was inspired by [PRSedm](https://github.com/sethsh7/PRSedm). No PRSedm code is included in snpkit.
