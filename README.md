# rust-varlap

[![Build & Release](https://github.com/jchchiu/rust-varlap/actions/workflows/release.yml/badge.svg)](https://github.com/jchchiu/rust-varlap/actions/workflows/release.yml)
[![Docs](https://github.com/jchchiu/rust-varlap/actions/workflows/docs.yml/badge.svg)](https://jchchiu.github.io/rust-varlap/)

rust-varlap is a quality control tool for genetic variants from high-throughput sequencing. It is a rewrite of [varlap](https://github.com/bjpop/varlap) from Python into Rust.

For each variant (or each base of a region), it computes alignment-context metrics from one or more BAM/CRAM files: edit distance, base and mapping quality, alignment length, clipping, strand distribution and allele counts. Results are written as CSV in varlap's output format.

📖 **Documentation:** <https://jchchiu.github.io/rust-varlap/>

* [Installation](#installation)
* [Usage](#usage)
* [Assumptions](#assumptions)
* [Changes from varlap](#changes-from-varlap)
* [Notes](#notes)
* [To-do](#to-do)

## Installation

The quickest way to install is the pre-built binary:

```bash
VERSION="v0.1.0-alpha.5"

wget https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}_linux-x86_64
chmod +x rust-varlap_${VERSION}_linux-x86_64
mv rust-varlap_${VERSION}_linux-x86_64 ~/.local/bin/rust-varlap

rust-varlap --version
```

While rust-varlap is in alpha, always pin a specific `VERSION`. GitHub's `releases/latest` URL doesn't resolve to pre-releases.

There are other ways to install it. The [installation guide](https://jchchiu.github.io/rust-varlap/1_getting_started/installation.html) has the full commands for each one:

| Method | Command / source |
|---|---|
| Pre-built binary | [Releases page](https://github.com/jchchiu/rust-varlap/releases): `linux-x86_64`, `linux-aarch64`, `macos-aarch64`, `macos-x86_64` |
| Docker | `docker run --rm -v $(pwd):/data ghcr.io/jchchiu/rust-varlap:${VERSION} rust-varlap --help` |
| Apptainer / Singularity | `apptainer pull docker://ghcr.io/jchchiu/rust-varlap:${VERSION}`, or download the `.sif` from the Releases page |
| crates.io | `cargo install rust-varlap --version 0.1.0-alpha.5 --locked` |
| Bioconda | Coming soon |
| From source | `cargo install --git https://github.com/jchchiu/rust-varlap --tag ${VERSION} --locked` |

Building from source, including through crates.io, needs Rust, a C toolchain, cmake and the htslib system libraries (zlib, bzip2, xz/lzma, curl, OpenSSL). The [installation guide](https://jchchiu.github.io/rust-varlap/1_getting_started/installation.html#build-from-source) has the packages for each OS.

## Usage

rust-varlap has two subcommands:

```
rust-varlap loci   --variants <VARIANTS> --varclass <snv|indel> --reads <READS>... --output <OUTPUT> [OPTIONS]
rust-varlap region --regions <BED>                               --reads <READS>... --output <OUTPUT> [OPTIONS]
```

* **`loci`** computes metrics for each variant in a VCF/CSV/TSV file.
* **`region`** computes metrics for every base in the regions of a BED file.

### Options shared by both modes

| Flag | Description |
|------|-------------|
| `-r, --reads <PATH>...` | One or more reads files (`.bam`, `.cram`) |
| `-o, --output <PATH>` | Output CSV (directory and filename) |
| `-f, --fasta <PATH>` | FASTA reference (**required** if any reads file is CRAM) |
| `--sample <NAME>` | Sample identifier |
| `--label <LABEL>...` | Label for each reads file (defaults to the filename) |
| `--gap <BP>` | Gap in base pairs that splits variants into bins (default: `100,000`) |
| `--merge` | Merge the per-reads-file CSVs into one |
| `-t, --threads <N>` | Number of reads files to process in parallel (default: `1`) |

### `loci` options

| Flag | Description |
|------|-------------|
| `-v, --variants <PATH>` | Variants file (`.vcf`, `.csv`, `.tsv`, plain or gzipped) |
| `-c, --varclass <CLASS>` | Variant class to analyse: `snv` or `indel` |
| `--alleles <MODE>` | How to treat missing REF/ALT in CSV/TSV input: `strict` (default), `partial` or `none`. VCFs always use `strict`. |

### `region` options

| Flag | Description |
|------|-------------|
| `-R, --regions <PATH>` | BED file of regions |

### Examples

SNVs from a VCF against a single BAM:

```bash
rust-varlap loci \
  -v sample.vcf \
  -r sample.bam \
  -c snv \
  -o results/output.csv
```

Indels from a TSV against tumour and normal BAMs, keeping variants with a missing REF/ALT:

```bash
rust-varlap loci \
  -v cohort.tsv \
  -r tumor.bam normal.bam \
  --label tumor normal \
  -c indel \
  --alleles partial \
  -o results/cohort.csv
```

CRAMs (which need a FASTA reference), merged output and two threads:

```bash
rust-varlap loci \
  --sample patient1 \
  -v variants.vcf.gz \
  -r tumor.cram normal.cram \
  --label tumor normal \
  -f reference.fasta \
  -c snv \
  -o results/output.csv \
  --merge \
  -t 2
```

Every base in a set of BED regions:

```bash
rust-varlap region \
  -R regions.bed \
  -r sample.bam \
  -o results/regions.csv
```

The [CLI reference](https://jchchiu.github.io/rust-varlap/2_usage/cli.html) describes every option.

## Assumptions

- **Sorted variants and regions:** the input variant files (VCF, CSV, TSV) and BED files are sorted by chromosome and then position, in ascending order.
- **No filtering:** variants are not filtered or interpreted.
- **Sorted reads:** the input reads files (BAM, CRAM) are coordinate-sorted.
- **Index files:** BAM, CRAM and FASTA files need their index files (`.bai`, `.crai` and `.fai`). A VCF needs one (`.tbi`) only if it is gzipped.

## Changes from varlap

- **Reads instead of variants:** the algorithm iterates over reads, not variants. Instead of calling pileup for each variant, it fetches the reads for groups of nearby variants. As a result, each read is fetched at most once per bin.
- **Binning:** variants are grouped into bins by chromosome. A new bin starts wherever consecutive variants are more than `--gap` apart.
  - The default gap of 100 kb balances performance between files with sparse variants and files with dense ones.
  - Sparse variants may run faster with a small gap, for example 5 × the mean read length.
  - Dense variants may run faster with a large gap, for example 1 Mb or more.
- **Subcommands:** the CLI has two subcommands, `loci` and `region`. Region mode takes a BED file.
- **Multiple reads files:** these are supported and can be processed in parallel with `--threads`. Their outputs are kept separate unless `--merge` is given.
- **Header labels:** the output header labels have changed, for example from `label + ' ' + 'ref avg nm'` to `label + ' ' + 'ref_avg_nm'`. This makes the header easier to split, because splitting on `' '` separates the reads label from the statistic.
- **New input formats:** rust-varlap accepts CRAM files and gzipped variant files.
- **Optional REF/ALT:** REF/ALT can be optional in CSV/TSV input (`--alleles partial|none`).
- **Truncated reads files:** if a reads file gives 10 or more consecutive read errors, rust-varlap treats it as truncated. It doesn't check for an EOF marker.
- **Not yet supported:** varlap's outlier mode.

## Notes

- **Multiallelic CSV rows:** CSV files don't support multiple alt alleles in one row (e.g. `A,T`). Split them into separate rows, one for the `A` alt and one for the `T` alt.

## To-do

- ~~Fix csv/tsv parsing~~
- ~~Add support for multiple chromosomes~~
- ~~Fix the output csv header~~ (see [Changes from varlap](#changes-from-varlap))
- ~~Improve error handling and add logging~~
- ~~Add region scanning mode (subcommands `loci` and `region`)~~
- ~~Make ref/alt input optional~~
- Add unit and integration tests
- Reduce region mode memory use, and implement multithreading for region mode
- Deal with paired reads (something similar to mosdepth?)
- Add an option to sort unsorted variants
- Support CRAM files whose standard FASTA reference can be fetched by the BAM reader
- Implement automatic outlier detection, like the original varlap
  - If an outlier is detected, also output an IGV command that opens the variant position
