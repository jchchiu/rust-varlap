# rust-varlap

rust-varlap is a rewrite of [varlap](https://github.com/bjpop/varlap) from python into rust.

* [Assumptions](#assumptions)
* [Changes](#changes)
* [Installation](#installation)
* [Usage](#usage)
* [Notes](#notes)

## Assumptions
- The input variant files (vcf, csv, tsv) are sorted by variant chromosome and position in ascending order
- Variants are not filtered or interpreted
- The input read files (bam, cram) are sorted by read chromosome and position in ascending order
- BAM, CRAM, FASTA need index files associated (.bai, .crai, .fai, respectively); for VCF, it is needed if gzipped (.tbi)

## Changes
- The algorithm now iterates over reads instead of variants
    - Instead of calling pileup over variants we are now fetching reads for variants within a given range
    - Variants are binned by chromosome, and further binned if variants are some `n` distance apart (gap)
    - From preliminary testing, we have set the default gap size as 100 kb (100,000 bp) as it balances performance between files with sparsely and densely populated variants
    - An optional hyperparameter has can be used to control this gap size between bins; generally speaking, sparsely populated variants throughout chromosomes may perform better if gap size is small (e.g. `5 x mean read length`); conversely, for densely populated variants all throughout the chromosome performance may be better if gap size is large (e.g. 1 Mb or even more)
- Multiple bam/cram parsing is supported but by default they are not merged; there is an optional flag to merge outputs at the end if wanted
- Changed header output field labels:
    - e.g. from [label + ' ' + 'ref avg nm'] to [label + ' ' + 'ref_avg_nm']
    - Easier to split header if necessary (split by ' ' will separate bam label and statistic field)
- CRAM files are now supported
- gzipped input variant files are now supported
- Region mode with bed files/Outliers mode is not currently supported

## Installation

### Option 1: Pre-built binary (Linux x86\_64)

Download the binary directly from the [Releases page](https://github.com/jchchiu/rust-varlap/releases).

```bash
VERSION="v0.1.0-alpha.4"

wget https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}_linux-x86_64
chmod +x rust-varlap_${VERSION}_linux-x86_64

# Optional: move to somewhere on your PATH
mv rust-varlap_${VERSION}_linux-x86_64 ~/.local/bin/rust-varlap

rust-varlap --version
```

---

### Option 2: Apptainer / Singularity `.sif`

Download the binary directly from the [Releases page](https://github.com/jchchiu/rust-varlap/releases).

```bash
VERSION="v0.1.0-alpha.4"

wget https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}.sif
```

**Run:**

```bash
# Directly via the runscript
./rust-varlap_${VERSION}.sif --version

# Or explicitly with apptainer
apptainer run rust-varlap_${VERSION}.sif --version

# Bind directories when data lives outside your home folder
apptainer exec \
  -B /scratch/$USER:/scratch/$USER \
  rust-varlap_${VERSION}.sif \
  rust-varlap --version
```

---

### Option 3: Build from source

Requires Rust ≥ 1.70 and a C toolchain.

**Install:**

```bash
cargo install \
  --git https://github.com/jchchiu/rust-varlap \
  --tag v0.1.0-alpha.4 \
  --locked
```

Places it in `~/.cargo/bin/rust-varlap` (which is on your `PATH` after a standard rustup install).

```bash
rust-varlap --version
```

## Usage

```
rust-varlap --variants <VARIANTS> --reads <READS>... --varclass <VARCLASS> --output <OUTPUT> [OPTIONS]
```

### Required Arguments

| Flag | Description |
|------|-------------|
| `-v, --variants <PATH>` | Path to variants file (`.vcf`, `.csv`, `.tsv`; optionally gzipped with `.gz`) |
| `-r, --reads <PATH>...` | Path to one or more reads files (`.bam`, `.cram`) |
| `-c, --varclass <CLASS>` | Variant class to analyze |
| `-o, --output <PATH>` | Path to output CSV (directory + filename) |

### Optional Arguments

| Flag | Description |
|------|-------------|
| `-f, --fasta <PATH>` | Path to FASTA reference (**required** if any reads file is CRAM) |
| `--sample <NAME>` | Sample identifier |
| `--label <LABEL>...` | Label(s) for reads file(s); defaults to filename if omitted |
| `--gap <BP>` | Bin size gap in base pairs (default: `100,000`) |
| `--merge` | Merge output CSVs when multiple BAMs are provided |

### Examples

Basic usage with a VCF and a single BAM:

```
rust-varlap \
  --variants sample.vcf \
  --reads sample.bam \
  --varclass snv \
  --output results/output.csv
```

Multiple BAMs with custom labels:

```
rust-varlap \
  --variants cohort.tsv \
  --reads tumor.bam normal.bam \
  --label tumor normal \
  --varclass indel \
  --output results/cohort.csv
```

Multiple CRAMs with custom labels, merged output (requires a FASTA reference):

```
rust-varlap \
  --sample patient1 \
  --variants variants.vcf.gz \
  --reads tumor.cram normal.cram  \
  --label tumor normal \
  --fasta reference.fasta \
  --varclass snv \
  --output results/output.csv \
  --merge
```

---

## Notes
- For csv files, multiple alt alleles are not supported (e.g. `A,T`). Please separate them into multiple rows (i.e. row for `A` alt, row for `T` alt)

## To-do
- ~~Fix csv/tsv parsing (figure out best idiomatic way in rust to reduce resources/minimize friction when attempting multithreading/duplicate code)~~
- ~~Add support for dealing with multiple chromosomes~~
- ~~Fix output csv header so it matches varlap exactly (?)~~
    - See changes for further info
- ~~Improve error handling and add logging~~
- Add unit tests/integration tests
- Add region scanning mode (split into subcommands, 'allele' and 'region')
- Implement multithreading
- Deal with paired reads (implement something similar to mosdepth?)
- Add option to sort variants if unsorted
- Make ref/alt input optional (Still output stats that do not need this info)
- Add function to deal with CRAM files that have standardized fasta reference files that can be accessed by bamreader
- Implement automatic outlier detection like original varlap?
    - If outlier detected can also output a command for IGV viewer that automatically opens the variant position of interest?
