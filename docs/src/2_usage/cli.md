# CLI

## `loci` mode

```
rust-varlap loci --variants <VARIANTS> --reads <READS> --varclass <VARCLASS> --output <OUTPUT> [OPTIONS]
```

### **Required Arguments**

#### `-v, --variants <PATH>`

Path to a variants file (`.vcf`, `.csv`, or `.tsv`) which can be plain or gzip-compressed (`.gz`).

Gzip compression is detected automatically by inspecting the file header (magic bytes), so the `.gz` extension is not required.

#### `-r, --reads <PATH>...`

Path to one or more reads files (`.bam`, `.cram`).

#### `-c, --varclass <CLASS>`

Variant class to analyze. Can choose between SNVs (`snv`) and INDELS (`indel`)

#### `-o, --output <PATH>`

Path to output CSV (directory + filename).

### Optional Arguments

#### `--sample <NAME>`

Sample identifier.

#### `--label <LABEL>...`

Label(s) for reads file(s). Defaults to the filename if omitted.

**Default:** filename if omitted (e.g. `sample` if `sample.bam`)

#### `--gap <BP>`

Bin size gap in base pairs. 

Generally speaking, sparsely populated variants throughout chromosomes may perform better if gap size is small (e.g. 5 x mean read length); conversely, for densely populated variants all throughout the chromosome performance may be better if gap size is large (e.g. 1 Mb or even more) 

From preliminary testing, we have set the default gap size as 100 kb (100,000 bp) as it balances performance between files with sparsely and densely populated variants

**Default:** `100,000`

#### `--merge`

Merge output CSVs when multiple BAMs are provided. The output will be the same as the original varlap.

**Default:** `false`

#### `-f, --fasta <PATH>`

Path to FASTA reference. **Required** if any reads file is CRAM.

#### `-t, --threads <N THREADS>`

Number of threads used if multithreading is wanted for multiple BAM/CRAM 

**Default:** `1`

#### `--alleles`

How to treat REF/ALT in variant CSV/TSV inputs:

| Mode      | Description |
|---------|-------------|
| `strict`  | Only variants with both REF/ALT will be analyzed |
| `partial` | Variants with empty REF or ALT or both will still be analyzed  |
| `none`    | For all variants REF/ALT are not used at all even if they are in the input |

For example, if a variant contains REF but no ALT, varlap will still output statistics for reads that support the REF.

Note that currently only partial alleles for REF/ALT with a length of 1 (i.e. SNVs) are supported for extra analysis. Otherwise, the variant will be skipped.

For VCFs, the mode is always set to `strict` to ensure properly formatted VCFs.

**Default:** `strict`

### Examples

Basic usage with a VCF and a single BAM analyzing SNVs:

```
rust-varlap loci \
  --variants sample.vcf \
  --reads sample.bam \
  --varclass snv \
  --output results/output.csv
```

Multiple BAMs with custom labels and including empty REF/ALT alleles analyzing INDELs:

```
rust-varlap loci \
  --variants cohort.tsv \
  --reads tumor.bam normal.bam \
  --label tumor normal \
  --varclass indel \
  --alleles partial \
  --output results/cohort.csv
```

Multiple CRAMs with custom labels, merged output (requires a FASTA reference):

```
rust-varlap loci \
  --sample patient1 \
  --variants variants.vcf.gz \
  --reads tumor.cram normal.cram \
  --label tumor normal \
  --fasta reference.fasta \
  --varclass snv \
  --output results/output.csv \
  --merge \
  --threads 2
```

## `region` mode

```
rust-varlap region --regions <REGIONS> --reads <READS> --output <OUTPUT> [OPTIONS]
```

### Required Arguments

#### `-R, --regions <PATH>`

Path to regions file ().

#### `-r, --reads <PATH>...`

Path to one or more reads files (`.bam`, `.cram`).

#### `-o, --output <PATH>`

Path to output CSV (directory + filename).