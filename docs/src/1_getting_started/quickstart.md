# Quick Start

## Loci analysis

To run varlap for variants in a given VCF along with its associated BAM to analyze SNVs:

```bash
rust-varlap loci \
  --variants sample.vcf \
  --reads sample.bam \
  --varclass snv \
  --output results/output.csv
```

For all SNVs in the VCF, varlap will calculate statistics for all reads in the BAM file that overlap each variant, and output these results to `results/output.csv`.

See [Installation](../1_getting_started/installation.md) to get set up, and the [CLI reference](../2_usage/cli.md) for all options.