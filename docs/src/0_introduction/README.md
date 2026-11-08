# Introduction

rust-varlap is a quality control tool for genetic variants arising from high-throughput DNA sequencing, where the variants have been called by aligning DNA sequencing reads to a reference genome. It is a rewrite/reimplementation of [varlap](https://github.com/bjpop/varlap) from Python into Rust.

Large-scale sequencing studies, particularly in cancer genomics, generate thousands of candidate variants, many arising from problematic alignment contexts such as low depth of coverage, sample contamination, or proximity to other variants. Manually inspecting these in tools such as IGV is effective for small variant sets but impractical at cohort scale.

## What it does

Varlap is an open-source tool that automates this process, computing per-variant alignment-context metrics from one or more BAM files, including:

- edit distance
- base and mapping quality
- alignment length
- clipping
- strand distribution

It has supported variant QC across thousands of samples in the Pan Prostate Cancer Group ([panprostate.org](https://panprostate.org)) and in published studies of breast tumour mutational signatures. Since it makes no assumptions about the origin of the variants or reads, it is also generalizable beyond human and cancer genomics to other organisms and variant-calling contexts.

## Why a Rust rewrite?

Varlap's original Python implementation queries reads independently for each variant position. This performs well for sparsely distributed somatic variants but scales poorly when variants are closely spaced or when scanning contiguous genomic windows, since overlapping reads are repeatedly refetched from the BAM file.

To address this bottleneck, varlap has been reimplemented in Rust ([jchchiu/rust-varlap](https://github.com/jchchiu/rust-varlap)), restructuring the core algorithm to iterate over reads rather than variant positions, so that each overlapping read is fetched at most once per genomic bin. This sets the foundation for substantially higher performance on dense or contiguous variant sets, while preserving varlap's original output format and metrics.

## Features

- **Two modes**: `loci` (per-variant metrics for a variants file) and `region` (metrics over regions). See [CLI](../2_usage/cli.md).
- **Flexible inputs**: variants as `.vcf`, `.csv` or `.tsv` (plain or gzip-compressed); reads as `.bam` or `.cram` (CRAM requires a FASTA reference).
- **Multiple reads files**: analyse several BAM/CRAM files (e.g. tumour and normal) in one run, with custom labels and optional merged output.
- **Multithreading**: use `--threads` to process multiple BAM/CRAM files in parallel.
- **Output compatibility**: results are written as CSV, preserving varlap's original output format and metrics.