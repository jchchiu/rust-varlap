use std::path::PathBuf;

use clap::Parser;

use crate::variant::VarClass;

#[derive(Debug, Parser)]
#[command(name = "rust-varlap")]
#[command(version = "0.1.0-alpha.1")]
#[command(about = "Quality control tool for genetic variants")]
pub struct Cli {
    /// Path to variants file [supported: vcf, csv, tsv; optionally gzipped (.gz)]
    #[arg(short, long)]
    pub variants: PathBuf,

    /// Path to reads file [supported: bam, cram]
    #[arg(short, long, num_args = 1..)]
    pub reads: Vec<PathBuf>,

    /// Variant class to analyze
    #[arg(short = 'c', long, value_enum)]
    pub varclass: VarClass,

    /// Path to csv output directory and filename of output
    #[arg(short, long)]
    pub output: PathBuf,

    /// Path to FASTA reference for CRAM input (required if reads is CRAM)
    #[arg(short, long)]
    pub fasta: Option<PathBuf>,

    /// Sample identifier
    #[arg(long)]
    pub sample: Option<String>,

    /// Label for reads file (defaults to reads filename)
    #[arg(long, num_args = 1..)]
    pub label: Vec<Option<String>>,

    /// Bin size gap in base pairs (defaults to 100,000)
    #[arg(long)]
    pub gap: Option<u64>,
}

pub fn parse() -> Cli {
    Cli::parse()
}
