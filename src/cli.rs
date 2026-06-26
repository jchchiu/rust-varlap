use std::path::PathBuf;

use clap::Parser;

use crate::variant::VarClass;

#[derive(Debug, Parser)]
#[command(name = "rust-varlap")]
#[command(version = "0.0.1")]
#[command(about = "Quality control tool for genetic variants")]
pub struct Cli {
    /// Input variant file path (can be vcf, csv, tsv and be gzipped (.gz))
    #[arg(short, long)]
    pub variants: PathBuf,

    /// Filepaths of read files (can be bam or cram)
    #[arg(short, long)]
    pub reads: PathBuf,

    /// Type of variants to consider. Options: snv, indel
    #[arg(long, value_enum)]
    pub varclass: VarClass,

    /// Filepath of where csv output should be stored and filename of output
    #[arg(short, long)]
    pub output: PathBuf,

    /// Optional sample identifier
    #[arg(long)]
    pub sample: Option<String>,

    /// Optional label for bam files (if not provided will default to name of bam file)
    #[arg(long)]
    pub label: Option<String>,

    /// Required for CRAM: Filepath of FASTA file associated with CRAM file
    #[arg(short, long)]
    pub fasta: Option<PathBuf>,
}

pub fn parse() -> Cli {
    Cli::parse()
}