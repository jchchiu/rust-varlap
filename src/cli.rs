use std::path::PathBuf;

use clap::Parser;

use crate::variant::VarClass;

#[derive(Debug, Parser)]
#[command(name = "rust-varlap")]
#[command(version = "0.0.1")]
#[command(about = "Quality control tool for genetic variants")]
pub struct Cli {
    /// Filepath of variants (Supported: vcf, vcf.gz, csv, tsv)
    #[arg(short, long)]
    pub variants: PathBuf,

    /// Filepath of reads (Supported: bam, cram)
    #[arg(short, long)]
    pub reads: PathBuf,

    /// Class of variants to consider for analysis (Options: snv, indel)
    #[arg(long, value_enum)]
    pub varclass: VarClass,

    /// Filepath of where csv output should be stored and filename of output
    #[arg(short, long)]
    pub output: PathBuf,

    /// Optional sample identifier
    #[arg(long)]
    pub sample: Option<String>,

    /// Optional label for reads file (if not provided will default to name of reads file)
    #[arg(long)]
    pub label: Option<String>,

    /// Filepath of FASTA file associated with CRAM file (Required if reads is CRAM)
    #[arg(short, long)]
    pub fasta: Option<PathBuf>,

    /// Optional bin size gap (if not provided will default to 100kb (100,000bp))
    #[arg(long)]
    pub gap: Option<u64>,
}

pub fn parse() -> Cli {
    Cli::parse()
}