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
    pub variant_file: PathBuf,

    /// Filepaths of BAM files
    #[arg(short, long)]
    pub bam_file: PathBuf,

    /// Type of variants to consider. Options: snv, indel
    #[arg(long, value_enum)]
    pub varclass: VarClass,

    /// Filepath of where csv output should be stored
    #[arg(short, long)]
    pub output_path: String,

    /// Optional sample identifier
    #[arg(long)]
    pub sample: Option<String>,

    /// Optional label for bam files (if not provided will default to name of bam file)
    #[arg(long)]
    pub label: Option<String>,

    /// Required for CRAM: Filepath of FASTA file associated with CRAM file
    #[arg(short, long)]
    pub fasta_file: Option<PathBuf>,
}

pub fn parse() -> Cli {
    Cli::parse()
}