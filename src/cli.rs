use clap::Parser;
use std::path::PathBuf;

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
    #[arg(short = 'o', long)]
    pub csv_path: String,

    /// Optional sample identifier
    #[arg(long)]
    pub sample: Option<String>,
}

pub fn parse() -> Cli {
    Cli::parse()
}