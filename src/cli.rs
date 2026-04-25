use clap::Parser;

use crate::variant::VarClass;

#[derive(Debug, Parser)]
#[command(name = "rust-varlap")]
#[command(version = "0.0.1")]
#[command(about = "Quality control tool for genetic variants")]
pub struct Cli {
    /// Input VCF file path
    #[arg(short, long)]
    pub vcf: String,

    /// Filepaths of BAM files
    #[arg(short, long)]
    pub bams: String,

    /// Type of variants to consider. Options: snv, indel
    #[arg(long, value_enum)]
    pub varclass: VarClass,

    /// Filepath of where csv output should be stored
    #[arg(short, long)]
    pub csv_path: String,

    /// Optional sample identifier
    #[arg(long)]
    pub sample: Option<String>,
}

pub fn parse() -> Cli {
    Cli::parse()
}