use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

use crate::variant::VarClass;

#[derive(Debug, Parser)]
#[command(name = "rust-varlap")]
#[command(version = "0.1.0-alpha.5")]
#[command(about = "Quality control tool for genetic variants")]
pub struct Cli {
    #[command(subcommand)]
    pub mode: Mode,
}

#[derive(Debug, Subcommand)]
pub enum Mode {
    /// Run the application in Loci mode
    Loci {
        #[command(flatten)]
        common: CommonArgs,

        #[command(flatten)]
        args: LociArgs,
    },

    /// Run the application in Region mode
    Region {
        #[command(flatten)]
        common: CommonArgs,

        #[command(flatten)]
        args: RegionArgs,
    },
}

#[derive(Debug, Args)]
pub struct CommonArgs {
    /// Path to reads file [supported: bam, cram]
    #[arg(short, long, num_args = 1..)]
    pub reads: Vec<PathBuf>,

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

    /// Merge output csvs if multiple BAMs are used
    #[arg(long)]
    pub merge: bool,

    /// Number of threads for multithreading reads file
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(Debug, Args)]
pub struct LociArgs {
    /// Path to variants file [supported: vcf, csv, tsv; optionally gzipped (.gz)]
    #[arg(short, long)]
    pub variants: PathBuf,

    /// Variant class to analyze
    #[arg(short = 'c', long, value_enum)]
    pub varclass: VarClass,
}

#[derive(Debug, Args)]
pub struct RegionArgs {
    /// Path to bed regions file
    #[arg(short = 'R', long)]
    pub regions: PathBuf,

}

pub fn parse() -> Cli {
    Cli::parse()
}
