use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("unsupported variants file format")]
    UnsupportedVariantsFormat  {
        filename: PathBuf,
        extension: String,
    },

    #[error("unsupported reads file format")]
    UnsupportedReadsFormat  {
        filename: PathBuf,
        extension: String,
    },

    #[error("missing file extension")]
    MissingExtension {
        filename: PathBuf,
    },

    #[error("invalid gzipped filename")]
    InvalidGzipName {
        filename: PathBuf,
    },

    #[error("CRAM file requires a reference FASTA")]
    MissingCramReference,

    #[error("I/O error")]
    Io(#[from] std::io::Error),
}

impl AppError {
    // Exit codes:
    //  1: File I/O Error
    //  2: Command Line Error
    //  3: File format error
    pub fn exit_code(&self) -> i32 {
        match self {
            AppError::Io(_) => 1,
            AppError::UnsupportedVariantsFormat  { .. } => 3,
            AppError::UnsupportedReadsFormat  { .. } => 3,
            AppError::MissingExtension { .. } => 3,
            AppError::InvalidGzipName { .. } => 3,
            AppError::MissingCramReference => 3,
        }
    }
}

pub fn print_error(program: &str, err: &AppError) {
    match err {
        AppError::UnsupportedVariantsFormat  {
            filename,
            extension,
        } => {
            eprintln!("{program} ERROR: unsupported file format");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Detected extension: {extension}");
            eprintln!("Supported formats are:");
            eprintln!("  .vcf");
            eprintln!("  .vcf.gz");
            eprintln!("  .csv");
            eprintln!("  .tsv");
        }

        AppError::UnsupportedReadsFormat  {
            filename,
            extension,
        } => {
            eprintln!("{program} ERROR: unsupported file format");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Detected extension: {extension}");
            eprintln!("Supported formats are:");
            eprintln!("  .bam");
            eprintln!("  .cram");
        }

        AppError::MissingExtension { filename } => {
            eprintln!("{program} ERROR: missing file extension");
            eprintln!("Input file: {}", filename.display());
            // eprintln!("Expected one of:");
            // eprintln!("  .vcf");
            // eprintln!("  .vcf.gz");
            // eprintln!("  .csv");
            // eprintln!("  .tsv");
        }

        AppError::InvalidGzipName { filename } => {
            eprintln!("{program} ERROR: invalid gzipped filename");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Expected something like:");
            eprintln!("  variants.vcf.gz");
        }

        AppError::MissingCramReference => {
            eprintln!("{program} ERROR: missing reference FASTA");
            eprintln!("CRAM files require a reference FASTA.");
            eprintln!("Specify one with --fasta-file <FASTA>.");
        }

        AppError::Io(err) => {
            eprintln!("{program} ERROR: file I/O error");
            eprintln!("{err}");
        }
    }
}