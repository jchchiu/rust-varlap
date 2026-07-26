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

    #[error("variants file is not sorted")]
    UnsortedVariants  {
        chromosome: String,
        error_pos: u64,
        previous_pos: u64,
    },

    #[error("missing variants file extension")]
    MissingVariantsExtension {
        filename: PathBuf,
    },

    #[error("missing reads file extension")]
    MissingReadsExtension {
        filename: PathBuf,
    },

    #[error("invalid gzipped filename")]
    InvalidGzipName {
        filename: PathBuf,
    },

    #[error("missing required CSV/TSV Header")]
    MissingDelimitedHeader {
        fields: String,
        headers: String,
    },

    // #[error("invalid VCF Header")]
    // InvalidVcfHeader {
    //     header: String,
    // },

    #[error("missing reference FASTA required for CRAM file")]
    MissingCramReference,

    #[error("reference sequence not found")]
    MissingReferenceSequence {
        chromosome: String,
    },
}

impl AppError {
    // Exit codes:
    //  1: File I/O Error
    //  2: Command Line Error
    //  3: File format error
    pub fn exit_code(&self) -> i32 {
        match self {
            AppError::UnsupportedVariantsFormat  { .. } => 3,
            AppError::UnsupportedReadsFormat  { .. } => 3,
            AppError::UnsortedVariants { ..} => 3,
            AppError::MissingVariantsExtension { .. } => 3,
            AppError::MissingReadsExtension { .. } => 3,
            AppError::InvalidGzipName { .. } => 3,
            AppError::MissingDelimitedHeader { .. } => 3,
            AppError::MissingCramReference => 3,
            AppError::MissingReferenceSequence { .. } => 3,
        }
    }
}

pub fn print_error(program: &str, err: &AppError) {
    match err {
        AppError::UnsupportedVariantsFormat  {
            filename,
            extension,
        } => {
            eprintln!("{program} ERROR: unsupported variants file format");
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
            eprintln!("{program} ERROR: unsupported reads file format");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Detected extension: {extension}");
            eprintln!("Supported formats are:");
            eprintln!("  .bam");
            eprintln!("  .cram");
        }

        AppError::UnsortedVariants {
            chromosome,
            error_pos,
            previous_pos
        } => {
            eprintln!("{program} ERROR: variants are not sorted");
            eprintln!("At chromosome: {}", chromosome);
            eprintln!("Variant at position {} comes before {}", error_pos, previous_pos);
            eprintln!("Please sort positions in ascending order before running again");
        }

        AppError::MissingVariantsExtension { filename } => {
            eprintln!("{program} ERROR: variants missing file extension");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Expected one of:");
            eprintln!("  .vcf");
            eprintln!("  .vcf.gz");
            eprintln!("  .csv");
            eprintln!("  .tsv");
        }

        AppError::MissingReadsExtension { filename } => {
            eprintln!("{program} ERROR: reads missing file extension");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Expected one of:");
            eprintln!("  .bam");
            eprintln!("  .cram");
        }

        AppError::InvalidGzipName { filename } => {
            eprintln!("{program} ERROR: invalid gzipped filename");
            eprintln!("Input file: {}", filename.display());
            eprintln!("Expected something like:");
            eprintln!("  variants.vcf.gz");
        }

        AppError::MissingDelimitedHeader { fields, headers } => {
            eprintln!("{program} ERROR: invalid csv/tsv header");
            eprintln!("Missing required header");
            eprintln!("Expected one of: {}", fields);
            eprintln!("Found headers: {}", headers);
        }

        AppError::MissingCramReference => {
            eprintln!("{program} ERROR: missing reference FASTA");
            eprintln!("CRAM files require the exact reference FASTA that was used for it's creation.");
            eprintln!("Specify one with --fasta <FASTA>.");
        }

        AppError::MissingReferenceSequence { chromosome } => {
            eprintln!("{program} ERROR: reference sequence not found");
            eprintln!("Chromosome: {}", chromosome);
            eprintln!("The alignment file does not contain this reference sequence.");
            eprintln!("Check that the BAM/CRAM and variant file use the same reference genome.");
        }
    }
}