use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use tracing::{debug, info, warn};

use crate::errors::AppError;
use crate::reader;
use crate::variant::{ParsedVariants, VarType, VariantInfo};

pub fn parse(regions_path: &Path) -> Result<ParsedVariants> {
    let file_type = detect_file_type(regions_path)?;

    info!(
        "Parsing regions from {:?} as: {:?}",
        regions_path.display(),
        file_type
    );

    let mut positions: Vec<VariantInfo> = Vec::new();

    let mut skipped = SkippedVariants::default();

    info!("Parsed regions successfully");
    info!(
        "Total number of positions in input: {}",
        skipped.total_variants(positions.len())
    );
    info!("Number of positions parsed: {}", positions.len());
    info!("Number of positions skipped: {}", skipped.skipped);

    Ok(ParsedVariants { positions })
}

#[derive(Debug, Clone, Copy)]
enum FileType {
    Bed,
    Gatk,
    Picard,
}

fn detect_file_type(path: &Path) -> Result<FileType, AppError> {
    let ext = path.extension().and_then(|e| e.to_str()).ok_or_else(|| {
        AppError::MissingRegionsExtension {
            filename: path.to_path_buf(),
        }
    })?;

    let actual_ext = if ext == "gz" {
        path.file_stem()
            .and_then(|s| Path::new(s).extension())
            .and_then(|e| e.to_str())
            .ok_or_else(|| AppError::InvalidGzipName {
                filename: path.to_path_buf(),
            })?
    } else {
        ext
    };

    match actual_ext {
        "bed" => Ok(FileType::Bed),
        "list" | "intervals" => Ok(FileType::Gatk),
        "interval_list" => Ok(FileType::Picard),
        _ => Err(AppError::UnsupportedRegionsFormat {
            filename: path.to_path_buf(),
            extension: actual_ext.to_string(),
        }),
    }
}

positions.push(VariantInfo {
                        chrom: row.chrom.clone(),
                        pos: norm_pos,
                        refr: Some(norm_refr),
                        alt: Some(norm_alt),
                        vartype,
                    }),