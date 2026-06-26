use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use tracing::{debug, info, warn};

use crate::errors::AppError;
use crate::features::{LocusFeatures, LocusFeaturesIndel, LocusFeaturesSnv};
use crate::variant::{Variant, VarClass, VarType, ChromBucket, ParsedVariants};

pub fn parse(variants_path: &Path, varclass: &VarClass) -> Result<ParsedVariants> {
    let file_type = detect_file_type(variants_path)
        .with_context(||
            format!(
                "Failed to detect variants file type for '{}'",
                variants_path.display()
            )
        )?;

    let reader = check_valid_gzip(variants_path)
        .with_context(||
            format!("Failed to open reader for {}", variants_path.display()))?;

    info!(
        "Parsing variants from {} as {:?}",
        variants_path.display(),
        file_type
    );

    // let mut variants = VecDeque::new();

    let mut buckets: Vec<ChromBucket> = Vec::new();

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.with_context(|| {
            format!("Failed reading line {} from {}", line_no + 1, variants_path.display())
        })?;

        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = match file_type {
            FileType::Vcf | FileType::Tsv => line.split_whitespace().collect(),
            FileType::Csv => line.split(',').collect(),
        };

        if fields.len() < 5 {
            warn!("Skipping ipnut row at line {}: {}", line_no + 1, line);
            continue;
        }

        let chrom = fields[0].to_string();

        let pos = match fields[1].parse::<u64>() {
            Ok(p) => p,
            Err(err) => {
                warn!(
                    "Invalid POS at line {}: {} ({err})",
                    line_no + 1,
                    line
                );
                continue;
            }
        };

        let refr = fields[3].to_string();

        for alt in fields[4].split(',') {
            let vartype = get_var_type(&refr, alt);

            if is_acceptable_variant(varclass, &vartype, &refr, alt) {
                let variant = match varclass {
                    VarClass::Snv => Variant {
                        chrom: chrom.clone(),
                        pos,
                        refr: refr.clone(),
                        alt: alt.to_string(),
                        vartype,
                        features: LocusFeatures::Snv(LocusFeaturesSnv::default()),
                    },
                    VarClass::Indel => Variant {
                        chrom: chrom.clone(),
                        pos,
                        refr: refr.clone(),
                        alt: alt.to_string(),
                        vartype,
                        features: LocusFeatures::Indel(LocusFeaturesIndel::default()),
                    },
                };
                // variants.push_back(variant);

                // Get unique chromosomes and their counts for variants addded to queue
                // NOTE: Variants file MUST be in sorted ascending order
                // We do not need to get the index as we are popping the queue when iterating over variants
                if let Some(last) = buckets.last_mut() {
                    if last.chrom == chrom {
                        last.variants.push_back(variant);
                        continue;
                    }
                }

                buckets.push(ChromBucket {
                    chrom: chrom.clone(),
                    variants: VecDeque::from([variant]),
                });
            } else {
                debug!(
                    "Skipped invalid variant at line {}: chrom={} pos={} ref={} alt={}",
                    line_no + 1,
                    chrom,
                    pos,
                    refr,
                    alt
                );
            }
        }
    }

    // info!("Parsed {} variants from {}", variants.len(), variants_path.display());
    Ok(ParsedVariants { chroms: buckets })
}

#[derive(Debug, Clone, Copy)]
enum FileType {
    Vcf,
    Csv,
    Tsv,
}

fn detect_file_type(path: &Path) -> Result<FileType, AppError> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .ok_or_else(|| AppError::MissingExtension {
            filename: path.to_path_buf(),
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
        "vcf" => Ok(FileType::Vcf),
        "csv" => Ok(FileType::Csv),
        "tsv" => Ok(FileType::Tsv),
        _ => Err(AppError::UnsupportedVariantsFormat {
            filename: path.to_path_buf(),
            extension: actual_ext.to_string(),
        }),
    }
}

fn get_var_type(refr: &str, alt: &str) -> VarType {
    if refr.len() == 1 && alt.len() == 1 {
        VarType::Snv
    } else if refr.len() > alt.len() {
        VarType::Del
    } else if refr.len() < alt.len() {
        VarType::Ins
    } else {
        warn!("Cannot determine the type of variant with ref:{} and alt:{}", refr, alt);
        VarType::Unknown
    }
}

// USING GZIP HANDLER AS PER RUSTQC
// https://github.com/seqeralabs/RustQC/blob/main/src/io.rs
// REWRITE MYSELF LATER

// https://cseweb.ucsd.edu/classes/sp22/cse223B-a/tribbler/flate2/read/struct.MultiGzDecoder.html
// Need to use multigzdecoder

/// Gzip magic bytes: the first two bytes of any gzip-compressed file.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

fn check_valid_gzip(file_path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(file_path)
        .with_context(|| format!("Failed to open variants file '{}'", file_path.display()))?;

    // Read the first two bytes to check for gzip magic number
    let mut magic = [0u8; 2];
    let bytes_read = file
        .read(&mut magic)
        .with_context(|| format!("Failed to read header from '{}'", file_path.display()))?;

    // Seek back to the beginning so the reader starts from byte 0
    file.seek(SeekFrom::Start(0))
        .with_context(|| format!("Failed to reseek start of file '{}'", file_path.display()))?;

    if bytes_read >= 2 && magic == GZIP_MAGIC {
        debug!("Detected gzip-compressed input: {}", file_path.display());
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        debug!("Detected plain-text input: {}", file_path.display());
        Ok(Box::new(BufReader::new(file)))
    }
}

// fn is_valid_vcf_header_line(line: &str) -> bool {
//     let expected = ["#CHROM", "POS", "ID", "REF", "ALT"];
//     line.split_whitespace().take(5).eq(expected)
// }

fn is_acceptable_variant(
    varclass: &VarClass,
    vartype: &VarType,
    refr: &str,
    alt: &str,
    // max_indel_size: u32,
) -> bool {
    if !is_only_dna_bases(refr) || !is_only_dna_bases(alt) {
        false
    } else if !is_desired_type(varclass, vartype) {
        false
    // } else if !is_within_max_size(varclass, max_indel_size, refr, alt) {
    //     false
    } else if matches!(varclass, VarClass::Indel) && !is_valid_indel(refr, alt) {
        false
    } else {
        true
    }
}

fn is_only_dna_bases(sequence: &str) -> bool {
    sequence
        .chars()
        .all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'T' | 'G' | 'C'))
}

fn is_desired_type(varclass: &VarClass, vartype: &VarType) -> bool {
    match varclass {
        VarClass::Snv => matches!(vartype, VarType::Snv),
        VarClass::Indel => matches!(vartype, VarType::Ins | VarType::Del),
    }
}

// fn is_within_max_size(
//     varclass,
//     max_indel_size,
//     refr,
//     alt
// ) -> bool {

// }

fn is_valid_indel(refr: &str, alt: &str) -> bool {
    match refr.len().cmp(&alt.len()) {
        std::cmp::Ordering::Equal => false,
        std::cmp::Ordering::Less => !refr.is_empty() && alt.starts_with(refr),
        std::cmp::Ordering::Greater => !alt.is_empty() && refr.starts_with(alt),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_dna_bases() {
        // Valid DNA base lowercase
        assert!(is_only_dna_bases("a"));
        // Valid DNA base uppercase
        assert!(is_only_dna_bases("A"));
        // Invalid DNA base
        assert!(!is_only_dna_bases("K"));
        // Invalid DNA string
        assert!(!is_only_dna_bases("string"));
    }
}