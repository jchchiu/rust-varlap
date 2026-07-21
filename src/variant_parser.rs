use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};
use csv::{ReaderBuilder, StringRecord};
use flate2::read::MultiGzDecoder;
use tracing::{debug, info, warn};

use crate::errors::AppError;
use crate::variant::{VariantInfo, VarClass, VarType, ParsedVariants};

pub fn parse(variants_path: &Path, varclass: &VarClass) -> Result<ParsedVariants> {
    let file_type = detect_file_type(variants_path)?;

    info!(
        "Parsing variants from {:?} as: {:?}",
        variants_path.display(),
        file_type
    );

    let mut variants: Vec<VariantInfo> = Vec::new();

    let mut skipped = SkippedVariants::default();

    match file_type {
        FileType::Vcf => parse_vcf(variants_path, varclass, &mut variants, &mut skipped)?,
        FileType::Csv | FileType::Tsv => parse_delimited(variants_path, file_type, varclass, &mut variants, &mut skipped)?,
    };

    info!("Parsed variants successfully");
    info!("Total number of variants in input: {}", skipped.total_variants(variants.len()));
    info!("Number of variants parsed: {}", variants.len());
    info!("Number of variants skipped: {}", skipped.skipped);
    
    Ok(ParsedVariants { variants })
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
        .ok_or_else(|| AppError::MissingVariantsExtension {
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

#[derive(Default, Debug, Clone)]
struct SkippedVariants {
    skipped: usize,
}

impl SkippedVariants {
    fn total_variants(&self, parsed_variants: usize) -> usize{
        self.skipped + parsed_variants
    }
}

#[derive(Debug, Clone)]
struct VariantRow {
    chrom: String,
    pos: u64,
    refr: String,
    alts: Vec<String>,
    line_no: usize,
}

enum MaybeGzipReader {
    Plain(File),
    Gzip(MultiGzDecoder<File>),
}

impl Read for MaybeGzipReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            MaybeGzipReader::Plain(file) => file.read(buf),
            MaybeGzipReader::Gzip(decoder) => decoder.read(buf),
        }
    }
}

fn open_variant_input(path: &Path) -> Result<MaybeGzipReader> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open variants file '{}'", path.display()))?;

    if is_gzip(&mut file)? {
        Ok(MaybeGzipReader::Gzip(MultiGzDecoder::new(file)))
    } else {
        Ok(MaybeGzipReader::Plain(file))
    }
}

fn parse_vcf(
    file_path: &Path,
    varclass: &VarClass,
    variants: &mut Vec<VariantInfo>,
    skipped: &mut SkippedVariants,
) -> Result<()> {
    let input = open_variant_input(file_path)?;
    let reader = BufReader::new(input);

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.with_context(|| {
            format!("Failed reading line {}", line_no + 1)
        })?;

        // if line.starts_with("##") {
        //     is_valid_vcf_header_line(&line)?;
        // }

        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 5 {
            warn!("Skipping input row at line {}: {}", line_no + 1, line);
            skipped.skipped += 1;
            continue;
        }

        let row = VariantRow {
            chrom: fields[0].to_string(),
            pos: match fields[1].parse() {
                Ok(p) => p,
                Err(err) => {
                    warn!("Invalid POS at line {}: {} ({err})", line_no + 1, line);
                    skipped.skipped += 1;
                    continue;
                }
            },
            refr: fields[3].to_string(),
            alts: fields[4].split(',').map(|s| s.to_string()).collect(),
            line_no: line_no + 1,
        };

        process_variant_row(&row, varclass, variants, skipped)?;
    }

    Ok(())
}

fn delimiter_for(file_type: FileType) -> u8 {
    match file_type {
        FileType::Csv => b',',
        FileType::Tsv => b'\t',
        FileType::Vcf => unreachable!(),
    }
}

fn get_header_index(headers: &StringRecord, fields: &[&str]) -> Result<usize, AppError> {
    for field in fields {
        if let Some(idx) = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case(field)) {
            return Ok(idx);
        }
    }

    Err(AppError::MissingDelimitedHeader {
        fields: fields.join(", "),
        headers: headers.iter().collect::<Vec<&str>>().join(", "),
    })
}

fn parse_delimited(
    file_path: &Path,
    file_type: FileType,
    varclass: &VarClass,
    variants: &mut Vec<VariantInfo>,
    skipped: &mut SkippedVariants,
) -> Result<()> {
    let input = open_variant_input(file_path)?;

    let mut csv_reader = ReaderBuilder::new() 
        .delimiter(delimiter_for(file_type))
        .has_headers(true)
        .flexible(true)
        .from_reader(input);

    let headers = csv_reader
        .headers()
        .with_context(|| format!("Failed to read headers from '{}'", file_path.display()))?
        .clone();

    let chrom_idx = get_header_index(&headers, &["chrom", "chr", "#chrom"])?;
    let pos_idx = get_header_index(&headers, &["pos", "position"])?;
    let ref_idx = get_header_index(&headers, &["ref", "refr", "reference"])?;
    let alt_idx = get_header_index(&headers, &["alt", "alts", "alternate"])?;

    for (record_no, record_result) in csv_reader.records().enumerate() {
        let record = record_result.with_context(|| {
            format!("Failed reading record {}", record_no + 2)
        })?;

        let line_no = record_no + 2;

        let row = VariantRow {
            chrom: match record.get(chrom_idx) {
                Some(v) if !v.trim().is_empty() => v.trim().to_string(),
                _ => {
                    warn!("Missing CHROM at line {}", line_no);
                    skipped.skipped += 1;
                    continue;
                }
            },
            pos: match record.get(pos_idx).map(str::trim).unwrap_or("").parse() {
                Ok(p) => p,
                Err(err) => {
                    warn!("Invalid POS at line {}: {:?} ({err})", line_no, record);
                    skipped.skipped += 1;
                    continue;
                }
            },
            refr: match record.get(ref_idx) {
                Some(v) if !v.trim().is_empty() => v.trim().to_string(),
                _ => {
                    warn!("Missing REF at line {}", line_no);
                    skipped.skipped += 1;
                    continue;
                }
            },
            alts: match record.get(alt_idx) {
                Some(v) if !v.trim().is_empty() => {
                    v.split(',').map(|s| s.trim().to_string()).collect()
                }
                _ => {
                    warn!("Missing ALT at line {}", line_no);
                    skipped.skipped += 1;
                    continue;
                }
            },
            line_no,
        };

        process_variant_row(&row, varclass, variants, skipped)?;
    }

    Ok(())
}

fn process_variant_row(
    row: &VariantRow, 
    varclass: &VarClass, 
    variants: &mut Vec<VariantInfo>,
    skipped: &mut SkippedVariants,
)-> Result<()> {
    for alt in &row.alts {
        let vartype = get_var_type(&row.refr, alt);

        if is_acceptable_variant(varclass, &vartype, row, alt) {
            let variant = match varclass {
                VarClass::Snv => VariantInfo {
                    chrom: row.chrom.clone(),
                    pos: row.pos,
                    refr: row.refr.clone(),
                    alt: alt.to_string(),
                    vartype,
                },
                VarClass::Indel => VariantInfo {
                    chrom: row.chrom.clone(),
                    pos: row.pos,
                    refr: row.refr.clone(),
                    alt: alt.to_string(),
                    vartype,
                },
            };
            
            variants.push(variant);
        } else {
            skipped.skipped += 1;
        }
    }

    Ok(())
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

fn is_gzip(file: &mut File) -> Result<bool> {
    let mut magic = [0; 2];
    let bytes = file
        .read(&mut magic)
        .with_context(|| format!("Failed to read first two bytes from '{:?}'", file))?;

    // Go back to the start of the file
    file.seek(SeekFrom::Start(0))
        .with_context(|| format!("Failed to reseek start of file '{:?}'", file))?;

    Ok(bytes >= 2 && magic == GZIP_MAGIC)
}

// fn is_valid_vcf_header_line(line: &str) -> Result<(), AppError> {
//     let expected = ["#CHROM", "POS", "ID", "REF", "ALT"];
//     if line.split('\t').take(5).eq(expected){
//         return Err(AppError::InvalidVcfHeader {
//             header: line.split('\t').collect(),
//         })
//     } else {
//         Ok(())
//     }
// }

fn is_acceptable_variant(
    varclass: &VarClass,
    vartype: &VarType,
    row: &VariantRow,
    alt: &str,
    // max_indel_size: u32,
) -> bool {
    if !is_only_dna_bases(&row.refr) || !is_only_dna_bases(alt) {
        debug!(
            "Skipped invalid variant at line {}: chrom={} pos={} \n 
            ref={} or alt={} contains non DNA bases (a, c, t, g)",
            row.line_no + 1,
            row.chrom,
            row.pos,
            row.refr,
            alt,
        );
        
        false
    } else if !is_desired_type(varclass, vartype) {
        debug!(
            "Skipped invalid variant at line {}: chrom={} pos={} ref={} alt={} \n 
            varclass={:?} and vartype={:?} do not match",
            row.line_no + 1,
            row.chrom,
            row.pos,
            row.refr,
            alt,
            varclass,
            vartype,
        );

        false
    // } else if !is_within_max_size(varclass, max_indel_size, refr, alt) {
    //     false
    } else if matches!(varclass, VarClass::Indel) && !is_valid_indel(&row.refr, alt) {
        debug!(
            "Skipped invalid variant at line {}: chrom={} pos={} vartype={:?} \n 
            ref={} is not a valid indel",
            row.line_no + 1,
            row.chrom,
            row.pos,
            vartype,
            row.refr,
        );
        
        false
    } else {
        true
    }
}

fn is_only_dna_bases(sequence: &str) -> bool {
    !sequence.is_empty()
        && sequence
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'T' | b'G' | b'C'))
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
    fn is_only_dna_bases_cases() {
        let cases = [
            // Valid uppercase/lowercase
            ("A", true),
            ("a", true),
            ("AaTtGgCc", true),

            // Invalid bases
            ("N", false),
            ("?", false),
            ("A-T", false),
            ("A1T", false),
            ("", false),
        ];

        for (sequence, expected) in cases {
            assert_eq!(
                is_only_dna_bases(sequence),
                expected,
                "sequence={sequence:?}"
            );
        }
    }

    #[test]
    fn is_desired_type_cases() {
        let cases = [
            (VarClass::Snv,   VarType::Snv, true),
            (VarClass::Snv,   VarType::Ins, false),
            (VarClass::Snv,   VarType::Del, false),
            (VarClass::Indel, VarType::Snv, false),
            (VarClass::Indel, VarType::Ins, true),
            (VarClass::Indel, VarType::Del, true),
        ];

        for (varclass, vartype, expected) in cases {
            assert_eq!(
                is_desired_type(&varclass, &vartype),
                expected,
                "varclass={varclass:?}, vartype={vartype:?}"
            );
        }
    }

    #[test]
    fn is_valid_indel_cases() {
        let cases = [
            // Valid insertions
            ("A", "AT", true),
            ("AT", "ATG", true),

            // Valid deletions
            ("AT", "A", true),
            ("ATG", "AT", true),

            // Equal length
            ("A", "T", false),
            ("AT", "GC", false),

            // Invalid insertion
            ("A", "GA", false),
            ("AT", "GAT", false),

            // Invalid deletion
            ("GA", "A", false),
            ("GAT", "AT", false),

            // Empty alleles
            ("", "A", false),
            ("A", "", false),
        ];

        for (refr, alt, expected) in cases {
            assert_eq!(
                is_valid_indel(refr, alt),
                expected,
                "REF={refr:?}, ALT={alt:?}"
            );
        }
    }

    #[test]
    fn is_acceptable_variant_cases() {
        let cases = [
            // Valid SNV
            (VarClass::Snv, VarType::Snv, "A", "T", true),
            // Valid insertion
            (VarClass::Indel, VarType::Ins, "A", "AT", true),
            // Valid deletion
            (VarClass::Indel, VarType::Del, "AT", "A", true),
            // Invalid reference base
            (VarClass::Snv, VarType::Snv, "N", "T", false),
            // Invalid alternate base
            (VarClass::Snv, VarType::Snv, "A", "N", false),
            // Wrong type
            (VarClass::Snv, VarType::Ins, "A", "AT", false),
            // Malformed insertion
            (VarClass::Indel, VarType::Ins, "A", "GA", false),
            // Malformed deletion
            (VarClass::Indel, VarType::Del, "GA", "A", false),
        ];

        for (varclass, vartype, refr, alt, expected) in cases {
            let row = VariantRow {
                chrom: "chr1".to_string(),
                pos: 1,
                refr: refr.to_string(),
                alts: vec![alt.to_string()],
                line_no: 0,
            };

            assert_eq!(
                is_acceptable_variant(&varclass, &vartype, &row, alt),
                expected,
                "varclass={varclass:?}, vartype={vartype:?}, REF={refr:?}, ALT={alt:?}"
            );
        }
    }
}