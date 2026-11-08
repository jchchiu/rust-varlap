use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use tracing::{info, warn};

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

    let mut parsed = RegionsParsed::default();

    match file_type {
        FileType::Bed => parse_bed(regions_path, &mut positions, &mut parsed)?,
        // To be implemented
        FileType::Gatk => unreachable!(), //parse_gatk(regions_path, &mut positions, &mut skipped)?,
        FileType::Picard => unreachable!(), //parse_picard(regions_path, &mut positions, &mut skipped)?,
    };

    info!("Parsed regions successfully");
    info!(
        "Total number of regions in input: {}",
        parsed.total_regions()
    );
    info!("Number of regions processed: {}", parsed.processed);
    info!("Number of regions skipped: {}", parsed.skipped);
    info!("Number of positions parsed: {}", positions.len());
    
    Ok(ParsedVariants { variants: positions })
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

#[derive(Debug, Clone)]
struct RegionRow {
    chrom: String,
    start: u64,
    end: u64,
}

#[derive(Default, Debug, Clone)]
struct RegionsParsed {
    processed: usize,
    skipped: usize,
}

impl RegionsParsed {
    fn total_regions(&self) -> usize {
        self.skipped + self.processed
    }
}

// ===================================================================
// Parsing BEDs
// ===================================================================

fn parse_bed(
    file_path: &Path,
    positions: &mut Vec<VariantInfo>,
    parsed: &mut RegionsParsed,
) -> Result<()> {
    // Fix error statement for open input later
    let input = reader::open_variant_input(file_path)?;
    let reader = BufReader::new(input);

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.with_context(|| format!("Failed reading line {}", line_no + 1))?;

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 3 {
            warn!("Skipping input row at line {}: {}", line_no + 1, line);
            parsed.skipped += 1;
            continue;
        }

        // BED is 0-based, start inclusive, end exclusive
        let start = match fields[1].parse::<u64>() {
                Ok(p) => p + 1,
                Err(err) => {
                    warn!("Invalid start POS at line {}: {} ({err})", line_no + 1, line);
                    parsed.skipped += 1;
                    continue;
                }
            };

        let end = match fields[2].parse::<u64>() {
                Ok(p) => p + 1,
                Err(err) => {
                    warn!("Invalid end POS at line {}: {} ({err})", line_no + 1, line);
                    parsed.skipped += 1;
                    continue;
                }
            };

        // Can change this later to is_acceptable_region function check
        let row = if (start > 0) && (end > 0) && (end > start) {
            RegionRow {
                chrom: fields[0].to_string(),
                start,
                end,
            }
        } else {
            warn!("Start/end pos is > 0 or end pos > start pos at line {}: {}", line_no + 1, line);
            parsed.skipped += 1;
            continue;            
        };

        process_region_row(&row, positions)?;
        parsed.processed += 1;
    }

    Ok(())
}

// ===================================================================
// Processing Region Rows
// ===================================================================

fn process_region_row(
    row: &RegionRow,
    positions: &mut Vec<VariantInfo>,
) -> Result<()> {
    for i in row.start..row.end {
        positions.push(VariantInfo {
            // NOTE: waste here to clone, can probably map chromosomes to ints
            //       or as a shared reference
            chrom: row.chrom.clone(),
            pos: i,
            refr: None,
            alt: None,
            vartype: VarType::Unknown,
        })
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn parse_bed_line(line: &str) -> (Vec<u64>, RegionsParsed) {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "{line}").unwrap();

        let mut positions = Vec::new();
        let mut parsed = RegionsParsed::default();
        parse_bed(file.path(), &mut positions, &mut parsed).unwrap();

        (positions.iter().map(|p| p.pos).collect(), parsed)
    }

    #[test]
    fn parse_bed_coordinate_cases() {
        let cases = [
            // BED is 0-based half-open; positions are 1-based inclusive
            ("chr1\t0\t10", (1..=10).collect::<Vec<u64>>()),
            // Single base region
            ("chr1\t4\t5", vec![5]),
            // Two base region
            ("chr1\t4\t6", vec![5, 6]),
            // Extra BED columns are ignored
            ("chr1\t4\t6\tname\t0\t+", vec![5, 6]),
            // Empty region (end == start) is skipped
            ("chr1\t5\t5", vec![]),
            // End before start is skipped
            ("chr1\t6\t5", vec![]),
            // End of 0 is skipped (no underflow)
            ("chr1\t0\t0", vec![]),
        ];

        for (line, expected) in cases {
            let (positions, parsed) = parse_bed_line(line);

            assert_eq!(positions, expected, "line={line:?}");
            assert_eq!(
                parsed.processed + parsed.skipped,
                1,
                "line={line:?}"
            );
            assert_eq!(
                parsed.skipped,
                usize::from(expected.is_empty()),
                "line={line:?}"
            );
        }
    }

    #[test]
    fn parse_bed_invalid_row_cases() {
        let cases = [
            // Too few columns
            "chr1\t0",
            // Non-numeric start
            "chr1\tx\t10",
            // Non-numeric end
            "chr1\t0\tx",
            // Negative coordinates
            "chr1\t-1\t10",
        ];

        for line in cases {
            let (positions, parsed) = parse_bed_line(line);

            assert!(positions.is_empty(), "line={line:?}");
            assert_eq!(parsed.skipped, 1, "line={line:?}");
            assert_eq!(parsed.processed, 0, "line={line:?}");
        }
    }

    #[test]
    fn parse_bed_region_positions_have_no_alleles() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "chrM\t99\t100").unwrap();

        let mut positions = Vec::new();
        let mut parsed = RegionsParsed::default();
        parse_bed(file.path(), &mut positions, &mut parsed).unwrap();

        assert_eq!(positions.len(), 1);
        assert_eq!(positions[0].chrom, "chrM");
        assert_eq!(positions[0].pos, 100);
        assert!(positions[0].refr.is_none());
        assert!(positions[0].alt.is_none());
        assert_eq!(positions[0].vartype, VarType::Unknown);
    }
}
