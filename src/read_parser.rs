use std::collections::VecDeque;
use std::fs::File;
use std::path::Path;
// use std::rc::Rc;

use anyhow::{Context, Result};
use csv::{Writer, WriterBuilder};
use rust_htslib::bam::{IndexedReader, Read, Record};
use tracing::{debug, info};

use crate::errors::AppError;
use crate::output::{write_header, write_variant_row};
use crate::variant::{BinnedVariants, VarClass, Variant, VariantBin};

pub fn parse(
    binned_variants: &mut BinnedVariants,
    reads_path: &Path,
    csv_path: &Path,
    sample: Option<&str>,
    label: Option<&str>,
    varclass: &VarClass,
    fasta_path: Option<&Path>,
    // ADD vector of [chrom/first variant index/length] here
) -> Result<()> {
    let mut reader = open_indexed_reader(reads_path, fasta_path)?;

    let mut csv_writer = WriterBuilder::new()
        .has_headers(false)
        .buffer_capacity(crate::output::WRITE_BUF_SIZE)
        .from_path(csv_path)
        .with_context(|| format!("Failed to create output CSV '{}'", csv_path.display()))?;

    // Write dynamic header based on BAM/CRAM filename
    write_header(&mut csv_writer, reads_path, label, varclass)
        .context("Could not CSV write header")?;

    for bin in &mut binned_variants.bins {
        process_bin(bin, &mut reader, &mut csv_writer, sample)?;

        csv_writer
            .flush()
            .context("Failed to flush output CSV to disk")?;
    }

    info!("Parsed reads successfully");
    info!("CSV output can be found at: {:?}", csv_path.display());

    Ok(())
}

fn open_indexed_reader(reads_path: &Path, fasta_path: Option<&Path>) -> Result<IndexedReader> {
    let file_type = detect_file_type(reads_path)?;

    let mut reader = IndexedReader::from_path(reads_path)
        .with_context(|| format!("Failed to open reads file '{}'", reads_path.display()))?;

    match file_type {
        FileType::Bam => {}
        FileType::Cram => {
            if fasta_path.is_some() {
                let fasta = fasta_path.ok_or(AppError::MissingCramReference)?;
                reader.set_reference(fasta).with_context(|| {
                    format!("Failed to set CRAM reference to '{}'", fasta.display())
                })?;
            }
        }
    }

    info!(
        "Parsing reads from {:?} as: {:?}",
        reads_path.display(),
        file_type
    );

    Ok(reader)
}

fn process_bin(
    bin: &mut VariantBin,
    reader: &mut IndexedReader,
    csv_writer: &mut Writer<File>,
    sample: Option<&str>,
) -> Result<()> {
    let chrom_info = get_chrom_info(&bin.variants);

    reader
        .fetch((&bin.chrom, chrom_info.min_pos - 1, chrom_info.max_pos))
        .with_context(|| {
            format!(
                "Failed to fetch region {}:{}-{}",
                bin.chrom, chrom_info.min_pos, chrom_info.max_pos,
            )
        })?;

    debug!(
        "Processing bin from chromosome: {} spanning {}-{}, number of variants {}",
        bin.chrom,
        chrom_info.min_pos,
        chrom_info.max_pos,
        bin.variants.len(),
    );

    let ref_seq_len = get_ref_len(reader, &bin.chrom)?;

    for read_result in reader.rc_records() {
        let record = read_result.context("Failed getting read from reads file")?;

        if skip_read_check(&record) {
            continue;
        }

        let read_start = record.pos() as u64;
        let read_end = record.cigar().end_pos() as u64;

        while let Some(var) = bin.variants.front() {
            // Pop and write the variant features if the start of the read is > than the variant position
            if (read_start + 1) > var.info.pos {
                let pos_normalized = var.get_pos_normalized(ref_seq_len);
                write_variant_row(csv_writer, var, pos_normalized, sample)?;
                bin.variants.pop_front();
            } else {
                break;
            }
        }

        for var in &mut bin.variants {
            let zero_based_pos = var.info.pos - 1;

            if zero_based_pos >= read_start && zero_based_pos < read_end {
                var.count_locus_features(&record, zero_based_pos);
            } else {
                break;
            }
        }
    }

    // Pop and write any remaining variants
    while let Some(var) = bin.variants.pop_front() {
        let pos_normalized = var.get_pos_normalized(ref_seq_len);
        write_variant_row(csv_writer, &var, pos_normalized, sample)?;
    }

    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FileType {
    Bam,
    Cram,
}

pub fn detect_file_type(path: &Path) -> Result<FileType, AppError> {
    let ext = path.extension().and_then(|e| e.to_str()).ok_or_else(|| {
        AppError::MissingReadsExtension {
            filename: path.to_path_buf(),
        }
    })?;

    match ext {
        "bam" => Ok(FileType::Bam),
        "cram" => Ok(FileType::Cram),
        _ => Err(AppError::UnsupportedReadsFormat {
            filename: path.to_path_buf(),
            extension: ext.to_string(),
        }),
    }
}

#[derive(Debug, Clone)]
struct ChromInfo {
    min_pos: u64,
    max_pos: u64,
}

// Get the min/max of a given chromosome so that it can be passed into fetch
// NOTE: This assumes that the variants from one chromosome only
fn get_chrom_info(variants: &VecDeque<Variant>) -> ChromInfo {
    let first = variants
        .front()
        .expect("Internal error: chromosome contains no variants.
                 INVARIANT BROKEN: a bin with a given chromosome should ALWAYS contain at least one variant");

    let last = variants
        .back()
        .expect("Internal error: chromosome contains no variants.
                 INVARIANT BROKEN: a bin with a given chromosome should ALWAYS contain at least one variant");

    ChromInfo {
        min_pos: first.info.pos,
        max_pos: last.info.pos,
    }
}

fn get_ref_len(bam_reader: &IndexedReader, chrom: &str) -> Result<u64> {
    let header = bam_reader.header();

    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .context("BAM header contains an invalid reference name")?;
        if name == chrom {
            return header
                .target_len(tid)
                .with_context(|| format!("Reference '{}' found, but has no length", chrom));
        }
    }

    Err(AppError::MissingReferenceSequence {
        chromosome: chrom.to_owned(),
    }
    .into())
}

fn skip_read_check(read: &Record) -> bool {
    // Check if read is orphan pair as this is skipped in the origial varlap pileup call (ignore_orphans=True)
    if read.is_paired() && !read.is_proper_pair() {
        return true;
    }

    // Settings equivalent to stepper='samtools'?
    // See -ff at https://www.htslib.org/doc/samtools-mpileup.html#DESCRIPTION
    if read.is_unmapped()
        || read.is_secondary()
        || read.is_quality_check_failed()
        || read.is_duplicate()
    {
        return true;
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::{VariantInfo, VarType};
    use rust_htslib::bam::record::Record;
    use std::path::PathBuf;

    // Detect file type tests
    #[test]
    fn detect_file_type_cases() {
        let cases = [
            ("reads.bam", FileType::Bam),
            ("reads.cram", FileType::Cram),
        ];

        for (filename, expected) in cases {
            let path = PathBuf::from(filename);
            assert_eq!(detect_file_type(&path).unwrap(), expected, "filename={filename:?}");
        }
    }

    #[test]
    fn detect_file_type_missing_extension() {
        let path = PathBuf::from("reads");
        assert!(matches!(
            detect_file_type(&path),
            Err(AppError::MissingReadsExtension { .. })
        ));
    }

    #[test]
    fn detect_file_type_unsupported_extension() {
        let path = PathBuf::from("reads.sam");
        assert!(matches!(
            detect_file_type(&path),
            Err(AppError::UnsupportedReadsFormat { .. })
        ));
    }

    // Get chrom info tests
    fn variant(chrom: &str, pos: u64, vartype: VarType) -> VariantInfo {
        VariantInfo {
            chrom: chrom.to_string(),
            pos,
            vartype,
            refr: "A".into(),
            alt: "T".into(),
        }
    }

    #[test]
    fn get_chrom_info_cases() {
        let cases = [
            ("single variant", vec![100], (100, 100)),
            ("two variants, ascending", vec![100, 250], (100, 250)),
            ("many variants, only first/last matter", vec![100, 150, 175, 300], (100, 300)),
        ];

        for (label, positions, (expected_min, expected_max)) in cases {
            let infos: Vec<VariantInfo> = positions
                .into_iter()
                .map(|pos| variant("chr1", pos, VarType::Snv))
                .collect();

            let variants: VecDeque<Variant> = infos
                .iter()
                .map(crate::binning::make_variant_features)
                .collect();

            let info = get_chrom_info(&variants);

            assert_eq!(info.min_pos, expected_min, "case: {label}");
            assert_eq!(info.max_pos, expected_max, "case: {label}");
        }
    }

    #[test]
    #[should_panic(expected = "INVARIANT BROKEN")]
    fn get_chrom_info_panics_on_empty() {
        let variants: VecDeque<Variant> = VecDeque::new();
        let _ = get_chrom_info(&variants);
    }

    // Skip read check tests
    fn make_read(
        paired: bool,
        proper_pair: bool,
        unmapped: bool,
        secondary: bool,
        qc_fail: bool,
        duplicate: bool,
    ) -> Record {
        let mut record = Record::new();

        if paired {
            record.set_paired();
        } else {
            record.unset_paired();
        }

        if proper_pair {
            record.set_proper_pair();
        } else {
            record.unset_proper_pair();
        }

        if unmapped {
            record.set_unmapped();
        } else {
            record.unset_unmapped();
        }

        if secondary {
            record.set_secondary();
        } else {
            record.unset_secondary();
        }

        if qc_fail {
            record.set_quality_check_failed();
        } else {
            record.unset_quality_check_failed();
        }

        if duplicate {
            record.set_duplicate();
        } else {
            record.unset_duplicate();
        }

        record
    }

    #[test]
    fn skip_read_check_cases() {
        // (paired, proper_pair, unmapped, secondary, qc_fail, duplicate, expected_skip)
        let cases = [
            ("single-end, clean", (false, false, false, false, false, false), false),
            ("paired, proper pair, clean", (true, true, false, false, false, false), false),
            ("paired but not proper pair (orphan)", (true, false, false, false, false, false), true),
            ("unmapped", (false, false, true, false, false, false), true),
            ("secondary alignment", (false, false, false, true, false, false), true),
            ("failed QC", (false, false, false, false, true, false), true),
            ("duplicate", (false, false, false, false, false, true), true),
            ("proper pair but also unmapped", (true, true, true, false, false, false), true),
        ];

        for (label, (paired, proper_pair, unmapped, secondary, qc_fail, duplicate), expected) in cases
        {
            let read = make_read(paired, proper_pair, unmapped, secondary, qc_fail, duplicate);
            assert_eq!(skip_read_check(&read), expected, "case: {label}");
        }
    }
}
