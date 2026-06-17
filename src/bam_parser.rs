use std::collections::VecDeque;
use std::path::Path;
use std::rc::Rc;
use csv::{WriterBuilder};
use rust_htslib::bam::{Read, IndexedReader, Record};
use anyhow::{Context, Result, bail};

use crate::variant::{Variant, VarClass};
use crate::output::{write_variant_row, write_header};

// NOTE: For now the algorithm only parses inputs which have a single chromosome only
// This is temporary depending on how we want to multithread
pub fn parse_region(
    variants: &mut VecDeque<Variant>,
    reads_path: &Path,
    csv_path: &str,
    sample: Option<&str>,
    label: Option<&str>,
    varclass: &VarClass,
    fasta_path: Option<&Path>,
    // ADD vector of [chrom/first variant index/length] here
) -> Result<(), Box<dyn std::error::Error>> {
    let file_type = detect_file_type(reads_path)
        .with_context(|| format!("Failed to detect file type for {}", reads_path.display()))?;
    
    let mut csv_writer = WriterBuilder::new()
        .has_headers(false)
        .from_path(csv_path)?;

    // Write dynamic header based on BAM/CRAM filename
    write_header(&mut csv_writer, &reads_path, label, &varclass)?;

    // ADD LOOP FOR CHROM HERE 

    // Get min/max position of variants in a given chromosome 
    //  to fetch only reads that are in this region
    // NOTE: This assumes that the variants are of only one chromosome
    let chrom_info = get_chrom_info(&variants)
        .ok_or("Could not determine chromosome min/max")?;

    let mut reader = IndexedReader::from_path(reads_path)?;

    match file_type {
        FileType::Bam => {}
        FileType::Cram => {
            if let Some(path) = fasta_path {
                reader
                    .set_reference(path)
                    .with_context(|| format!("Failed to set CRAM reference to: {}", path.display()))?;
            } else {
                return Err("A reference FASTA is required to read CRAM files".into());
            }
        },
    };

    reader.fetch((&chrom_info.chrom, chrom_info.min_pos - 1, chrom_info.max_pos))?;

    let ref_seq_len = get_ref_len(&reader, &chrom_info.chrom)?;

    for read_result in reader.rc_records() {
        let record = read_result?;
        
        if skip_read_check(&record) {
            continue;
        }

        let read_start = record.pos() as u64;

        loop {
            let should_pop = match variants.front() {
                Some(var) => (read_start + 1) > var.pos,
                None => false,
            };

            if should_pop {
                if let Some(var) = variants.pop_front() {
                    let pos_fraction = var.get_pos_fraction(ref_seq_len);
                    write_variant_row(&mut csv_writer, &var, pos_fraction, sample)?;
                }
            } else {
                break;
            }
        }

        for var in &mut *variants {
            let zero_based_pos = var.pos - 1;
            let read_end = record.cigar().end_pos() as u64;

            if zero_based_pos >= read_start && zero_based_pos < read_end {
                var.count_locus_features(&record, zero_based_pos);
            } else {
                break;
            }
        }
    }

    while let Some(var) = variants.pop_front() {
        let pos_fraction = var.get_pos_fraction(ref_seq_len);
        write_variant_row(&mut csv_writer, &var, pos_fraction, sample)?;
    }
    csv_writer.flush()?;

    Ok(()) 
}

#[derive(Debug, Clone, Copy)]
enum FileType {
    Bam,
    Cram,
}

fn detect_file_type(path: &Path) -> Result<FileType> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .context("Missing file extension")?;

    match ext {
        "bam" => Ok(FileType::Bam),
        "cram" => Ok(FileType::Cram),
        _ => bail!("Unsupported file format: {}", ext),
    }
}

#[derive(Debug, Clone)]
struct ChromInfo {
    chrom: String,
    min_pos: u64,
    max_pos: u64,
}

// NOTE: This assumes that the variants from only one chromosome
fn get_chrom_info(variants: &VecDeque<Variant>) -> Option<ChromInfo> {
    let first = variants.front()?;

    Some(ChromInfo {
        chrom: first.chrom.clone(),
        min_pos: first.pos,
        max_pos: variants.back()?.pos,
    })
}

fn get_ref_len(
    bam_reader: &IndexedReader,
    chrom: &str,
) -> Result<u64, Box<dyn std::error::Error>> {
    let header = bam_reader.header();

    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))?;
        if name == chrom {
            return header
                .target_len(tid)
                .ok_or_else(|| format!("Reference '{}' found, but has no length", chrom).into());
        }
    }

    Err(format!("Could not find reference '{}' in BAM header", chrom).into())
}

fn skip_read_check(read: &Rc<Record>) -> bool {
    // Check if read is orphan pair as this is skipped in the origial varlap pileup call (ignore_orphans=True)
    if read.is_paired() && !read.is_proper_pair() {
        return true;
    }

    // Settings equivalent to stepper='samtools'? 
    // See -ff at https://www.htslib.org/doc/samtools-mpileup.html#DESCRIPTION
    if read.is_unmapped() 
        || read.is_secondary() 
        || read.is_quality_check_failed() 
        || read.is_duplicate() {
        return true;
    }

    false
}
