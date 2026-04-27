use std::collections::VecDeque;
use std::rc::Rc;
use csv::Writer;
use rust_htslib::bam::{Read, IndexedReader, Record};

mod cli;
mod vcf_parser;
mod variant;
mod features;
mod output;

use crate::variant::{Variant, VarClass};
use crate::output::{write_variant_row, CSV_HEADER_SNV, CSV_HEADER_INDEL};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let vcf_path = "test_data/chrM_heavy_stress.vcf";
    // let bam_path: &str = "test_data/chrM_heavy_stress.sorted.bam";
    // let varclass = "SNV";
    // let csv_path = "test_data/rust-varlap.output.csv";
    // let sample = "";

    let args = cli::parse();
    
    let mut variants = vcf_parser::parse(&args.vcf, &args.varclass)?;

    let (region_chrom, min_pos, max_pos) = 
        get_vcf_min_max(&variants).ok_or("Could not determine VCF min/max")?;

    println!("Region Chromosome: {}, Min Pos: {}, Max Pos: {}", &region_chrom, min_pos, max_pos);

    // VARCLASS INPUT TEMP FIX FOR CSV HEADER
    process_bam_region(&mut variants, &args.bams, &region_chrom, min_pos, max_pos, &args.csv_path, args.sample.as_deref(), &args.varclass)?;

    Ok(())
}

fn get_vcf_min_max(variants: &VecDeque<Variant>) -> Option<(String, u64, u64)> {
    let first = variants.front()?;
    let chrom = first.chrom.clone();

    let min_pos = first.pos;
    let max_pos = variants.back()?.pos;

    Some((chrom, min_pos, max_pos))
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

fn process_bam_region(
    variants: &mut VecDeque<Variant>,
    bam_path: &str,
    region_chrom: &str, 
    min_pos: u64, 
    max_pos: u64,
    csv_path: &str,
    sample: Option<&str>,
    // FOR TEMP HEADER FIX
    varclass: &VarClass,
    //
) -> Result<(), Box<dyn std::error::Error>> {
    let mut bam_reader = IndexedReader::from_path(bam_path)?;
    bam_reader.fetch((region_chrom, min_pos - 1, max_pos))?;

    let mut csv_writer = Writer::from_path(csv_path)?;

    // TEMP HEADER CSV FIX
    match varclass {
        VarClass::Snv =>  csv_writer.write_record(CSV_HEADER_SNV)?,
        VarClass::Indel =>  csv_writer.write_record(CSV_HEADER_INDEL)?,
    }
    //

    let ref_seq_len = get_ref_len(&bam_reader, &region_chrom)?;

    for read_result in bam_reader.rc_records() {
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