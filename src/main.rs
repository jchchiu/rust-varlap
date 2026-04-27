use std::fs::File;
use std::cmp;
use std::collections::VecDeque;
use std::error::Error;
use std::rc::Rc;
use csv::Writer;
use serde::Serialize;
use rust_htslib::bam::{Read, IndexedReader, Record};
use rust_htslib::bam::record::{Aux, Cigar};

mod cli;
mod vcf_parser;
mod variant;

use crate::variant::{Variant, VarClass, VarType};

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

// enum LocusFeatures and match?
#[derive(Debug, Clone)]
enum LocusFeatures {
    Snv(LocusFeaturesSNV),
    Indel(LocusFeaturesINDEL),
}

impl LocusFeatures {
    fn normalized_row(&self) -> NormalizedLocusFeaturesRow {
        match self {
            LocusFeatures::Snv(f) => f.common.normalized_row(),
            LocusFeatures::Indel(f) => f.common.normalized_row(),
        }
    }
}

#[derive(Debug, Clone, Default)]
struct CommonLocusFeatures {
    ref_read_features: ReadFeatures,
    alt_read_features: ReadFeatures,
    all_read_features: ReadFeatures,
}

impl CommonLocusFeatures {
    fn normalized_row(&self) -> NormalizedLocusFeaturesRow {
        let r = self.ref_read_features.normalized();
        let a = self.alt_read_features.normalized();
        let all = self.all_read_features.normalized();

        NormalizedLocusFeaturesRow {
            ref_nm: r.nm,
            ref_base_qual: r.base_qual,
            ref_map_qual: r.map_qual,
            ref_align_len: r.align_len,
            ref_clipping: r.clipping,
            ref_indel: r.indel,
            ref_forward_strand: r.forward_strand,
            ref_reverse_strand: r.reverse_strand,
            ref_supplementary: r.supplementary,
            ref_normalised_read_position: r.normalised_read_position,

            alt_nm: a.nm,
            alt_base_qual: a.base_qual,
            alt_map_qual: a.map_qual,
            alt_align_len: a.align_len,
            alt_clipping: a.clipping,
            alt_indel: a.indel,
            alt_forward_strand: a.forward_strand,
            alt_reverse_strand: a.reverse_strand,
            alt_supplementary: a.supplementary,
            alt_normalised_read_position: a.normalised_read_position,

            all_nm: all.nm,
            all_base_qual: all.base_qual,
            all_map_qual: all.map_qual,
            all_align_len: all.align_len,
            all_clipping: all.clipping,
            all_indel: all.indel,
            all_forward_strand: all.forward_strand,
            all_reverse_strand: all.reverse_strand,
            all_supplementary: all.supplementary,
            all_normalised_read_position: all.normalised_read_position,
        }
    }    
}

#[derive(Debug, Clone, Default)]
struct LocusFeaturesSNV {
    base_counts: BaseCountsSNV,
    common: CommonLocusFeatures,
}

impl LocusFeaturesSNV {
    fn count(
        &mut self,
        read: &Rc<Record>,
        base: Option<u8>,
        refr: char,
        alt: char,
        query_pos: Option<u32>,
    ) {
        if let Some(base_u8) = base {
            let base_char = base_u8 as char;

            self.base_counts.count(base_char);

            if base_char == refr {
                self.common.ref_read_features.count(read, query_pos);
            } else if base_char == alt {
                self.common.alt_read_features.count(read, query_pos);
            }
        }

        self.common.all_read_features.count(&read, query_pos);
    }
}

#[derive(Debug, Clone)]
struct IndelEvent {
    indel_type: VarType,
    start: u64,
    end: u64,
    bases: String,
}

#[derive(Debug, Clone, Default)]
struct LocusFeaturesINDEL {
    overlapping_indels_count: u64,
    common: CommonLocusFeatures,
}

impl LocusFeaturesINDEL {
    fn count(
        &mut self,
        read: &Rc<Record>,
        refr: &str,
        alt: &str,
        ref_pos: u64,
        query_pos: Option<u32>,
        indel_type: &VarType,
    ) {
        let start = self.get_indel_start_coord(ref_pos, refr, alt);
        let size = refr.len().abs_diff(alt.len()) as u64;
        let end = start + size - 1;

        let overlapping_indels = self.indels_overlapping_variant(&read, start, end);
        self.overlapping_indels_count += overlapping_indels.len() as u64;

        let mut read_supports_alt = false;

        let bases = match indel_type {
            VarType::Del => Some(String::new()),
            VarType::Ins => Some(alt[1..].to_string()),
            _ => None,
        };

        for event in overlapping_indels.iter() {
            if event.indel_type == *indel_type {
                if (event.start == start) && (event.end == end) {
                    if matches!(indel_type, VarType::Del) ||           // FIX UNWRAP HERE
                    (matches!(indel_type, VarType::Ins) && event.bases == bases.clone().unwrap()) {
                        read_supports_alt = true;
                        break;
                    }
                }
            }
        }

        let mut read_supports_ref = false;
        if overlapping_indels.is_empty() && let Some(qpos) = query_pos {
            let seq_bytes = read.seq().as_bytes();
            let read_bases = match indel_type {
                VarType::Ins => {
                    Some((seq_bytes[qpos as usize] as char)
                    .to_string()
                    .to_ascii_uppercase())
                },
                VarType::Del => {
                    Some(String::from_utf8(
                        seq_bytes[qpos as usize .. (qpos + (size as u32) + 1) as usize].to_vec()
                    )
                    .unwrap()
                    .to_ascii_uppercase())
                },
                _ => None,
            };
            
            // FIX UNWRAP HERE
            if refr == read_bases.unwrap() {
                read_supports_ref = true;
            }
        }

        if read_supports_ref {
            self.common.ref_read_features.count(read, query_pos);
        } else if read_supports_alt {
            self.common.alt_read_features.count(read, query_pos);
        }

        self.common.all_read_features.count(&read, query_pos);
    }

    // get the genome coordinates of where an INDEL variant will actually
    // start, as opposed to the location of where the variant is reported
    // the starting position must take into account the context bases that
    // are given when the variant is reported in the VCF file.
    // Note: if we assume normalized biallelic variants, could just do (pos + 1)
    fn get_indel_start_coord(&self, pos: u64, refr: &str, alt: &str) -> u64 {
        let shortest_len = cmp::min(refr.len(), alt.len());
        pos + shortest_len as u64
    }

    // True if the intervals of two indels overlap
    fn interval_overlaps(&self, start1: u64, end1: u64, start2: u64, end2: u64) -> bool {
        !((end1 < start2) || (start1 > end2))
    }

    // Determine the allele in the read at the locus of an INDEL variant
    fn indels_overlapping_variant(
        &self,
        read: &Rc<Record>,
        var_start: u64,
        var_end: u64,
    ) -> Vec<IndelEvent> {
        let mut read_pos: u32 = 0;
        let mut ref_pos = read.pos() as u32;
        let mut result = Vec::new();

        // See https://samtools.github.io/hts-specs/SAMv1.pdf page 8 for how CIGAR consumes
        for c in read.cigar().iter() {
            match *c {
                // Consumes both reference and query
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    ref_pos += len;
                    read_pos += len;
                }
                // Only consumes query
                Cigar::Ins(len) => {
                    let this_start = ref_pos as u64;
                    let this_end = this_start + len as u64 - 1;

                    if self.interval_overlaps(var_start, var_end, this_start, this_end) {
                        let seq_bytes = read.seq().as_bytes();
                        let inserted_bases = String::from_utf8(
                            seq_bytes[read_pos as usize .. (read_pos + len) as usize].to_vec()
                        )
                        .unwrap()
                        .to_ascii_uppercase();

                        result.push(IndelEvent {
                            indel_type: VarType::Ins,
                            start: this_start,
                            end: this_end,
                            bases: inserted_bases,
                        });
                    }

                    read_pos += len;
                }
                // Only consumes reference
                Cigar::Del(len) => {
                    let this_start = ref_pos as u64;
                    let this_end = this_start + len as u64 - 1;

                    if self.interval_overlaps(var_start, var_end, this_start, this_end) {
                        result.push(IndelEvent {
                            indel_type: VarType::Del,
                            start: this_start,
                            end: this_end,
                            bases: String::new(),
                        });
                    }

                    ref_pos += len;
                }
                // Only consumes reference
                Cigar::RefSkip(len) => {
                    ref_pos += len;
                }
                // Only consumes query
                Cigar::SoftClip(len) => {
                    read_pos += len;
                }
                // Consumes neither reference nor query
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }        

        result
    }

}

impl LocusFeaturesINDEL {
    fn stats(&self) -> INDELCountsStats {
        let depth = self.common.all_read_features.num_reads;
        let ref_count = self.common.ref_read_features.num_reads;
        let alt_count = self.common.alt_read_features.num_reads;
        let other_count = depth - ref_count - alt_count;
        let alt_vaf = if depth > 0 {
            alt_count as f64 / depth as f64
        } else {
            0.0
        };

        INDELCountsStats {
            depth,
            ref_count,
            alt_count,
            other_count,
            alt_vaf,
            overlapping_indels_count: self.overlapping_indels_count,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
struct INDELCountsStats {
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    other_count: u32,
    alt_vaf: f64,
    overlapping_indels_count: u64,
}

#[derive(Debug, Clone, Default)]
struct BaseCountsSNV {
    a: u32,
    c: u32,
    g: u32,
    t: u32,
    n: u32,
}

impl BaseCountsSNV {
    fn count(&mut self, base: char) {
        match base {
            'A' => self.a += 1,
            'C' => self.c += 1,
            'G' => self.g += 1,
            'T' => self.t += 1,
            'N' => self.n += 1,
            _ => eprintln!("Warning: Base does not match: {}", base),
        }
    }

    fn count_for_base(&self, base: char) -> u32 {
        match base {
            'A' => self.a,
            'C' => self.c,
            'G' => self.g,
            'T' => self.t,
            'N' => self.n,
            _ => 0,
        }
    }

    // fn depth(&self) -> u32 {
    //     self.a + self.c + self.g + self.t + self.n
    // }

    fn stats(&self, refr: char, alt: char) -> BaseCountsSNVStats {
        let depth = self.a + self.c + self.g + self.t + self.n;
        let ref_count = self.count_for_base(refr);
        let alt_count = self.count_for_base(alt);
        let alt_vaf = if depth > 0 {
            alt_count as f64 / depth as f64
        } else {
            0.0
        };

        BaseCountsSNVStats {
            depth,
            ref_count,
            alt_count,
            alt_vaf,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
struct BaseCountsSNVStats {
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    alt_vaf: f64,
}

#[derive(Debug, Clone, Default)]
struct ReadFeatures {
    nm: u32,
    base_qual: u32,
    map_qual: u32,
    align_len: u32,
    clipping: u32,
    indel: u32,
    forward_strand: u32,
    reverse_strand: u32,
    supplementary: u32,
    normalised_read_position: f64, 
    num_reads: u32,
}

impl ReadFeatures {
    fn count(
        &mut self,
        read: &Rc<Record>,
        query_pos: Option<u32>,
    ) {
        // Instead of counting num of reads, can create a function that sums foward and reverse strand?
        self.num_reads += 1;
        let query_len = read.seq_len() as usize;

        if let Some(qpos) = query_pos {
            if query_len > 0 {
                self.normalised_read_position += qpos as f64 / query_len as f64;
            }

            let qpos_usize = qpos as usize;
            if qpos_usize < read.qual().len() {
                let pos_qual = read.qual()[qpos_usize];
                self.base_qual += pos_qual as u32;
            }
        }

        self.align_len += self.query_alignment_length(&read) as u32;
        self.map_qual += read.mapq() as u32;

        for c in read.cigar().iter() {
            match *c {
                Cigar::Ins(len) | Cigar::Del(len) => self.indel += len as u32,
                Cigar::SoftClip(len) | Cigar::HardClip(len) => self.clipping += len as u32,
                _ => {}
            }
        }

        if let Ok(aux) = read.aux(b"NM") {
            match aux {
                Aux::I8(v) => self.nm += v as u32,
                Aux::U8(v) => self.nm += v as u32,
                Aux::I16(v) => self.nm += v as u32,
                Aux::U16(v) => self.nm += v as u32,
                Aux::I32(v) => self.nm += v as u32,
                Aux::U32(v) => self.nm += v as u32,
                _ => {}
            }
        }

        if read.is_reverse() {
            self.reverse_strand += 1;
        } else {
            self.forward_strand += 1;
        }
        if read.is_supplementary() {
            self.supplementary += 1;
        }

    }

    fn query_alignment_length(&self, record: &Rc<Record>) -> u32 {
        let mut len = 0;
        for c in record.cigar().iter() {
            match *c {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Ins(l) => len += l,
                _ => {}
            }
        }
        len
    }

    fn normalized(&self) -> NormalizedReadFeatures {
        if self.num_reads > 0 {
            let n = self.num_reads as f64;
            NormalizedReadFeatures {
                nm: Some(self.nm as f64 / n),
                base_qual: Some(self.base_qual as f64 / n),
                map_qual: Some(self.map_qual as f64 / n),
                align_len: Some(self.align_len as f64 / n),
                clipping: Some(self.clipping as f64 / n),
                indel: Some(self.indel as f64 / n),
                forward_strand: Some(self.forward_strand as f64 / n),
                reverse_strand: Some(self.reverse_strand as f64 / n),
                supplementary: Some(self.supplementary as f64 / n),
                normalised_read_position: Some(self.normalised_read_position / n),
            }
        } else {
            NormalizedReadFeatures::default()
        }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize)]
struct NormalizedLocusFeaturesRow {
    ref_nm: Option<f64>,
    ref_base_qual: Option<f64>,
    ref_map_qual: Option<f64>,
    ref_align_len: Option<f64>,
    ref_clipping: Option<f64>,
    ref_indel: Option<f64>,
    ref_forward_strand: Option<f64>,
    ref_reverse_strand: Option<f64>,
    ref_supplementary: Option<f64>,
    ref_normalised_read_position: Option<f64>,

    alt_nm: Option<f64>,
    alt_base_qual: Option<f64>,
    alt_map_qual: Option<f64>,
    alt_align_len: Option<f64>,
    alt_clipping: Option<f64>,
    alt_indel: Option<f64>,
    alt_forward_strand: Option<f64>,
    alt_reverse_strand: Option<f64>,
    alt_supplementary: Option<f64>,
    alt_normalised_read_position: Option<f64>,

    all_nm: Option<f64>,
    all_base_qual: Option<f64>,
    all_map_qual: Option<f64>,
    all_align_len: Option<f64>,
    all_clipping: Option<f64>,
    all_indel: Option<f64>,
    all_forward_strand: Option<f64>,
    all_reverse_strand: Option<f64>,
    all_supplementary: Option<f64>,
    all_normalised_read_position: Option<f64>,
}

// Use option or not? Because in python script if there are no values it defaults to ''
// Compared to default which would initialize it as 0.0
#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
struct NormalizedReadFeatures {
    nm: Option<f64>,
    base_qual: Option<f64>,
    map_qual: Option<f64>,
    align_len: Option<f64>,
    clipping: Option<f64>,
    indel: Option<f64>,
    forward_strand: Option<f64>,
    reverse_strand: Option<f64>,
    supplementary: Option<f64>,
    normalised_read_position: Option<f64>,
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

// NOTE: SERDE NO LONGER WRITING HEADER BECAUSE OF NESTED: THINK OF FIX

fn write_variant_row(
    writer: &mut Writer<File>,
    var: &Variant,
    pos_fraction: f64,
    sample: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    writer.serialize(OutputRow::from_variant(var, pos_fraction, sample))?;
    Ok(())
}

#[derive(serde::Serialize)]
struct OutputReadFeatures {
    ref_nm: Option<f64>,
    ref_base_qual: Option<f64>,
    ref_map_qual: Option<f64>,
    ref_align_len: Option<f64>,
    ref_clipping: Option<f64>,
    ref_indel: Option<f64>,
    ref_forward_strand: Option<f64>,
    ref_reverse_strand: Option<f64>,
    ref_supplementary: Option<f64>,
    ref_normalised_read_position: Option<f64>,

    alt_nm: Option<f64>,
    alt_base_qual: Option<f64>,
    alt_map_qual: Option<f64>,
    alt_align_len: Option<f64>,
    alt_clipping: Option<f64>,
    alt_indel: Option<f64>,
    alt_forward_strand: Option<f64>,
    alt_reverse_strand: Option<f64>,
    alt_supplementary: Option<f64>,
    alt_normalised_read_position: Option<f64>,

    all_nm: Option<f64>,
    all_base_qual: Option<f64>,
    all_map_qual: Option<f64>,
    all_align_len: Option<f64>,
    all_clipping: Option<f64>,
    all_indel: Option<f64>,
    all_forward_strand: Option<f64>,
    all_reverse_strand: Option<f64>,
    all_supplementary: Option<f64>,
    all_normalised_read_position: Option<f64>,
}

impl From<NormalizedLocusFeaturesRow> for OutputReadFeatures {
    fn from(rf: NormalizedLocusFeaturesRow) -> Self {
        Self {
            ref_nm: rf.ref_nm,
            ref_base_qual: rf.ref_base_qual,
            ref_map_qual: rf.ref_map_qual,
            ref_align_len: rf.ref_align_len,
            ref_clipping: rf.ref_clipping,
            ref_indel: rf.ref_indel,
            ref_forward_strand: rf.ref_forward_strand,
            ref_reverse_strand: rf.ref_reverse_strand,
            ref_supplementary: rf.ref_supplementary,
            ref_normalised_read_position: rf.ref_normalised_read_position,

            alt_nm: rf.alt_nm,
            alt_base_qual: rf.alt_base_qual,
            alt_map_qual: rf.alt_map_qual,
            alt_align_len: rf.alt_align_len,
            alt_clipping: rf.alt_clipping,
            alt_indel: rf.alt_indel,
            alt_forward_strand: rf.alt_forward_strand,
            alt_reverse_strand: rf.alt_reverse_strand,
            alt_supplementary: rf.alt_supplementary,
            alt_normalised_read_position: rf.alt_normalised_read_position,

            all_nm: rf.all_nm,
            all_base_qual: rf.all_base_qual,
            all_map_qual: rf.all_map_qual,
            all_align_len: rf.all_align_len,
            all_clipping: rf.all_clipping,
            all_indel: rf.all_indel,
            all_forward_strand: rf.all_forward_strand,
            all_reverse_strand: rf.all_reverse_strand,
            all_supplementary: rf.all_supplementary,
            all_normalised_read_position: rf.all_normalised_read_position,
        }
    }
}

#[derive(serde::Serialize)]
enum OutputRow<'a> {
    Snv(OutputRowSNV<'a>),
    Indel(OutputRowINDEL<'a>),
}

impl<'a> OutputRow<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64, sample: Option<&'a str>) -> Self {
        match &var.features {
            LocusFeatures::Snv(_) => OutputRow::Snv(OutputRowSNV::from_variant(var, pos_fraction, sample)),
            LocusFeatures::Indel(_) => OutputRow::Indel(OutputRowINDEL::from_variant(var, pos_fraction, sample)),
        }
    }
}

#[derive(serde::Serialize)]
struct OutputRowSNV<'a> {
    chrom: &'a str,
    pos: u64,
    #[serde(rename = "ref")]
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalised: f64,
    sample: Option<&'a str>,
    depth: u32,

    count_a: u32,
    count_t: u32,
    count_g: u32,
    count_c: u32,
    count_n: u32,

    ref_count: u32,
    alt_count: u32,
    alt_vaf: f64,
    
    read_features: OutputReadFeatures,
}

// Make branching for snv/indel; make common function for read_features
impl<'a> OutputRowSNV<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64, sample: Option<&'a str>) -> Self {
        let bcs = var.base_counts_stats().expect("Could not get base count statistics.");
        let read_features = OutputReadFeatures::from(var.normalized_row());

        let f = match &var.features {
            LocusFeatures::Snv(f) => f,
            LocusFeatures::Indel(_) => panic!("Expected SNV features"),
        };

        Self {
            chrom: &var.chrom,
            pos: var.pos,
            refr: &var.refr,
            alt: &var.alt,
            vartype: &var.vartype.as_str(),
            pos_normalised: pos_fraction,
            sample: sample,
            depth: bcs.depth,
            
            count_a: f.base_counts.a,
            count_t: f.base_counts.t,
            count_g: f.base_counts.g,
            count_c: f.base_counts.c,
            count_n: f.base_counts.n,

            ref_count: bcs.ref_count,
            alt_count: bcs.alt_count,
            alt_vaf: bcs.alt_vaf,

            read_features,
        }
    }
}

#[derive(serde::Serialize)]
struct OutputRowINDEL<'a> {
    chrom: &'a str,
    pos: u64,
    #[serde(rename = "ref")]
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalised: f64,
    sample: Option<&'a str>,
    depth: u32,

    ref_count: u32,
    alt_count: u32,
    other_count: u32,
    alt_vaf: f64,
    overlapping_indels_count: u64,

    read_features: OutputReadFeatures,
}

// Make branching for snv/indel; make common function for read_features
impl<'a> OutputRowINDEL<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64, sample: Option<&'a str>) -> Self {
        let read_features = OutputReadFeatures::from(var.normalized_row());

        let stats = var.indel_stats().expect("Could not get indel count statistics.");

        Self {
            chrom: &var.chrom,
            pos: var.pos,
            refr: &var.refr,
            alt: &var.alt,
            vartype: &var.vartype.as_str(),
            pos_normalised: pos_fraction,
            sample: sample,
            depth: stats.depth,

            ref_count: stats.ref_count,
            alt_count: stats.alt_count,
            other_count: stats.other_count,
            alt_vaf: stats.alt_vaf,
            overlapping_indels_count: stats.overlapping_indels_count,

            read_features,
        }
    }
}

// TEMP HEADER FIX
const CSV_HEADER_SNV: &[&str] = &[
    "chrom",
    "pos",
    "ref",
    "alt",
    "vartype",
    "pos_normalised",
    "sample",
    "depth",
    "count_a",
    "count_t",
    "count_g",
    "count_c",
    "count_n",
    "ref_count",
    "alt_count",
    "alt_vaf",
    "ref_nm",
    "ref_base_qual",
    "ref_map_qual",
    "ref_align_len",
    "ref_clipping",
    "ref_indel",
    "ref_forward_strand",
    "ref_reverse_strand",
    "ref_supplementary",
    "ref_normalised_read_position",
    "alt_nm",
    "alt_base_qual",
    "alt_map_qual",
    "alt_align_len",
    "alt_clipping",
    "alt_indel",
    "alt_forward_strand",
    "alt_reverse_strand",
    "alt_supplementary",
    "alt_normalised_read_position",
    "all_nm",
    "all_base_qual",
    "all_map_qual",
    "all_align_len",
    "all_clipping",
    "all_indel",
    "all_forward_strand",
    "all_reverse_strand",
    "all_supplementary",
    "all_normalised_read_position",
];

const CSV_HEADER_INDEL: &[&str] = &[
    "chrom",
    "pos",
    "ref",
    "alt",
    "vartype",
    "pos_normalised",
    "sample",
    "depth",
    "ref_count",
    "alt_count",
    "alt_vaf",
    "other_count",
    "overlapping_indels_count",
    "ref_nm",
    "ref_base_qual",
    "ref_map_qual",
    "ref_align_len",
    "ref_clipping",
    "ref_indel",
    "ref_forward_strand",
    "ref_reverse_strand",
    "ref_supplementary",
    "ref_normalised_read_position",
    "alt_nm",
    "alt_base_qual",
    "alt_map_qual",
    "alt_align_len",
    "alt_clipping",
    "alt_indel",
    "alt_forward_strand",
    "alt_reverse_strand",
    "alt_supplementary",
    "alt_normalised_read_position",
    "all_nm",
    "all_base_qual",
    "all_map_qual",
    "all_align_len",
    "all_clipping",
    "all_indel",
    "all_forward_strand",
    "all_reverse_strand",
    "all_supplementary",
    "all_normalised_read_position",
];
//