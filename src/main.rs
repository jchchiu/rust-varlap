use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::VecDeque;
use std::error::Error;
use std::rc::Rc;
use csv::Writer;
use serde::Serialize;
use rust_htslib::bam::{Read, IndexedReader, Record};
use rust_htslib::bam::record::{Aux, Cigar};

#[derive(serde::Serialize)]
struct OutputRow<'a> {
    chrom: &'a str,
    pos: u64,
    #[serde(rename = "ref")]
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalised: f64,

    count_a: u32,
    count_c: u32,
    count_g: u32,
    count_t: u32,
    count_n: u32,

    #[serde(flatten)]
    base_counts_stats: BaseCountsStats,

    #[serde(flatten)]
    read_features: NormalizedLocusFeaturesRow,
}

impl<'a> OutputRow<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64) -> Self {
        Self {
            chrom: &var.chrom,
            pos: var.pos,
            refr: &var.refr,
            alt: &var.alt,
            vartype: &var.vartype.as_str(),
            pos_normalised: pos_fraction,

            count_a: var.counts.a,
            count_c: var.counts.c,
            count_g: var.counts.g,
            count_t: var.counts.t,
            count_n: var.counts.n,

            base_counts_stats: var.base_counts_stats().expect("Could not get base count statistics."),

            read_features: var.read_features.normalized_row(),
        }
    }
}

fn write_header_row(writer: &mut Writer<File>) -> Result<(), Box<dyn Error>> {
    writer.write_record(&[
        "chrom", "pos", "ref", "alt", "vartype", "pos_normalised", "depth",
        "A", "T", "G", "C", "N", "ref_count", "alt_count", "alt_vaf",
    ])?;
    Ok(())
}

fn write_variant_row(
    writer: &mut Writer<File>,
    var: &Variant,
    base_stats: &BaseCountsStats,
    pos_fraction: f64,
) -> Result<(), Box<dyn Error>> {
    writer.write_record(&[
        &var.chrom,
        &var.pos.to_string(),
        &var.refr,
        &var.alt,
        &var.vartype.as_str().to_string(),
        &pos_fraction.to_string(),
        &base_stats.depth.to_string(),
        &var.counts.a.to_string(),
        &var.counts.t.to_string(),
        &var.counts.g.to_string(),
        &var.counts.c.to_string(),
        &var.counts.n.to_string(),
        &base_stats.refr_count.to_string(),
        &base_stats.alt_count.to_string(),
        &base_stats.alt_vaf.to_string(),
    ])?;
    Ok(())
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

fn process_bam_region(
    variants: &mut VecDeque<Variant>,
    bam_path: &str,
    region_chrom: &str, 
    min_pos: u64, 
    max_pos: u64,
    csv_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut bam_reader = IndexedReader::from_path(bam_path)?;
    bam_reader.fetch((region_chrom, min_pos - 1, max_pos))?;

    let mut csv_writer = Writer::from_path(csv_path)?;
    write_header_row(&mut csv_writer)?;

    let ref_seq_len = get_ref_len(&bam_reader, &region_chrom)?;

    for read_result in bam_reader.rc_records() {
        let record = read_result?;
        let read_start = record.pos() as u64;
        let seq = record.seq();

        loop {
            let should_pop = match variants.front() {
                Some(var) => (read_start + 1) > var.pos,
                None => false,
            };

            if should_pop {
                if let Some(var) = variants.pop_front() {
                    let base_counts_stats = var.base_counts_stats().ok_or("error")?;
                    let pos_fraction = var.get_pos_fraction(ref_seq_len);
                    write_variant_row(&mut csv_writer, &var, &base_counts_stats, pos_fraction)?;
                }
            } else {
                break;
            }
        }

        for var in &mut *variants {
            let zero_based_pos = var.pos - 1;
            let read_end = read_start + record.seq_len() as u64;

            if zero_based_pos >= read_start && zero_based_pos < read_end {
                let base = seq[(zero_based_pos - read_start) as usize] as char;
                var.counts.count(base);
                let refr_char = var.refr.chars().next().ok_or("Could not get refr char")?;
                let alt_char = var.alt.chars().next().ok_or("Could not get alt char")?;
                var.read_features.count(&record, base, refr_char, alt_char, zero_based_pos);
            } else {
                break;
            }
        }
    }

    while let Some(var) = variants.pop_front() {
        let base_counts_stats = var.base_counts_stats().ok_or("error")?;
        let pos_fraction = var.get_pos_fraction(ref_seq_len);
        write_variant_row(&mut csv_writer, &var, &base_counts_stats, pos_fraction)?;
    }
    csv_writer.flush()?;

    Ok(()) 
}

fn varclass_matches(varclass: &str, vartype: &VarType) -> bool {
    match varclass.to_ascii_uppercase().as_str() {
        "SNV" => matches!(vartype, VarType::Snv | VarType::Unknown),
        "INDEL" => matches!(vartype, VarType::Ins | VarType::Del | VarType::Unknown),
        _ => false,
    }
}

fn vcf_reader(file_path: &str, varclass: &str) -> Result<VecDeque<Variant>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

	let mut variants = VecDeque::new();

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("#") {
            continue;
            }
        
        let fields: Vec<&str> = line.split_whitespace().collect();
        
        if fields.len() >= 5 {
            let chrom = fields[0].to_string();

            let pos = match fields[1].parse::<u64>() {
                Ok(p) => p,
                Err(_) => {
                    eprintln!("Warning: invalid POS, skipping row: {}", line);
                    continue;
                }
            };

            let refr = fields[3].to_string();

            for alt in fields[4].split(',') {
                let vartype = get_var_type(&refr, alt);

                if !varclass_matches(varclass, &vartype) {
                    continue;
                }

                variants.push_back(Variant {
                    chrom: chrom.clone(),
                    pos,
                    refr: refr.clone(),
                    alt: alt.to_string(),
                    vartype,
                    counts: BaseCounts::default(),
                    read_features: LocusFeaturesSNV::default(),
                });
            }
        } else {
            eprintln!("Warning: Skipping input row: {}", line);
        }
    }
			
	Ok(variants)
}

#[derive(Debug, Clone, Copy, Serialize)]
struct BaseCountsStats {
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    alt_vaf: f64,
}

#[derive(Debug, Clone, Default)]
struct BaseCounts {
    a: u32,
    c: u32,
    g: u32,
    t: u32,
    n: u32,
}

impl BaseCounts {
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

    fn depth(&self) -> u32 {
        self.a + self.c + self.g + self.t + self.n
    }

    fn stats(&self, refr: char, alt: char) -> BaseCountsStats {
        let depth = self.depth();
        let ref_count = self.count_for_base(refr);
        let alt_count = self.count_for_base(alt);
        let alt_vaf = if depth > 0 {
            alt_count as f64 / depth as f64
        } else {
            0.0
        };

        BaseCountsStats {
            depth,
            ref_count,
            alt_count,
            alt_vaf,
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
        query_pos: u64,
    ) {
        self.num_reads += 1;
        let query_len = read.seq_len();
        if query_len > 0 {
            self.normalised_read_position += query_pos as f64 / query_len as f64;
        }
        // Can maybe to Some(pos_qual) = read.qual().get() if want to return an Option for safety
        let pos_qual = read.qual()[query_pos as usize];
        self.base_qual += pos_qual as u32;

        self.align_len += query_len as u32;
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

#[derive(Debug, Clone, Default)]
struct LocusFeaturesSNV {
    ref_read_features: ReadFeatures,
    alt_read_features: ReadFeatures,
    all_read_features: ReadFeatures,
}

impl LocusFeaturesSNV {
    fn count(&mut self, read: &Rc<Record>, base: char, refr: char, alt: char, query_pos: u64) {
        if base == refr {
            self.ref_read_features.count(&read, query_pos);
        } else if base == alt{
            self.alt_read_features.count(&read, query_pos);
        }
        self.all_read_features.count(&read, query_pos);
    }

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

#[derive(Debug, Clone)]
struct Variant {
	chrom: String,
	pos: u64,
	refr: String,
	alt: String,
    vartype: VarType,
	counts: BaseCounts,
    read_features: LocusFeaturesSNV,
}

impl Variant {
    fn base_counts_stats(&self) -> Option<BaseCountsStats> {
        let refr_char = self.refr.chars().next()?;
        let alt_char = self.alt.chars().next()?;
        Some(self.counts.stats(refr_char, alt_char))
    }

    fn get_pos_fraction(&self, ref_seq_len: u64) -> f64 {
        self.pos as f64/ ref_seq_len as f64
    }
}

fn get_vcf_min_max(variants: &VecDeque<Variant>) -> Option<(String, u64, u64)> {
    let first = variants.front()?;
    let chrom = first.chrom.clone();

    let min_pos = first.pos;
    let max_pos = variants.back()?.pos;

    Some((chrom, min_pos, max_pos))
}

#[derive(Debug, Clone, Copy)]
enum VarType {
    Snv,
    Del,
    Ins,
    Unknown,
}

impl VarType {
    fn as_str(&self) -> &'static str {
        match self {
            VarType::Snv => "SNV",
            VarType::Del => "DEL",
            VarType::Ins => "INS",
            VarType::Unknown => "UNKNOWN",            
        }
    }
}

fn get_var_type(refr: &str, alt: & str) -> VarType {
    if refr.len() == 1 && alt.len() == 1 {
        VarType::Snv
    } else if refr.len() > alt.len() {
        VarType::Del
    } else if refr.len() < alt.len() {
        VarType::Ins
    } else {
        eprintln!(
            "Warning: Cannot determine the type of variant with ref: {} and alt: {}",
            refr, alt
        );
        VarType::Unknown
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let vcf_path = "/home/jch/genomic-data/chrM_heavy_stress.vcf";
    let bam_path = "/home/jch/genomic-data/chrM_heavy_stress.sorted.bam";
    let varclass = "SNV";
    let csv_path = "/home/jch/git/rust-varlap/test_data/rust-varlap.output.csv";
    
    let mut variants = vcf_reader(&vcf_path, &varclass)?;

    let (region_chrom, min_pos, max_pos) = 
        get_vcf_min_max(&variants).ok_or("Could not determine VCF min/max")?;

    println!("Region Chromosome: {}, Min Pos: {}, Max Pos: {}", &region_chrom, min_pos, max_pos);

    process_bam_region(&mut variants, &bam_path, &region_chrom, min_pos, max_pos, &csv_path)?;

    Ok(())
}
