use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::VecDeque;
use std::error::Error;
use csv::Writer;
use rust_htslib::bam::{Read, IndexedReader};

fn write_header_row(writer: &mut Writer<File>) -> Result<(), Box<dyn Error>> {
    writer.write_record(&[
        "chrom",
        "pos",
        "ref",
        "alt",
        "vartype",
        "pos_normalised",
        "depth",
        "A",
        "T",
        "G",
        "C",
        "N",
        "ref_count",
        "alt_count",
        "alt_vaf",
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

fn print_header_row() -> Result<(), Box<dyn Error>> {
    println!("chrom\tpos\tref\talt\tvartype\tpos_normalised\tdepth\tA\tT\tG\tC\tN\tref_count\talt_count\talt_vaf");
    Ok(())
}

fn print_variant_row(var: &Variant, base_stats: &BaseCountsStats, pos_fraction: f64) -> Result<(), Box<dyn Error>> {
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        var.chrom,
        var.pos,
        var.refr,
        var.alt,
        var.vartype.as_str(),
        pos_fraction,
        base_stats.depth,
        var.counts.a,
        var.counts.t,
        var.counts.g,
        var.counts.c,
        var.counts.n,
        base_stats.refr_count,
        base_stats.alt_count,
        base_stats.alt_vaf,
    );
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
    //let header = bam_reader.header().to_owned();

    bam_reader.fetch((region_chrom, min_pos - 1, max_pos))?;

    let mut csv_writer = Writer::from_path(csv_path)?;

    write_header_row(&mut csv_writer)?;

    print_header_row()?;

    let ref_seq_len = get_ref_len(&bam_reader, &region_chrom)?;

    for read_result in bam_reader.rc_records() {
        let record = read_result?;

        //let tid = record.tid();
        //let read_chrom = String::from_utf8(header.tid2name(tid as u32).to_vec())?;

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
                    print_variant_row(&var, &base_counts_stats, pos_fraction)?;
                    write_variant_row(&mut csv_writer, &var, &base_counts_stats, pos_fraction)?;
                }
            } else {
                break;
            }
        }

        for var in &mut *variants {
            //if variant.chrom != read_chrom {
            //    break;
            //}
            
            let zero_based_pos = var.pos - 1;
            let read_end = read_start + record.seq_len() as u64;

            if zero_based_pos >= read_start && zero_based_pos < read_end {
                let base = seq[(zero_based_pos - read_start) as usize] as char;
                var.counts.increment(base);
            } else {
                break;
            }
        }
    }

    while let Some(var) = variants.pop_front() {
        let base_counts_stats = var.base_counts_stats().ok_or("error")?;
        let pos_fraction = var.get_pos_fraction(ref_seq_len);
        print_variant_row(&var, &base_counts_stats, pos_fraction)?;
        write_variant_row(&mut csv_writer, &var, &base_counts_stats, pos_fraction)?;
    }

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
                });
            }
        } else {
            eprintln!("Warning: Skipping input row: {}", line);
        }
    }
			
	Ok(variants)
}

#[derive(Debug, Clone, Copy)]
struct BaseCountsStats {
    depth: u32,
    refr_count: u32,
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
    fn increment(&mut self, base: char) {
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
        let refr_count = self.count_for_base(refr);
        let alt_count = self.count_for_base(alt);
        let alt_vaf = if depth > 0 {
            alt_count as f64 / depth as f64
        } else {
            0.0
        };

        BaseCountsStats {
            depth,
            refr_count,
            alt_count,
            alt_vaf,
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
