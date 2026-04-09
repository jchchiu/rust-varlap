use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::collections::VecDeque;
use std::error::Error;

fn vcf_reader(file_path: &str) -> Result<VecDeque<Variant>, Box<dyn Error>> {
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
                let vartype = get_var_type(&refr, &alt);

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

#[derive(Debug, Clone, Default)]
struct BaseCounts {
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
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
    let file_path = "test_data/vars.vcf";
    
    let variants = vcf_reader(file_path)?;

    let (region_chrom, min_pos, max_pos) = 
        get_vcf_min_max(&variants).ok_or("Could not determine VCF min/max")?;

    println!("Region Chromosome: {}, Min Pos: {}, Max Pos: {}", region_chrom, min_pos, max_pos);
    println!("{:?}", variants);

    Ok(())
}
