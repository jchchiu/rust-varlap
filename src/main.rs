use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::collections::VecDeque;
use std::error::Error;

fn vcf_reader(file_path: &str) -> Result<VecDeque<VariantReader>, Box<dyn Error>> {
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
            let pos = fields[1].parse::<u64>();
            let refr = fields[3].to_string();
            let alt = fields[4].to_string();
            
            variants.push_back(VariantReader {
                chrom: chrom,
                pos: pos?,
                refr: refr,
                alt: alt,
                counts: BaseCounts::default(),
            });
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
            _ => self.n += 1,
        }
    }
}

#[derive(Debug, Clone)]
struct VariantReader {
	chrom: String,
	pos: u64,
	refr: String,
	alt: String,
	counts: BaseCounts,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_path = "test_data/small_test.vcf";
    
    let variants = vcf_reader(file_path)?;

    println!("{:?}", variants);

    Ok(())
}
