use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::VecDeque;
use std::error::Error;

use crate::variant::{VarClass, VarType};
use crate::{Variant};
use crate::features::{LocusFeatures, LocusFeaturesIndel, LocusFeaturesSnv};

pub fn parse(file_path: &str, varclass: &VarClass) -> Result<VecDeque<Variant>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

	let mut variants = VecDeque::new();

    for line_result in reader.lines() {
        let line = line_result?;

        if line.starts_with("#") {
            continue;
        }

        // if line.starts_with("#") && !is_valid_vcf_header_line(&line) {
        //     match {
        //         Ok() => continue,
        //         Err(error) => panic!("Invalid VCF: {error:?}"),
        //     };
        // }
        
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

                match varclass {
                    VarClass::Snv => {
                        if matches!(vartype, VarType::Snv) {
                            variants.push_back(Variant {
                                chrom: chrom.clone(),
                                pos,
                                refr: refr.clone(),
                                alt: alt.to_string(),
                                vartype,
                                features: LocusFeatures::Snv(LocusFeaturesSnv::default()),
                            });
                        }
                    }
                    VarClass::Indel => {
                        if matches!(vartype, VarType::Ins | VarType::Del) {
                            variants.push_back(Variant {
                                chrom: chrom.clone(),
                                pos,
                                refr: refr.clone(),
                                alt: alt.to_string(),
                                vartype,
                                features: LocusFeatures::Indel(LocusFeaturesIndel::default()),
                            });
                        }
                    }
                }
            }
        } else {
            eprintln!("Warning: Skipping input row: {}", line);
        }
    }
			
	Ok(variants)
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

// fn is_valid_vcf_header_line(line: &str) -> bool {
//     let expected = ["#CHROM", "POS", "ID", "REF", "ALT"];
//     line.split_whitespace().take(5).eq(expected.into_iter())
// }