use std::fs::File;
use std::path::PathBuf;
use std::io::{BufRead, BufReader, Read, Seek};
use std::collections::VecDeque;
use std::error::Error;
use flate2::read::MultiGzDecoder;

use crate::variant::{VarClass, VarType};
use crate::{Variant};
use crate::features::{LocusFeatures, LocusFeaturesIndel, LocusFeaturesSnv};

pub fn parse(file_path: &PathBuf, varclass: &VarClass) -> Result<VecDeque<Variant>, Box<dyn Error>> {
    let reader = check_gzip(file_path)?;

	let mut variants = VecDeque::new();

    for line_result in reader.lines() {
        let line = line_result?;

        if line.starts_with("#") {
            continue;
        }

        // if line.starts_with("#CHROM") && !is_valid_vcf_header_line(&line) {
        //     break;
        // } else {
        //     continue
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

                if is_acceptable_variant(&varclass, &vartype, &refr, &alt){
                    match varclass {
                        VarClass::Snv => {
                            variants.push_back(Variant {
                                chrom: chrom.clone(),
                                pos,
                                refr: refr.clone(),
                                alt: alt.to_string(),
                                vartype,
                                features: LocusFeatures::Snv(LocusFeaturesSnv::default()),
                            });
                        }
                        VarClass::Indel => {
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

// USING GZIP HANDLER AS PER RUSTQC
// https://github.com/seqeralabs/RustQC/blob/main/src/io.rs
// REWRITE MYSELF LATER

// https://cseweb.ucsd.edu/classes/sp22/cse223B-a/tribbler/flate2/read/struct.MultiGzDecoder.html
// Need to use multigzdecoder

/// Gzip magic bytes: the first two bytes of any gzip-compressed file.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

fn check_gzip(file_path: &PathBuf) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    let mut file = File::open(file_path).unwrap();

    // Read the first two bytes to check for gzip magic number
    let mut magic = [0u8; 2];
    let bytes_read = file
        .read(&mut magic)
        .unwrap();

    // Seek back to the beginning so the reader starts from byte 0
    file.seek(std::io::SeekFrom::Start(0)).unwrap();

    if bytes_read >= 2 && magic == GZIP_MAGIC {
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

// fn is_valid_vcf_header_line(line: &str) -> bool {
//     let expected = ["#CHROM", "POS", "ID", "REF", "ALT"];
//     line.split_whitespace().take(5).eq(expected)
// }

fn is_acceptable_variant(
    varclass: &VarClass,
    vartype: &VarType,
    refr: &str,
    alt: &str,
    // max_indel_size: u32,
) -> bool {
    if !is_only_dna_bases(refr) || !is_only_dna_bases(alt) {
        false
    } else if !is_desired_type(varclass, vartype) {
        false
    // } else if !is_within_max_size(varclass, max_indel_size, refr, alt) {
    //     false
    } else if matches!(varclass, VarClass::Indel) && !is_valid_indel(refr, alt) {
        false
    } else {
        true
    }
}

fn is_only_dna_bases(sequence: &str) -> bool {
    sequence
        .chars()
        .all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'T' | 'G' | 'C'))
}

fn is_desired_type(varclass: &VarClass, vartype: &VarType) -> bool {
    match varclass {
        VarClass::Snv => matches!(vartype, VarType::Snv),
        VarClass::Indel => matches!(vartype, VarType::Ins | VarType::Del),
    }
}

// fn is_within_max_size(
//     varclass,
//     max_indel_size,
//     refr,
//     alt
// ) -> bool {

// }

fn is_valid_indel(refr: &str, alt: &str) -> bool {
    match refr.len().cmp(&alt.len()) {
        std::cmp::Ordering::Equal => false,
        std::cmp::Ordering::Less => {
            !refr.is_empty() && alt.starts_with(refr)
        }
        std::cmp::Ordering::Greater => {
            !alt.is_empty() && refr.starts_with(alt)
        }
    }
}