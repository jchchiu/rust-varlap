use std::collections::VecDeque;
// use std::path::Path;

// use anyhow::{Context, Result};
use anyhow::Result;
// use rust_htslib::bam::Reader;
// use rust_htslib::bam::{Record, Read};
use tracing::{info};

use crate::errors::AppError;
use crate::features::{LocusFeatures, LocusFeaturesIndel, LocusFeaturesSnv};
// use crate::read_parser::{detect_file_type, FileType};
use crate::variant::{BinnedVariants, ParsedVariants, VarType, Variant, VariantInfo, VariantBin};

pub fn bin<'a>(
    parsed_variants: &'a ParsedVariants,
    // reads_path: &Path,
    // fasta_path: Option<&Path>,
    gap: Option<u64>,
) -> Result<BinnedVariants<'a>> {
    let mut bins: Vec<VariantBin> = Vec::new();

    // let file_type = detect_file_type(reads_path)
    //     .with_context(|| format!("Failed to detect reads file type for '{}'", reads_path.display()))?;

    // let mut reader = Reader::from_path(reads_path)
    //     .with_context(|| format!("Failed to open reads file '{}'", reads_path.display()))?;

    // match file_type {
    //     FileType::Bam => {}
    //     FileType::Cram => {
    //         let fasta = fasta_path.ok_or(AppError::MissingCramReference)?;

    //         reader
    //             .set_reference(fasta)
    //             .with_context(|| format!("Failed to set CRAM reference to '{}'", fasta.display()))?;
    //     },
    // };

    // let mean_read_len = estimate_mean_read_len(&mut reader, 1000)?;
    // let max_hyperparameter = 5.0;
    // let max_gap = (max_hyperparameter * mean_read_len).round() as u64;

    // println!("Calculated Mean Read Len: {:?}
    //           Set Max Hyperparameter:   {:?}
    //           Calculated Max Gap:       {:?}"
    //           , mean_read_len, max_hyperparameter, max_gap);

    // let default_gap = 100000u64;

    // Use the gap provided, or set the bin size gap to default of 100kb
    let gap = gap.unwrap_or(DEFAULT_GAP);

    info!("Set gap between bins as {:?}bp", gap);
    
    for variant in &parsed_variants.variants {
        let should_append = if let Some(last_bin) = bins.last() {
            if last_bin.chrom != variant.chrom {
                false
            } else if let Some(last_variant) = last_bin.variants.back() {
                let distance = get_variant_distance(variant, &last_variant)?;

                distance <= gap
            } else {
                true
            }
        } else {
            false
        };

        let variant_features = make_variant_features(variant);

        if should_append {
            bins.last_mut()
                .expect("Internal error: cannot find any initialized bins.
                         INVARIANT BROKEN: 'if' boolean check already assumes that a previous bin has been found")
                .variants
                .push_back(variant_features);
        } else {
            bins.push(VariantBin {
                chrom: variant.chrom.clone(),
                variants: VecDeque::from([variant_features]),
            });
        }
    }

    info!("Created {:?} bins for variants successfully", bins.len());

    Ok(BinnedVariants { bins })
}

const DEFAULT_GAP: u64 = 100_000;

fn get_variant_distance(
    current: &VariantInfo,
    previous: &Variant,
) -> Result<u64, AppError> {
    current.pos
        .checked_sub(previous.info.pos)
        .ok_or_else(|| AppError::UnsortedVariants {
            chromosome: current.chrom.clone(),
            error_pos: current.pos,
            previous_pos: previous.info.pos,
        })
}

// fn estimate_mean_read_len(
//     reader: &mut Reader,
//     n: usize,
// ) -> Result<f64> {
//     let mut record = Record::new();
//     let mut total = 0usize;
//     let mut count = 0usize;

//     while count < n {
//         match reader.read(&mut record) {
//             Some(Ok(())) => {
//                 total += record.seq_len();
//                 count += 1;
//             }
//             Some(Err(e)) => return Err(e.into()),
//             None => break,
//         }
//     }

//     Ok(total as f64 / count as f64)
// }

fn make_variant_features<'a>(variant: &'a VariantInfo) -> Variant<'a> {
    match variant.vartype {
        VarType::Snv => Variant {
            info: variant,
            features: LocusFeatures::Snv(LocusFeaturesSnv::default()),
        },
        VarType::Del | VarType::Ins => Variant {
            info: variant,
            features: LocusFeatures::Indel(LocusFeaturesIndel::default()),
        },
        VarType::Unknown => unreachable!(),
    }
}