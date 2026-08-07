use std::collections::VecDeque;
// use std::path::Path;

// use anyhow::{Context, Result};
use anyhow::Result;
// use rust_htslib::bam::Reader;
// use rust_htslib::bam::{Record, Read};
use tracing::{debug, info};

use crate::errors::AppError;
use crate::features::{LocusFeatures, LocusFeaturesIndel, LocusFeaturesSnv};
// use crate::read_parser::{detect_file_type, FileType};
use crate::variant::{BinnedVariants, ParsedVariants, VarType, Variant, VariantBin, VariantInfo};

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
    let gap = match gap {
        Some(gap) => {
            debug!("Using user-provided bin gap: {} bp", gap);
            gap
        }
        None => {
            // Can maybe change it so that a value of -1 has no gap? (just parses whole chromosome)
            debug!("Using default bin gap: {} bp", DEFAULT_GAP);
            DEFAULT_GAP
        }
    };

    for variant in &parsed_variants.variants {
        let should_append = if let Some(last_bin) = bins.last() {
            if last_bin.chrom != variant.chrom {
                false
            } else if let Some(last_variant) = last_bin.variants.back() {
                let distance = get_variant_distance(variant, last_variant)?;

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

    info!(
        "Number of bins created for variants successfully: {}",
        bins.len()
    );

    Ok(BinnedVariants { bins })
}

const DEFAULT_GAP: u64 = 100_000;

fn get_variant_distance(current: &VariantInfo, previous: &Variant) -> Result<u64, AppError> {
    current
        .pos
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

#[cfg(test)]
mod tests {
    use super::*;

    fn variant(chrom: &str, pos: u64, vartype: VarType) -> VariantInfo {
        VariantInfo {
            chrom: chrom.to_string(),
            pos,
            vartype,
            refr: "A".into(),
            alt: "T".into(),
        }
    }

    fn parsed(variants: Vec<VariantInfo>) -> ParsedVariants {
        ParsedVariants { variants }
    }

    // Variant feeature tests
    #[test]
    fn make_variant_features_creates_snv() {
        let v = variant("chr1", 100, VarType::Snv);

        let result = make_variant_features(&v);

        assert!(matches!(result.features, LocusFeatures::Snv(_)));
        assert_eq!(result.info.pos, 100);
    }

    #[test]
    fn make_variant_features_creates_indel_for_insertions() {
        let v = variant("chr1", 100, VarType::Ins);

        let result = make_variant_features(&v);

        assert!(matches!(result.features, LocusFeatures::Indel(_)));
    }

    #[test]
    fn make_variant_features_creates_indel_for_deletions() {
        let v = variant("chr1", 100, VarType::Del);

        let result = make_variant_features(&v);

        assert!(matches!(result.features, LocusFeatures::Indel(_)));
    }

    // Variant distance tests
    #[test]
    fn variant_distance_cases() {
        let cases = [
            // (previous_pos, current_pos, expected_distance)
            ("returns difference", 100, 150, 50),
            ("zero when same position", 100, 100, 0),
        ];

        for (label, previous_pos, current_pos, expected) in cases {
            let previous_variant = variant("chr1", previous_pos, VarType::Snv);
            let previous = make_variant_features(&previous_variant);
            let current = variant("chr1", current_pos, VarType::Snv);

            let distance = get_variant_distance(&current, &previous).unwrap();
            assert_eq!(
                distance, expected,
                "case: {label}, previous_pos={previous_pos}, current_pos={current_pos}"
            );
        }
    }

    #[test]
    fn variant_distance_errors_when_unsorted() {
        let previous_variant = variant("chr1", 200, VarType::Snv);
        let previous = make_variant_features(&previous_variant);
        let current = variant("chr1", 100, VarType::Snv);

        let err = get_variant_distance(&current, &previous).unwrap_err();

        assert!(matches!(
            err,
            AppError::UnsortedVariants {
                chromosome,
                error_pos,
                previous_pos,
            }
            if chromosome == "chr1"
                && error_pos == 100
                && previous_pos == 200
        ));
    }

    // Bin tests
    #[test]
    fn bin_cases() {
        // (label, variants as (chrom, pos, VarType), gap, expected bins as (chrom, variant_count))
        let cases: Vec<(&str, Vec<(&str, u64, VarType)>, Option<u64>, Vec<(&str, usize)>)> = vec![
            ("empty", vec![], Some(10), vec![]),
            ("single variant", vec![("chr1", 100, VarType::Snv)], Some(10), vec![("chr1", 1)]),
            (
                "within gap merges into one bin",
                vec![
                    ("chr1", 100, VarType::Snv),
                    ("chr1", 105, VarType::Snv),
                    ("chr1", 110, VarType::Snv),
                ],
                Some(10),
                vec![("chr1", 3)],
            ),
            (
                "outside gap splits into two bins",
                vec![("chr1", 100, VarType::Snv), ("chr1", 120, VarType::Snv)],
                Some(10),
                vec![("chr1", 1), ("chr1", 1)],
            ),
            (
                "chromosome change always starts a new bin",
                vec![("chr1", 100, VarType::Snv), ("chr2", 101, VarType::Snv)],
                Some(1000),
                vec![("chr1", 1), ("chr2", 1)],
            ),
            (
                "gap exactly equal to threshold stays in same bin",
                vec![("chr1", 100, VarType::Snv), ("chr1", 110, VarType::Snv)],
                Some(10),
                vec![("chr1", 2)],
            ),
            (
                "no gap passed uses default of 100_000",
                vec![("chr1", 100, VarType::Snv), ("chr1", 1_000_000, VarType::Snv)],
                None,
                vec![("chr1", 1), ("chr1", 1)],
            ),
        ];

        for (label, variants, gap, expected_bins) in cases {
            let variants: Vec<_> = variants
                .into_iter()
                .map(|(chrom, pos, vt)| variant(chrom, pos, vt))
                .collect();
            let parsed = parsed(variants);

            let bins = bin(&parsed, gap).unwrap();

            assert_eq!(bins.bins.len(), expected_bins.len(), "case: {label}");
            for (i, (expected_chrom, expected_count)) in expected_bins.iter().enumerate() {
                assert_eq!(&bins.bins[i].chrom, expected_chrom, "case: {label}, bin {i}");
                assert_eq!(
                    bins.bins[i].variants.len(),
                    *expected_count,
                    "case: {label}, bin {i}"
                );
            }
        }
    }
}
