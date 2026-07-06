use std::collections::VecDeque;

use anyhow::{Result};
// use tracing::{debug, info, warn};

// use crate::errors::AppError;
use crate::features::{LocusFeatures, LocusFeaturesIndel, LocusFeaturesSnv};
use crate::variant::{BinnedVariants, ParsedVariants, VarType, Variant, VariantBin};

pub fn bin<'a>(parsed_variants: &'a ParsedVariants) -> Result<BinnedVariants<'a>> {
    let mut bins: Vec<VariantBin> = Vec::new();
            // Get unique chromosomes and their counts for variants addded to queue
            // NOTE: Variants file MUST be sorted in ascending order
            // We do not need to get the index as we are popping the queue when iterating over variants
    
    for variant in &parsed_variants.variants {
        if let Some(last) = bins.last_mut()
            && last.chrom == variant.chrom {
                last.variants.push_back(match variant.vartype {
                    VarType::Snv => Variant {
                        info: variant,
                        features: LocusFeatures::Snv(LocusFeaturesSnv::default()),
                    },
                    VarType::Del | VarType::Ins => Variant {
                        info: variant,
                        features: LocusFeatures::Indel(LocusFeaturesIndel::default()),
                    },
                    VarType::Unknown => unreachable!(),
                });
                
                continue;
            }

        bins.push(VariantBin {
            chrom: variant.chrom.clone(),
            variants: VecDeque::from([match variant.vartype {
                VarType::Snv => Variant {
                    info: variant,
                    features: LocusFeatures::Snv(LocusFeaturesSnv::default()),
                },
                VarType::Del | VarType::Ins => Variant {
                    info: variant,
                    features: LocusFeatures::Indel(LocusFeaturesIndel::default()),
                },
                VarType::Unknown => unreachable!(),
            }]),
        });
    }

    Ok(BinnedVariants { bins })
}
