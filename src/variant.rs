use std::collections::VecDeque;

use clap::ValueEnum;
use rust_htslib::bam::{Record};

use crate::features::{LocusFeatures, NormalizedLocusFeaturesRow, AlleleCountsSnvStats, AlleleCountsIndelStats};

// NOTE: 
// If we want to process multiple BAMs, we should move the features out
//  of the variant struct (maybe VariantInfo and VariantFeat structs)
//  where VariantFeat has an immutable lifetime borrow of VariantInfo
#[derive(Debug, Clone)]
pub struct Variant {
	pub chrom: String,
	pub pos: u64,
	pub refr: String,
	pub alt: String,
    pub vartype: VarType,
    pub features: LocusFeatures,
}

impl Variant {
    pub fn base_counts_stats(&self) -> Option<AlleleCountsSnvStats> {
        let ref_char = self.refr.chars().next()?;
        let alt_char = self.alt.chars().next()?;

        match &self.features {
            LocusFeatures::Snv(f) => Some(f.base_counts.stats(ref_char, alt_char)),
            LocusFeatures::Indel(_) => None,
        }
    }

    pub fn indel_stats(&self) -> Option<AlleleCountsIndelStats> {
        match &self.features {
            LocusFeatures::Snv(_) => None,
            LocusFeatures::Indel(f) => Some(f.stats()),
        }
    }

    pub fn count_locus_features(&mut self, read: &Record, ref_pos: u64) {
        let qpos = self.ref_pos_to_query_pos(read, ref_pos);

        match &mut self.features {
            LocusFeatures::Snv(f) => {
                if let (Some(refr_char), Some(alt_char)) =
                    (self.refr.chars().next(), self.alt.chars().next())
                {
                    f.count(read, refr_char, alt_char, qpos);
                }
            }
            LocusFeatures::Indel(f) => {
                f.count(read, &self.refr, &self.alt, ref_pos, qpos, &self.vartype);
            }
        }
    }

    pub fn normalized_row(&self) -> NormalizedLocusFeaturesRow {
        self.features.normalized_row()
    }

    pub fn get_pos_fraction(&self, ref_seq_len: u64) -> f64 {
        self.pos as f64/ ref_seq_len as f64
    }

    pub fn ref_pos_to_query_pos (&self, read: &Record, target_pos: u64) -> Option<u32> {
        let cigar = read.cigar();
        Some(cigar.read_pos(target_pos as u32, false, false).ok()?)?
    }
}

// User supplied 
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum VarClass {
    Snv,
    Indel,
}

// Determined using base and reference of variant
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VarType {
    Snv,
    Del,
    Ins,
    Unknown,
}

impl VarType {
    pub fn as_str(&self) -> &'static str {
        match self {
            VarType::Snv => "SNV",
            VarType::Del => "DEL",
            VarType::Ins => "INS",
            VarType::Unknown => "UNKNOWN",            
        }
    }
}

#[derive(Debug, Clone)]
pub struct ChromBucket {
    pub chrom: String,
    pub variants: VecDeque<Variant>,
}

#[derive(Debug, Clone)]
pub struct ParsedVariants  {
    pub chroms: Vec<ChromBucket>,
}
