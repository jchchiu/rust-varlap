use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result};
use csv::Writer;
use serde::Serialize;
use std::error::Error;

use crate::features::{LocusFeatures, NormalizedLocusFeaturesRow};
use crate::variant::{Variant, VarClass};

pub fn write_header(
    writer: &mut Writer<File>,
    reads_path: &Path,
    label: Option<&str>,
    varclass: &VarClass,
) -> Result<()> {
    let label_prefix = match label {
        Some(label) => label,
        None => {
            reads_path
                .file_name()
                .context("Invalid or missing file name")?
                .to_str()
                .context("File name contains invalid UTF-8")?
        },
    };

    let dynamic_headers = match varclass {
        VarClass::Snv => HEADER_FIELDS_SNV,
        VarClass::Indel => HEADER_FIELDS_INDEL,
    };

    let labelled_dynamic: Vec<String> = dynamic_headers
        .iter()
        .map(|field| format!("{label_prefix} {field}"))
        .collect();

    writer.write_record(
        HEADER_FIELDS_SHARED
            .iter()
            .copied()
            .chain(labelled_dynamic.iter().map(String::as_str)),
    )?;

    Ok(())
}

pub fn write_variant_row(
    writer: &mut Writer<File>,
    var: &Variant,
    pos_fraction: f64,
    sample: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    writer.serialize(OutputRow::from_variant(var, pos_fraction, sample))?;

    Ok(())
}

#[derive(Serialize)]
struct OutputReadFeatures {
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

impl From<NormalizedLocusFeaturesRow> for OutputReadFeatures {
    fn from(rf: NormalizedLocusFeaturesRow) -> Self {
        Self {
            ref_nm: rf.ref_nm,
            ref_base_qual: rf.ref_base_qual,
            ref_map_qual: rf.ref_map_qual,
            ref_align_len: rf.ref_align_len,
            ref_clipping: rf.ref_clipping,
            ref_indel: rf.ref_indel,
            ref_forward_strand: rf.ref_forward_strand,
            ref_reverse_strand: rf.ref_reverse_strand,
            ref_supplementary: rf.ref_supplementary,
            ref_normalised_read_position: rf.ref_normalised_read_position,

            alt_nm: rf.alt_nm,
            alt_base_qual: rf.alt_base_qual,
            alt_map_qual: rf.alt_map_qual,
            alt_align_len: rf.alt_align_len,
            alt_clipping: rf.alt_clipping,
            alt_indel: rf.alt_indel,
            alt_forward_strand: rf.alt_forward_strand,
            alt_reverse_strand: rf.alt_reverse_strand,
            alt_supplementary: rf.alt_supplementary,
            alt_normalised_read_position: rf.alt_normalised_read_position,

            all_nm: rf.all_nm,
            all_base_qual: rf.all_base_qual,
            all_map_qual: rf.all_map_qual,
            all_align_len: rf.all_align_len,
            all_clipping: rf.all_clipping,
            all_indel: rf.all_indel,
            all_forward_strand: rf.all_forward_strand,
            all_reverse_strand: rf.all_reverse_strand,
            all_supplementary: rf.all_supplementary,
            all_normalised_read_position: rf.all_normalised_read_position,
        }
    }
}

#[derive(Serialize)]
enum OutputRow<'a> {
    Snv(OutputRowSNV<'a>),
    Indel(OutputRowINDEL<'a>),
}

impl<'a> OutputRow<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64, sample: Option<&'a str>) -> Self {
        match &var.features {
            LocusFeatures::Snv(_) => OutputRow::Snv(OutputRowSNV::from_variant(var, pos_fraction, sample)),
            LocusFeatures::Indel(_) => OutputRow::Indel(OutputRowINDEL::from_variant(var, pos_fraction, sample)),
        }
    }
}

#[derive(Serialize)]
struct OutputRowSNV<'a> {
    chrom: &'a str,
    pos: u64,
    #[serde(rename = "ref")]
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalised: f64,
    sample: Option<&'a str>,
    depth: u32,

    count_a: u32,
    count_t: u32,
    count_g: u32,
    count_c: u32,
    count_n: u32,

    ref_count: u32,
    alt_count: u32,
    alt_vaf: f64,
    
    read_features: OutputReadFeatures,
}

// Make branching for snv/indel; make common function for read_features
impl<'a> OutputRowSNV<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64, sample: Option<&'a str>) -> Self {
        let bcs = var.base_counts_stats().expect("Could not get base count statistics.");
        let read_features = OutputReadFeatures::from(var.normalized_row());

        let f = match &var.features {
            LocusFeatures::Snv(f) => f,
            LocusFeatures::Indel(_) => panic!("Expected SNV features"),
        };

        Self {
            chrom: &var.chrom,
            pos: var.pos,
            refr: &var.refr,
            alt: &var.alt,
            vartype: &var.vartype.as_str(),
            pos_normalised: pos_fraction,
            sample: sample,
            depth: bcs.depth,
            
            count_a: f.base_counts.a,
            count_t: f.base_counts.t,
            count_g: f.base_counts.g,
            count_c: f.base_counts.c,
            count_n: f.base_counts.n,

            ref_count: bcs.ref_count,
            alt_count: bcs.alt_count,
            alt_vaf: bcs.alt_vaf,

            read_features,
        }
    }
}

#[derive(Serialize)]
struct OutputRowINDEL<'a> {
    chrom: &'a str,
    pos: u64,
    #[serde(rename = "ref")]
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalised: f64,
    sample: Option<&'a str>,
    depth: u32,

    ref_count: u32,
    alt_count: u32,
    other_count: u32,
    alt_vaf: f64,
    overlapping_indels_count: u64,

    read_features: OutputReadFeatures,
}

// Make branching for snv/indel; make common function for read_features
impl<'a> OutputRowINDEL<'a> {
    fn from_variant(var: &'a Variant, pos_fraction: f64, sample: Option<&'a str>) -> Self {
        let read_features = OutputReadFeatures::from(var.normalized_row());

        let stats = var.indel_stats().expect("Could not get indel count statistics.");

        Self {
            chrom: &var.chrom,
            pos: var.pos,
            refr: &var.refr,
            alt: &var.alt,
            vartype: &var.vartype.as_str(),
            pos_normalised: pos_fraction,
            sample: sample,
            depth: stats.depth,

            ref_count: stats.ref_count,
            alt_count: stats.alt_count,
            other_count: stats.other_count,
            alt_vaf: stats.alt_vaf,
            overlapping_indels_count: stats.overlapping_indels_count,

            read_features,
        }
    }
}

// Header fields
const HEADER_FIELDS_SHARED: &[&str] = &[
    "chrom",
    "pos",
    "ref",
    "alt",
    "vartype",
    "pos_normalised",
    "sample",
];

const HEADER_FIELDS_SNV: &[&str] = &[
    "depth",
    "count_a",
    "count_t",
    "count_g",
    "count_c",
    "count_n",
    "ref_count",
    "alt_count",
    "alt_vaf",
    "ref_avg_nm",
    "ref_avg_base_qual",
    "ref_avg_map_qual",
    "ref_avg_align_len",
    "ref_avg_clipping",
    "ref_avg_indel",
    "ref_forward_strand",
    "ref_reverse_strand",
    "ref_supplementary",
    "ref_normalised_read_position",
    "alt_avg_nm",
    "alt_avg_base_qual",
    "alt_avg_map_qual",
    "alt_avg_align_len",
    "alt_avg_clipping",
    "alt_avg_indel",
    "alt_forward_strand",
    "alt_reverse_strand",
    "alt_supplementary",
    "alt_normalised_read_position",
    "all_avg_nm",
    "all_avg_base_qual",
    "all_avg_map_qual",
    "all_avg_align_len",
    "all_avg_clipping",
    "all_avg_indel",
    "all_forward_strand",
    "all_reverse_strand",
    "all_supplementary",
    "all_normalised_read_position",
];

const HEADER_FIELDS_INDEL: &[&str] = &[
    "depth",
    "ref_count",
    "alt_count",
    "alt_vaf",
    "other_count",
    "overlapping_indels_count",
    "ref_avg_nm",
    "ref_avg_base_qual",
    "ref_avg_map_qual",
    "ref_avg_align_len",
    "ref_avg_clipping",
    "ref_avg_indel",
    "ref_forward_strand",
    "ref_reverse_strand",
    "ref_supplementary",
    "ref_normalised_read_position",
    "alt_avg_nm",
    "alt_avg_base_qual",
    "alt_avg_map_qual",
    "alt_avg_align_len",
    "alt_avg_clipping",
    "alt_avg_indel",
    "alt_forward_strand",
    "alt_reverse_strand",
    "alt_supplementary",
    "alt_normalised_read_position",
    "all_avg_nm",
    "all_avg_base_qual",
    "all_avg_map_qual",
    "all_avg_align_len",
    "all_avg_clipping",
    "all_avg_indel",
    "all_forward_strand",
    "all_reverse_strand",
    "all_supplementary",
    "all_normalised_read_position",
];