use std::fs::{File};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use csv::{ByteRecord, ReaderBuilder, Writer, WriterBuilder};
use serde::Serialize;
use tracing::info;

use crate::features::{LocusFeatures, NormalizedLocusFeaturesRow};
use crate::variant::{VarClass, Variant};

pub fn make_output_csv_paths(
    csv_path: &Path,
    reads_paths: &[PathBuf],
    labels: &[Option<String>],
) -> Result<Vec<PathBuf>> {
    let parent = csv_path.parent().unwrap_or_else(|| Path::new(""));

    let stem = csv_path
        .file_stem()
        .and_then(|s| s.to_str())
        .with_context(|| format!("Invalid file stem for '{}'", csv_path.display()))?;

    let extension = csv_path.extension().and_then(|s| s.to_str());

    let mut paths = Vec::with_capacity(reads_paths.len());

    for (i, reads_path) in reads_paths.iter().enumerate() {
        let suffix = match labels.get(i).and_then(|l| l.as_deref()) {
            Some(label) if !label.is_empty() => Some(label.to_string()),
            _ => reads_path
                .file_stem()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string()),
        };

        let path = match suffix {
            Some(suffix) => {
                let filename = match extension {
                    Some(ext) => format!("{stem}_{suffix}.{ext}"),
                    None => format!("{stem}_{suffix}"),
                };
                parent.join(filename)
            }
            None => csv_path.to_path_buf(),
        };

        paths.push(path);
    }

    Ok(paths)
}

pub fn write_header(
    writer: &mut Writer<File>,
    reads_path: &Path,
    label: Option<&str>,
    varclass: &VarClass,
) -> Result<()> {
    let label_prefix = match label {
        Some(label) => label,
        None => reads_path
            .file_name()
            .context("Invalid or missing file name")?
            .to_str()
            .context("File name contains invalid UTF-8")?,
    };

    let dynamic_headers = match varclass {
        VarClass::Snv => HEADER_FIELDS_SNV,
        VarClass::Indel => HEADER_FIELDS_INDEL,
    };

    let labelled_dynamic: Vec<String> = dynamic_headers
        .iter()
        .map(|field| format!("{label_prefix} {field}"))
        .collect();

    writer
        .write_record(
            HEADER_FIELDS_SHARED
                .iter()
                .copied()
                .chain(labelled_dynamic.iter().map(String::as_str)),
        )
        .context("Failed to write CSV header")?;

    Ok(())
}

pub fn write_variant_row(
    writer: &mut Writer<File>,
    var: &Variant,
    pos_normalized: f64,
    sample: Option<&str>,
) -> Result<()> {
    writer
        .serialize(OutputRow::from_variant(var, pos_normalized, sample))
        .with_context(|| {
            format!(
                "Failed to write variant '{}:{}' to output CSV",
                var.info.chrom, var.info.pos
            )
        })?;

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
    ref_normalized_read_position: Option<f64>,

    alt_nm: Option<f64>,
    alt_base_qual: Option<f64>,
    alt_map_qual: Option<f64>,
    alt_align_len: Option<f64>,
    alt_clipping: Option<f64>,
    alt_indel: Option<f64>,
    alt_forward_strand: Option<f64>,
    alt_reverse_strand: Option<f64>,
    alt_supplementary: Option<f64>,
    alt_normalized_read_position: Option<f64>,

    all_nm: Option<f64>,
    all_base_qual: Option<f64>,
    all_map_qual: Option<f64>,
    all_align_len: Option<f64>,
    all_clipping: Option<f64>,
    all_indel: Option<f64>,
    all_forward_strand: Option<f64>,
    all_reverse_strand: Option<f64>,
    all_supplementary: Option<f64>,
    all_normalized_read_position: Option<f64>,
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
            ref_normalized_read_position: rf.ref_normalized_read_position,

            alt_nm: rf.alt_nm,
            alt_base_qual: rf.alt_base_qual,
            alt_map_qual: rf.alt_map_qual,
            alt_align_len: rf.alt_align_len,
            alt_clipping: rf.alt_clipping,
            alt_indel: rf.alt_indel,
            alt_forward_strand: rf.alt_forward_strand,
            alt_reverse_strand: rf.alt_reverse_strand,
            alt_supplementary: rf.alt_supplementary,
            alt_normalized_read_position: rf.alt_normalized_read_position,

            all_nm: rf.all_nm,
            all_base_qual: rf.all_base_qual,
            all_map_qual: rf.all_map_qual,
            all_align_len: rf.all_align_len,
            all_clipping: rf.all_clipping,
            all_indel: rf.all_indel,
            all_forward_strand: rf.all_forward_strand,
            all_reverse_strand: rf.all_reverse_strand,
            all_supplementary: rf.all_supplementary,
            all_normalized_read_position: rf.all_normalized_read_position,
        }
    }
}

#[derive(Serialize)]
enum OutputRow<'a> {
    Snv(OutputRowSNV<'a>),
    Indel(OutputRowINDEL<'a>),
}

impl<'a> OutputRow<'a> {
    fn from_variant(var: &'a Variant, pos_normalized: f64, sample: Option<&'a str>) -> Self {
        match &var.features {
            LocusFeatures::Snv(_) => {
                OutputRow::Snv(OutputRowSNV::from_variant(var, pos_normalized, sample))
            }
            LocusFeatures::Indel(_) => {
                OutputRow::Indel(OutputRowINDEL::from_variant(var, pos_normalized, sample))
            }
        }
    }
}

#[derive(Serialize)]
struct OutputRowSNV<'a> {
    chrom: &'a str,
    pos: u64,
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalized: f64,
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

impl<'a> OutputRowSNV<'a> {
    fn from_variant(var: &'a Variant, pos_normalized: f64, sample: Option<&'a str>) -> Self {
        let bcs = var
            .base_counts_stats()
            .expect("Could not get base count statistics.");
        let read_features = OutputReadFeatures::from(var.normalized_row());

        let f = match &var.features {
            LocusFeatures::Snv(f) => f,
            LocusFeatures::Indel(_) => panic!("Expected SNV features"),
        };

        Self {
            chrom: &var.info.chrom,
            pos: var.info.pos,
            refr: &var.info.refr,
            alt: &var.info.alt,
            vartype: var.info.vartype.as_str(),
            pos_normalized,
            sample,
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
    refr: &'a str,
    alt: &'a str,
    vartype: &'a str,
    pos_normalized: f64,
    sample: Option<&'a str>,
    depth: u32,

    ref_count: u32,
    alt_count: u32,
    other_count: u32,
    alt_vaf: f64,
    overlapping_indels_count: u64,

    read_features: OutputReadFeatures,
}

impl<'a> OutputRowINDEL<'a> {
    fn from_variant(var: &'a Variant, pos_normalized: f64, sample: Option<&'a str>) -> Self {
        let read_features = OutputReadFeatures::from(var.normalized_row());

        let stats = var
            .indel_stats()
            .expect("Could not get indel count statistics.");

        Self {
            chrom: &var.info.chrom,
            pos: var.info.pos,
            refr: &var.info.refr,
            alt: &var.info.alt,
            vartype: var.info.vartype.as_str(),
            pos_normalized,
            sample,
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

const READ_BUF_SIZE: usize = 1024 * 1024; // 1 MiB per input file
const WRITE_BUF_SIZE: usize = 1024 * 1024;

// Horizontally merges `input_paths` into `output_path`. The first file's
// columns are kept in full; every subsequent file has its leading
// `HEADER_FIELDS_SHARED` columns dropped before being appended.
pub fn merge_output_csvs(input_paths: &[PathBuf], output_path: &Path) -> Result<()> {
    if input_paths.is_empty() {
        bail!("Internal error: cannot find any input paths.");
    }
 
    info!(
        "Merging {} output CSVs", input_paths.len());
 
    let n_shared = HEADER_FIELDS_SHARED.len();
 
    let mut readers: Vec<_> = input_paths
        .iter()
        .map(|p| -> Result<_> {
            let file = File::open(p)
                .with_context(|| format!("Failed to open input CSV {}", p.display()))?;
            Ok(ReaderBuilder::new()
                .has_headers(true)
                .buffer_capacity(READ_BUF_SIZE)
                .from_reader(file))
        })
        .collect::<Result<_>>()?;
 
    let out_file = File::create(output_path)
        .with_context(|| format!("Failed to create output CSV {}", output_path.display()))?;
    let mut wtr = WriterBuilder::new()
        .buffer_capacity(WRITE_BUF_SIZE)
        .from_writer(out_file);
 
    // Merge & write the header.
    let mut out_header = ByteRecord::new();
    for (i, rdr) in readers.iter_mut().enumerate() {
        let header = rdr
            .byte_headers()
            .with_context(|| format!("Failed to read header from {}", input_paths[i].display()))?;
        if i == 0 {
            out_header.extend(header.iter());
        } else {
            out_header.extend(header.iter().skip(n_shared));
        }
    }
    wtr.write_byte_record(&out_header)
        .context("Failed to write merged CSV header")?;
 
    // Stream rows in lockstep, reusing buffers to avoid per-row allocation.
    let mut records: Vec<ByteRecord> = vec![ByteRecord::new(); readers.len()];
    let mut out_record = ByteRecord::new();
    let mut row_num: u64 = 0;
 
    loop {
        let has_more = readers[0]
            .read_byte_record(&mut records[0])
            .with_context(|| format!("Failed reading row {row_num} from {}", input_paths[0].display()))?;
        if !has_more {
            break;
        }
        for (rdr, (rec, path)) in readers[1..]
            .iter_mut()
            .zip(records[1..].iter_mut().zip(input_paths[1..].iter()))
        {
            rdr.read_byte_record(rec)
                .with_context(|| format!("Failed reading row {row_num} from {}", path.display()))?;
        }
 
        out_record.clear();
        out_record.extend(records[0].iter());
        for rec in records[1..].iter() {
            out_record.extend(rec.iter().skip(n_shared));
        }
        wtr.write_byte_record(&out_record)
            .with_context(|| format!("Failed writing merged row {row_num}"))?;
 
        row_num += 1;
    }
 
    wtr.flush().context("Failed to flush merged CSV writer")?;
 
    info!(
        "Merging succesful. CSV output can be found at: {:?}",
        output_path.display()
    );
 
    // Clean up redundant per-file outputs only if writer flushe is successful
    for p in input_paths {
        std::fs::remove_file(p)
            .with_context(|| format!("Failed to remove original output CSV: {:?}", p.display()))?;
        info!("Removed original output CSV: {:?}", p.display());
    }
 
    Ok(())
}

// Header fields
const HEADER_FIELDS_SHARED: &[&str] = &[
    "chrom",
    "pos",
    "ref",
    "alt",
    "vartype",
    "pos_normalized",
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
    "ref_normalized_read_position",
    "alt_avg_nm",
    "alt_avg_base_qual",
    "alt_avg_map_qual",
    "alt_avg_align_len",
    "alt_avg_clipping",
    "alt_avg_indel",
    "alt_forward_strand",
    "alt_reverse_strand",
    "alt_supplementary",
    "alt_normalized_read_position",
    "all_avg_nm",
    "all_avg_base_qual",
    "all_avg_map_qual",
    "all_avg_align_len",
    "all_avg_clipping",
    "all_avg_indel",
    "all_forward_strand",
    "all_reverse_strand",
    "all_supplementary",
    "all_normalized_read_position",
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
    "ref_normalized_read_position",
    "alt_avg_nm",
    "alt_avg_base_qual",
    "alt_avg_map_qual",
    "alt_avg_align_len",
    "alt_avg_clipping",
    "alt_avg_indel",
    "alt_forward_strand",
    "alt_reverse_strand",
    "alt_supplementary",
    "alt_normalized_read_position",
    "all_avg_nm",
    "all_avg_base_qual",
    "all_avg_map_qual",
    "all_avg_align_len",
    "all_avg_clipping",
    "all_avg_indel",
    "all_forward_strand",
    "all_reverse_strand",
    "all_supplementary",
    "all_normalized_read_position",
];
