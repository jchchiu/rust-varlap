mod cli;
mod vcf_parser;
mod variant;
mod features;
mod bam_parser;
mod output;

use crate::variant::Variant;
use crate::bam_parser::process_bam_region;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let vcf_path = "test_data/chrM_heavy_stress.vcf";
    // let bam_path: &str = "test_data/chrM_heavy_stress.sorted.bam";
    // let varclass = "SNV";
    // let csv_path = "test_data/rust-varlap.output.csv";
    // let sample = "";

    let args = cli::parse();
    
    let mut variants = vcf_parser::parse(&args.vcf, &args.varclass)?;

    // println!("Region Chromosome: {}, Min Pos: {}, Max Pos: {}", &region_chrom, min_pos, max_pos);

    // VARCLASS INPUT TEMP FIX FOR CSV HEADER
    process_bam_region(&mut variants, &args.bams, &args.csv_path, args.sample.as_deref(), &args.varclass)?;

    Ok(())
}