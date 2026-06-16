mod cli;
mod variant_parser;
mod variant;
mod features;
mod bam_parser;
mod output;

use crate::variant::Variant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let vcf_path = "test_data/chrM_heavy_stress.vcf";
    // let bam_path: &str = "test_data/chrM_heavy_stress.sorted.bam";
    // let varclass = "SNV";
    // let csv_path = "test_data/rust-varlap.output.csv";
    // let sample = "";

    let args = cli::parse();
    
    let mut variants = variant_parser::parse(&args.variant_file, &args.varclass)?;

    // println!("Region Chromosome: {}, Min Pos: {}, Max Pos: {}", &region_chrom, min_pos, max_pos);

    // VARCLASS INPUT TEMP FIX FOR CSV HEADER
    bam_parser::parse_region(
        &mut variants, 
        &args.bam_file, 
        &args.csv_path, 
        args.sample.as_deref(),
        args.label.as_deref(), 
        &args.varclass, 
        args.fasta_file.as_deref()
    )?;

    Ok(())
}