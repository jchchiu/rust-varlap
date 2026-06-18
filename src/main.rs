mod cli;
mod variant_parser;
mod variant;
mod features;
mod bam_parser;
mod output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = cli::parse();
    
    let mut parsed_variants = variant_parser::parse(&args.variant_file, &args.varclass)?;

    bam_parser::parse_region(
        &mut parsed_variants, 
        &args.bam_file, 
        &args.csv_path, 
        args.sample.as_deref(),
        args.label.as_deref(), 
        &args.varclass, 
        args.fasta_file.as_deref()
    )?;

    Ok(())
}