mod cli;
mod variant_parser;
mod variant;
mod features;
mod bam_parser;
mod output;
mod errors;

use anyhow::Result;

use crate::errors::{AppError, print_error};

fn run() -> Result<()> {
    let args = cli::parse();

    let mut parsed_variants = 
        variant_parser::parse(&args.variant_file, &args.varclass)?;

    bam_parser::parse_region(
        &mut parsed_variants, 
        &args.bam_file, 
        &args.output_path, 
        args.sample.as_deref(),
        args.label.as_deref(), 
        &args.varclass, 
        args.fasta_file.as_deref()
    )?;

    Ok(())
}

fn main() {
    if let Err(err) = run() {
        let program = env!("CARGO_PKG_NAME");

        if let Some(app_err) = err.downcast_ref::<AppError>() {
            print_error(program, app_err);
            std::process::exit(app_err.exit_code());
        }

        // Everything else is an unexpected/runtime failure.
        eprintln!("{program} ERROR: {:#}", err);
        std::process::exit(1);
    }
}
