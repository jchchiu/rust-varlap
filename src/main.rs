mod bam_parser;
mod cli;
mod errors;
mod features;
mod output;
mod variant;
mod variant_parser;

use anyhow::Result;

use crate::errors::{AppError, print_error};

fn run() -> Result<()> {
    let args = cli::parse();

    let mut parsed_variants = 
        variant_parser::parse(&args.variants, &args.varclass)?;

    bam_parser::parse_region(
        &mut parsed_variants, 
        &args.reads, 
        &args.output, 
        args.sample.as_deref(),
        args.label.as_deref(), 
        &args.varclass, 
        args.fasta.as_deref()
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

        // Set everything else as a I/O error for now. FIX LATER
        eprintln!("{program} ERROR: {:#}", err);
        std::process::exit(1);
    }
}
