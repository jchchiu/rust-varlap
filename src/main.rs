mod binning;
mod cli;
mod errors;
mod features;
mod output;
mod read_parser;
mod variant;
mod variant_parser;

use anyhow::Result;
use rayon::prelude::*;
use tracing::{info, info_span};

use crate::errors::AppError;

fn run() -> Result<()> {
    let args = cli::parse();

    info!("-------------------------------------------------------------------------------");
    let parsed_variants = variant_parser::parse(&args.variants, &args.varclass)?;

    let output_paths = output::make_output_csv_paths(
        &args.output,
        &args.reads,
        &args.label,
    )?;

    // Initialize the thread pool for analyzing multiple BAM files in parallel
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;

    info!("-------------------------------------------------------------------------------");
    args.reads
        .par_iter()
        .enumerate()
        .try_for_each(|(i, reads)| -> Result<()> {
            let span = info_span!("sample", index = i, path = %reads.display());
            span.in_scope(|| -> Result<()> {

                let label = args.label.get(i).and_then(|l| l.as_deref());

                let mut binned_variants = binning::bin(&parsed_variants, args.gap)?;

                read_parser::parse(
                    &mut binned_variants,
                    reads,
                    &output_paths[i],
                    args.sample.as_deref(),
                    label,
                    &args.varclass,
                    args.fasta.as_deref(),
                )?;

                Ok(())
            })
        })?;

    if args.merge {
        info!("-------------------------------------------------------------------------------");
        output::merge_output_csvs(&output_paths, &args.output)?;
    }

    Ok(())
}

fn main() {
    tracing_subscriber::fmt().init();

    let program = env!("CARGO_PKG_NAME");

    info!("{} has started", program);

    if let Err(err) = run() {
        if let Some(app_err) = err.downcast_ref::<AppError>() {
            crate::errors::print_error(program, app_err);
            std::process::exit(app_err.exit_code());
        }

        // Set everything else as a I/O error for now. FIX LATER
        eprintln!("{program} ERROR: {:#}", err);
        std::process::exit(1);
    }

    info!("-------------------------------------------------------------------------------");
    info!("{} has completed successfully", program);
}
