use crate::cli::Predictosaurus;
use crate::graph::VariantGraph;
use anyhow::Result;
use clap::Parser;
mod cli;
mod graph;
mod transcription;
mod translation;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();

    let calls_file = args.calls;
    let observation_files = args.observations;
    let output_file = args.output;

    utils::create_output_dir(&output_file)?;

    VariantGraph::build(&calls_file, &observation_files, &output_file)?;

    Ok(())
}
