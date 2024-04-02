use crate::cli::Predictosaurus;
use crate::graph::VariantGraph;
use anyhow::Result;
use clap::Parser;
mod cli;
mod graph;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();

    let calls_file = args.calls;
    let observations_file = args.observations;
    let output_file = args.output;

    VariantGraph::build(&calls_file, &observations_file, &output_file)?;

    Ok(())
}
