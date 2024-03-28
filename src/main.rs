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
    let alignment_properties_file = args.alignment_properties;
    let output_file = args.output;

    let max_read_length =
        utils::AlignmentProperties::from_file(&alignment_properties_file)?.max_read_length;

    VariantGraph::build(
        &calls_file,
        &observations_file,
        max_read_length,
        &output_file,
    )?;

    Ok(())
}
