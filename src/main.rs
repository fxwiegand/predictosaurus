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
    let features_file = args.features;
    let output_file = args.output;

    utils::create_output_dir(&output_file)?;

    let variant_graph = VariantGraph::build(&calls_file, &observation_files)?;
    variant_graph.to_file(&output_file)?;

    // TODO: Parse the features file and iter transcripts per feature. For each transcript go through paths of the graph and translate each codon around a variant to an amino acid and calculate impact keeping track of the frameshift of the path.

    Ok(())
}
