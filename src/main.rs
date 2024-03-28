use crate::cli::Predictosaurus;
use crate::graph::VariantGraph;
use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use petgraph::dot::{Config, Dot};
use rust_htslib::bcf::{Read, Reader};

mod cli;
mod graph;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();

    let calls_file = args.calls;
    let observations_file = args.observations;
    let alignment_properties_file = args.alignment_properties;

    let _max_read_length =
        utils::AlignmentProperties::from_file(&alignment_properties_file)?.max_read_length;

    let variant_graph = VariantGraph::new(&calls_file, &observations_file)?;

    println!(
        "digraph {{ {:?} }}",
        Dot::with_config(&variant_graph.0, &[Config::GraphContentOnly])
    );

    Ok(())
}
