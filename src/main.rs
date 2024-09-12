use crate::cli::{Command, Predictosaurus};
use crate::graph::duck::{feature_graph, write_graphs};
use crate::graph::VariantGraph;
use crate::utils::bcf::get_targets;
use anyhow::{Context, Result};
use bio::io::gff;
use bio::io::gff::GffType;
use clap::Parser;
use std::collections::HashMap;

mod cli;
mod graph;
mod impact;
mod transcription;
mod translation;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();
    args.command.run()
}

impl Command {
    fn run(&self) -> Result<()> {
        match self {
            Command::Build {
                calls,
                observations,
                output,
            } => {
                utils::create_output_dir(output)?;
                let targets = get_targets(calls)?;
                let mut graphs = HashMap::new();
                for target in targets {
                    let variant_graph = VariantGraph::build(calls, observations, &target)?;
                    if !variant_graph.is_empty() {
                        graphs.insert(target, variant_graph);
                    }
                }
                write_graphs(graphs, output)?;
            }
            Command::Process {
                features,
                reference,
                graph,
                output,
            } => {
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;

                println!("Reading reference genome from {:?}", reference);
                let _reference_genome = utils::fasta::read_reference(reference);

                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    let target = record.seqname().to_string();
                    if let Ok(graph) = feature_graph(graph.to_owned(), target.to_string(), *record.start(), *record.end()) {
                        // TODO: Implement filtering of nodes and edges based on feature coordinates to create subgraph.
                        // Create Paths and their impacts for subgraph/feature and serialize them for usage in show command.
                    } else {
                        anyhow::bail!("No variant graph found for target {}", target);
                    }
                }
            }
            Command::Show {
                input,
                format,
                output,
            } => {
                unimplemented!("Show command not implemented")
            }
        }
        Ok(())
    }
}
