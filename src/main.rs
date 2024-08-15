use crate::cli::{Command, Predictosaurus};
use crate::graph::{write_graphs, VariantGraph};
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
                graph,
                output,
            } => {
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;
                let variant_graphs: HashMap<String, VariantGraph> =
                    serde_json::from_reader(std::fs::File::open(graph)?)
                        .context("Failed to read graph file")?;
                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    let target = record.seqname().to_string();
                    if let Some(variant_graph) = variant_graphs.get(&target) {
                        let subgraph =
                            variant_graph.subgraph(*record.start() as i64, *record.end() as i64);
                        let output_file_path =
                            output.join(format!("{}.json", record.attributes().get("ID").unwrap()));
                        subgraph.write(&output_file_path)?;
                    } else {
                        anyhow::bail!("No variant graph found for target {}", target);
                    }
                }
            }
            Command::Filter {
                input,
                reference,
                output,
            } => {
                println!("Reading reference genome from {:?}", reference);
                let _reference_genome = utils::fasta::read_reference(reference);
                unimplemented!("Filter command not implemented")
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
