use crate::cli::{Command, Format, Predictosaurus};
use crate::graph::duck::{create_paths, feature_graph, read_paths, write_graphs, write_paths};
use crate::graph::paths::{transcripts, Cds};
use crate::graph::VariantGraph;
use crate::show::{render_html_paths, render_tsv_paths, render_vl_paths};
use crate::utils::bcf::get_targets;
use crate::utils::create_output_dir;
use anyhow::{Context, Result};
use bio::bio_types::strand::Strand;
use bio::io::gff;
use bio::io::gff::GffType;
use clap::Parser;
use env_logger::Env;
use itertools::Itertools;
use log::{debug, info};
use rayon::prelude::*;
use std::collections::HashMap;

mod cli;
mod graph;
mod impact;
mod show;
mod transcription;
mod translation;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();
    let log_level = if args.verbose { "info" } else { "off" };
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();
    args.command.run()
}

impl Command {
    fn run(&self) -> Result<()> {
        match self {
            Command::Build {
                calls,
                observations,
                min_prob_present,
                output,
            } => {
                let targets = get_targets(calls)?;
                let mut graphs = HashMap::new();
                for target in targets {
                    info!("Building graph for target {}", target);
                    let variant_graph =
                        VariantGraph::build(calls, observations, &target, *min_prob_present)?;
                    info!(
                        "Finished building graph for target {} with {} nodes",
                        target,
                        variant_graph.graph.node_count()
                    );
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
                create_paths(output)?;
                info!("Reading reference genome from {:?}", reference);
                let reference_genome = utils::fasta::read_reference(reference);
                for transcript in transcripts(features)? {
                    info!("Processing transcript {}", transcript.name());
                    let weights = transcript.weights(graph, &reference_genome)?;
                    // -------- old ----------
                    if let Ok(graph) = feature_graph(
                        graph.to_owned(),
                        transcript.target.to_string(),
                        transcript.start()?,
                        transcript.end()?,
                    ) {
                        info!(
                            "Subgraph for transcript {} has {} nodes",
                            transcript.name(),
                            graph.graph.node_count()
                        );
                        let strand = record.strand().expect("Strand not found");
                        let phase: u8 = record.phase().clone().try_into().unwrap();
                        let weights = match strand {
                            Strand::Forward => Ok(graph
                                .paths()
                                .iter()
                                .map(|path| {
                                    path.weights(
                                        &graph,
                                        phase,
                                        reference_genome.get(&transcript.target).unwrap(),
                                        strand,
                                    )
                                    .unwrap()
                                })
                                .collect_vec()),
                            Strand::Reverse => Ok(graph
                                .reverse_paths()
                                .iter()
                                .map(|path| {
                                    path.weights(
                                        &graph,
                                        phase,
                                        &utils::fasta::reverse_complement(
                                            reference_genome.get(&transcript.target).unwrap(),
                                        ),
                                        strand,
                                    )
                                    .unwrap()
                                })
                                .collect_vec()),
                            Strand::Unknown => Err(anyhow::anyhow!(
                                "Strand is unknown for transcript {}",
                                transcript.name()
                            )),
                        };
                        let weights = weights?;
                        // ------------ old ------------
                    } else {
                        anyhow::bail!("No variant graph found for target {}", transcript.target);
                    }
                    info!(
                        "Writing {} paths for Transcript {}",
                        weights.len(),
                        transcript.name()
                    );
                    write_paths(output, weights, transcript)?;
                }
            }
            Command::Plot {
                input,
                format,
                output,
            } => {
                create_output_dir(output)?;
                let paths = read_paths(input)?;
                for (cds, paths) in paths {
                    match format {
                        Format::Html => {
                            render_html_paths(output, &paths, cds)?;
                        }
                        Format::Tsv => {
                            render_tsv_paths(output, &paths, cds)?;
                        }
                        Format::Vega => {
                            render_vl_paths(output, &paths, cds)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
