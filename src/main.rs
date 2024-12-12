use crate::cli::{Command, Format, Predictosaurus};
use crate::graph::duck::{create_paths, feature_graph, read_paths, write_graphs, write_paths};
use crate::graph::VariantGraph;
use crate::show::{render_html_paths, render_tsv_paths, render_vl_paths};
use crate::utils::bcf::get_targets;
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
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;

                info!("Reading reference genome from {:?}", reference);
                let reference_genome = utils::fasta::read_reference(reference);
                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    let target = record.seqname().to_string();
                    let cds_id = record
                        .attributes()
                        .get("ID")
                        .expect(
                            format!("No ID found for CDS in sequence {}", record.seqname())
                                .as_str(),
                        )
                        .to_string();
                    info!("Processing CDS {} in sequence {}", cds_id, target);
                    if let Ok(graph) = feature_graph(
                        graph.to_owned(),
                        target.to_string(),
                        *record.start(),
                        *record.end(),
                    ) {
                        info!(
                            "Subgraph for CDS {} with target {} has {} nodes",
                            cds_id,
                            target,
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
                                        reference_genome.get(&target).unwrap(),
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
                                        reference_genome.get(&target).unwrap(),
                                        strand,
                                    )
                                    .unwrap()
                                })
                                .collect_vec()),
                            Strand::Unknown => Err(anyhow::anyhow!(
                                "Strand is unknown for sequence {}",
                                record.seqname()
                            )),
                        };
                        write_paths(output, weights?, target, cds_id)?;
                    } else {
                        anyhow::bail!("No variant graph found for target {}", target);
                    }
                }
            }
            Command::Plot {
                input,
                format,
                output,
            } => {
                let paths = read_paths(input)?;
                for (feature, paths) in paths {
                    match format {
                        Format::Html => {
                            render_html_paths(output, &paths, feature)?;
                        }
                        Format::Tsv => {
                            render_tsv_paths(output, &paths, feature)?;
                        }
                        Format::Vega => {
                            render_vl_paths(output, &paths, feature)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
