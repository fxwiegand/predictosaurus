use crate::cli::{Command, Format, Predictosaurus};
use crate::graph::duck::{create_paths, feature_graph, read_paths, write_graphs, write_paths};
use crate::graph::paths::Cds;
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
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;

                info!("Reading reference genome from {:?}", reference);
                let reference_genome = utils::fasta::read_reference(reference);
                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    let cds = Cds::from_record(&record)?;
                    info!("Processing CDS {} in sequence {}", cds.feature, cds.target);
                    if let Ok(graph) = feature_graph(
                        graph.to_owned(),
                        cds.target.to_string(),
                        *record.start(),
                        *record.end(),
                    ) {
                        info!(
                            "Subgraph for CDS {} with target {} has {} nodes",
                            cds.feature,
                            cds.target,
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
                                        reference_genome.get(&cds.target).unwrap(),
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
                                        reference_genome.get(&cds.target).unwrap(),
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
                        let weights = weights?;
                        info!("Writing {} paths for CDS {}", weights.len(), cds.name());
                        write_paths(output, weights, cds)?;
                    } else {
                        anyhow::bail!("No variant graph found for target {}", cds.target);
                    }
                }
            }
            Command::Peptides {
                features,
                reference,
                graph,
                interval,
                output,
            } => {
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;
                let reference_genome = utils::fasta::read_reference(reference);
                let mut peptides = Vec::new();

                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    let target = record.seqname().to_string();
                    let cds_id = record.attributes().get("ID").unwrap().to_string();
                    if let Ok(graph) = feature_graph(
                        graph.to_owned(),
                        target.to_string(),
                        *record.start(),
                        *record.end(),
                    ) {
                        let strand = record.strand().expect("Strand not found");
                        let phase: u8 = record.phase().clone().try_into().unwrap();
                        let proteins = match strand {
                            Strand::Forward => Ok(graph
                                .paths()
                                .iter()
                                .map(|path| {
                                    path.protein(
                                        &graph,
                                        phase,
                                        reference_genome.get(&target).unwrap(),
                                        strand,
                                        *record.start() as usize,
                                        *record.end() as usize,
                                    )
                                    .unwrap()
                                })
                                .collect_vec()),
                            Strand::Reverse => Ok(graph
                                .reverse_paths()
                                .iter()
                                .map(|path| {
                                    path.protein(
                                        &graph,
                                        phase,
                                        reference_genome.get(&target).unwrap(),
                                        strand,
                                        *record.start() as usize,
                                        *record.end() as usize,
                                    )
                                    .unwrap()
                                })
                                .collect_vec()),
                            Strand::Unknown => Err(anyhow::anyhow!(
                                "Strand is unknown for sequence {}",
                                record.seqname()
                            )),
                        };
                        proteins?.iter().for_each(|protein| {
                            peptides.push(protein.peptides(interval.to_owned()));
                        });
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
