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
use itertools::Itertools;
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
                let reference_genome = utils::fasta::read_reference(reference);

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
                        create_paths(output)?;
                        write_paths(output, weights?, target, cds_id)?;
                    } else {
                        anyhow::bail!("No variant graph found for target {}", target);
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
            Command::Show {
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
