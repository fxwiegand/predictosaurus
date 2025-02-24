use crate::cli::{Command, Format, Predictosaurus};
use crate::graph::duck::{create_paths, feature_graph, read_paths, write_graphs, write_paths};
use crate::graph::paths::Cds;
use crate::graph::peptide::write_peptides;
use crate::graph::transcript::transcripts;
use crate::graph::VariantGraph;
use crate::show::{render_html_paths, render_tsv_paths, render_vl_paths};
use crate::utils::bcf::get_targets;
use crate::utils::create_output_dir;
use anyhow::{Context, Result};
use bio::bio_types::strand::Strand;
use bio::io::gff;
use bio::io::gff::GffType;
use bio::stats::{LogProb, Prob};
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
                    let variant_graph = VariantGraph::build(
                        calls,
                        observations,
                        &target,
                        LogProb::from(Prob(*min_prob_present)),
                    )?;
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
                    info!(
                        "Writing {} paths for Transcript {}",
                        weights.len(),
                        transcript.name()
                    );
                    write_paths(output, weights, transcript)?;
                }
            }
            Command::Peptides {
                features,
                reference,
                graph,
                interval,
                sample,
                output,
                events,
                min_event_prob,
                background_events,
                max_background_event_prob,
            } => {
                info!("Reading reference genome from {:?}", reference);
                let reference_genome = utils::fasta::read_reference(reference);
                let mut peptides = Vec::new();
                for transcript in transcripts(features)? {
                    info!("Processing transcript {}", transcript.name());
                    let transcript_peptides = transcript.peptides(
                        graph,
                        &reference_genome,
                        interval.clone(),
                        &sample,
                        events,
                        LogProb::from(Prob(*min_event_prob)),
                        background_events,
                        LogProb::from(Prob(*max_background_event_prob)),
                    )?;
                    peptides.extend(transcript_peptides);
                }
                info!("Writing peptides to {:?}", output);
                write_peptides(peptides, output)?;
            }
            Command::Plot {
                input,
                format,
                output,
            } => {
                create_output_dir(output)?;
                let paths = read_paths(input)?;
                for (transcript, paths) in paths {
                    match format {
                        Format::Html => {
                            render_html_paths(output, &paths, transcript)?;
                        }
                        Format::Tsv => {
                            render_tsv_paths(output, &paths, transcript)?;
                        }
                        Format::Vega => {
                            render_vl_paths(output, &paths, transcript)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
