use crate::cli::{Command, Format, Predictosaurus};
use crate::graph::duck::{
    create_paths, create_scores, feature_graph, read_scores, write_graphs, write_scores,
};
use crate::graph::paths::Cds;
use crate::graph::peptide::write_peptides;
use crate::graph::transcript::transcripts;
use crate::graph::VariantGraph;
use crate::show::{render_html_paths, render_scores, render_tsv_paths, render_vl_paths};
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
use rayon::ThreadPoolBuilder;
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
    if let Some(n) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to build thread pool");
    }
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
                let mut graphs = targets
                    .into_par_iter()
                    .filter_map(|target| {
                        info!("Building graph for target {target}");
                        let variant_graph = VariantGraph::build(
                            calls,
                            observations,
                            &target,
                            LogProb::from(Prob(*min_prob_present)),
                        )
                        .ok()?;
                        info!(
                            "Finished building graph for target {} with {} nodes",
                            target,
                            variant_graph.graph.node_count()
                        );
                        if !variant_graph.is_empty() {
                            Some((target, variant_graph))
                        } else {
                            None
                        }
                    })
                    .collect();
                write_graphs(graphs, output)?;
            }
            Command::Process {
                features,
                reference,
                graph,
                output,
            } => {
                create_scores(output)?;
                info!("Reading reference genome from {reference:?}");
                let reference_genome = utils::fasta::read_reference(reference);
                for transcript in transcripts(features, graph)? {
                    info!("Processing transcript {}", transcript.name());
                    let scores = transcript.scores(graph, &reference_genome)?;
                    info!(
                        "Writing scores for {} different haplotypes for transcript {}",
                        scores.len(),
                        transcript.name()
                    );
                    write_scores(output, &scores, transcript)?;
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
                info!("Reading reference genome from {reference:?}");
                let reference_genome = utils::fasta::read_reference(reference);
                let mut peptides = Vec::new();
                for transcript in transcripts(features, graph)? {
                    info!("Processing transcript {}", transcript.name());
                    let transcript_peptides = transcript.peptides(
                        graph,
                        &reference_genome,
                        interval.clone(),
                        sample,
                        events,
                        LogProb::from(Prob(*min_event_prob)),
                        background_events,
                        LogProb::from(Prob(*max_background_event_prob)),
                    )?;
                    peptides.extend(transcript_peptides);
                }
                info!("Writing peptides to {output:?}");
                write_peptides(peptides, output)?;
            }
            Command::Plot { input, output } => {
                if let Some(parent) = output.parent() {
                    create_output_dir(parent)?;
                }
                let scores = read_scores(input)?;
                render_scores(output, &scores)?;
            }
        }
        Ok(())
    }
}
