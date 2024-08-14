use crate::cli::{Command, Predictosaurus};
use crate::graph::VariantGraph;
use crate::utils::bcf::get_targets;
use anyhow::{Context, Result};
use bio::io::gff;
use bio::io::gff::GffType;
use clap::Parser;

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
                for target in targets {
                    let variant_graph = VariantGraph::build(calls, observations, &target)?;
                    variant_graph.write(&target, output)?;
                }
            }
            Command::Process { features, output } => {
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;
                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    println!(
                        "Feature: {} {}..{}",
                        record.seqname(),
                        record.start(),
                        record.end()
                    );
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
