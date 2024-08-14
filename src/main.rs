use crate::cli::{Command, Predictosaurus};
use crate::graph::VariantGraph;
use crate::utils::bcf::get_targets;
use anyhow::{anyhow, Context, Result};
use bio::bio_types::strand::Strand;
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
                utils::create_output_dir(&output)?;
                let targets = get_targets(&calls)?;
                for target in targets {
                    let variant_graph = VariantGraph::build(&calls, &observations, &target)?;
                    variant_graph.write(&output)?;
                }
                unimplemented!("Build command not implemented")
            }
            Command::Process { features, output } => {
                let mut feature_reader = gff::Reader::from_file(features, GffType::GFF3)
                    .context("Failed to open GFF file")?;
                for record in feature_reader
                    .records()
                    .filter_map(Result::ok)
                    .filter(|record| record.feature_type() == "CDS")
                {
                    unimplemented!("Process command not implemented")
                }
            }
            Command::Filter {
                input,
                reference,
                output,
            } => {
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
