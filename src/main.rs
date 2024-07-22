use crate::cli::Predictosaurus;
use crate::graph::VariantGraph;
use anyhow::Result;
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

    let calls_file = args.calls;
    let observation_files = args.observations;
    let features_file = args.features;
    let output_file = args.output;

    utils::create_output_dir(&output_file)?;

    // TODO: Load reference genome into memory and pass slice to impact calculation

    let mut feature_reader = gff::Reader::from_file(features_file, GffType::GFF3).unwrap();
    for record in feature_reader
        .records()
        .filter_map(Result::ok)
        .filter(|record| record.feature_type() == "CDS")
    {
        let variant_graph = VariantGraph::build(
            &calls_file,
            &observation_files,
            record.seqname(),
            *record.start() as i64,
            (record.end() - record.start()) as i64,
        )?;
        let phase = record.phase().as_u8();
        println!("{:?}", record);
        if let Some(phase) = phase {
            let paths = variant_graph.paths();
            for path in paths {
                println!("{:?}", path.impact(&variant_graph, phase, &[])?);
            }
        } else {
            println!("No phase found for CDS record");
        }
    }
    Ok(())
}
