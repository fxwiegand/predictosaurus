use crate::cli::Predictosaurus;
use crate::graph::VariantGraph;
use anyhow::Result;
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

    let calls_file = args.calls;
    let observation_files = args.observations;
    let features_file = args.features;
    let output_file = args.output;

    utils::create_output_dir(&output_file)?;

    println!("Reading reference genome from {:?}", args.reference);
    let reference_genome = utils::fasta::read_reference(&args.reference);

    let mut feature_reader = gff::Reader::from_file(features_file, GffType::GFF3).unwrap();
    for record in feature_reader
        .records()
        .filter_map(Result::ok)
        .filter(|record| record.feature_type() == "CDS")
    {
        let start = *record.start() as i64;
        let end = *record.end() as i64;
        println!(
            "Building variant graph for CDS at {}:{}-{}",
            record.seqname(),
            start,
            end
        );
        let variant_graph = VariantGraph::build(
            &calls_file,
            &observation_files,
            record.seqname(),
            start,
            end,
        )?;
        let phase: u8 = record.phase().clone().try_into().expect("Invalid phase");
        let strand = record.strand().expect("Strand not found");
        let forward_seq = reference_genome.get(record.seqname()).expect(
            format!(
                "Reference sequence {} not found in provided FASTA file",
                record.seqname()
            )
            .as_str(),
        );
        let ref_seq = match strand {
            Strand::Forward => forward_seq,
            Strand::Reverse => &{ utils::fasta::reverse_complement(forward_seq) },
            Strand::Unknown => {
                panic!("Strand is unknown")
            }
        };

        let paths = variant_graph.paths();
        for path in paths {
            println!(
                "Path with impact {} and weight {}",
                path.impact(&variant_graph, phase, &ref_seq,).unwrap(),
                path.weight(&variant_graph)
            );
            println!();
        }
    }
    Ok(())
}
