use crate::cli::Predictosaurus;
use crate::utils::bcf::extract_event_names;
use anyhow::Result;
use clap::Parser;
use rust_htslib::bcf::{Read, Reader};
use varlociraptor::calling::variants::preprocessing::read_observations;
use varlociraptor::utils::collect_variants::collect_variants;

mod cli;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();

    let calls_file = args.calls;
    let observations_file = args.observations;

    let mut calls_reader = Reader::from_path(&calls_file)?;
    let mut observations_reader = Reader::from_path(observations_file)?;

    let _event_names = extract_event_names(&calls_file);

    for (calls_record, observations_record) in
        calls_reader.records().zip(observations_reader.records())
    {
        let mut calls_record = calls_record?;
        let mut observations_record = observations_record?;

        let variants = collect_variants(&mut calls_record, false, None)?;
        // let observations = read_observations(&mut observations_record)?;

        dbg!(&variants);
        // dbg!(observations.pileup);
    }

    Ok(())
}
