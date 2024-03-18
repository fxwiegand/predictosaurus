use crate::cli::Predictosaurus;
use crate::utils::bcf::extract_event_names;
use anyhow::Result;
use clap::Parser;
use rust_htslib::bcf::{Read, Reader};

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
        let calls_record = calls_record?;
        let observations_record = observations_record?;

        // Process records from both files together
        println!("Calls record: {:?}", calls_record);
        println!("Observations record: {:?}", observations_record);
    }

    Ok(())
}
