use clap::Parser;
use std::path::PathBuf;

/// Uncertainty aware haplotype based genomic variant effect prediction
#[derive(Parser, Debug)]
#[clap(version, about)]
pub(crate) struct Predictosaurus {
    /// Path to the calls file
    #[clap(short, long)]
    pub(crate) calls: PathBuf,

    /// Path to the observations file
    #[clap(short, long)]
    pub(crate) observations: PathBuf,

    /// Path to the alignment properties file
    #[clap(short, long)]
    pub(crate) alignment_properties: PathBuf,
}
