use clap_derive::{Parser, Subcommand};
use serde::{Deserialize, Deserializer};
use std::path::PathBuf;
use std::str::FromStr;

/// Uncertainty aware haplotype based genomic variant effect prediction
#[derive(Parser, Debug)]
#[clap(version, about)]
pub(crate) struct Predictosaurus {
    #[clap(subcommand)]
    pub(crate) command: Command,
}

#[derive(Subcommand, Debug)]
pub(crate) enum Command {
    /// Build a full variant graph out of VCF files and store it.
    Build {
        /// Path to the calls file
        #[clap(short, long)]
        calls: PathBuf,

        /// One or more observation files in the format `sample=observations.vcf`. Make sure the sample names match the sample names in the calls file.
        #[clap(short, long)]
        observations: Vec<ObservationFile>,

        /// Output path for the stored variant graphs
        #[clap(short, long)]
        output: PathBuf,
    },

    /// Retrieve subgraphs for individual features from the given GFF file
    Process {
        /// Path to the gff file containing the features of interest.
        #[clap(short, long)]
        features: PathBuf,

        /// Path to the output file
        #[clap(short, long)]
        output: PathBuf,
    },

    /// Filter paths with specific impacts or allele frequencies
    Filter {
        /// Path to the input file
        #[clap(short, long)]
        input: PathBuf,

        /// Path to reference genome fasta file
        #[clap(short, long)]
        reference: PathBuf,

        /// Path to the output Parquet file
        #[clap(short, long)]
        output: PathBuf,
    },

    /// Create visualizations and output HTML, TSV, or Vega specs
    Show {
        /// Path to the input data file
        #[clap(short, long)]
        input: PathBuf,

        /// Output format (html, tsv, vega)
        #[clap(short, long)]
        format: String,

        /// Path to the output file
        #[clap(short, long)]
        output: PathBuf,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct ObservationFile {
    pub(crate) path: PathBuf,
    pub(crate) sample: String,
}

impl FromStr for ObservationFile {
    type Err = String;

    fn from_str(string: &str) -> Result<ObservationFile, Self::Err> {
        let (sample, path) = string.split_once('=').expect("Invalid observation file parameter format. Make sure to use the format `--observations sample=observations.vcf`");
        Ok(ObservationFile {
            sample: sample.to_string(),
            path: PathBuf::from(path),
        })
    }
}

impl<'de> Deserialize<'de> for ObservationFile {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let string = String::deserialize(deserializer)?;
        Ok(ObservationFile::from_str(&string).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_str_parses_valid_observation_file() {
        let input = "sample1=observations.vcf";
        let observation_file = ObservationFile::from_str(input).unwrap();
        assert_eq!(observation_file.sample, "sample1");
        assert_eq!(observation_file.path, PathBuf::from("observations.vcf"));
    }

    #[test]
    #[should_panic]
    fn from_str_fails_on_invalid_format() {
        let input = "invalid_format";
        let result = ObservationFile::from_str(input);
    }
}
