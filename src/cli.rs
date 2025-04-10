use clap_derive::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Deserializer};
use std::fmt::Display;
use std::path::PathBuf;
use std::str::FromStr;

/// Uncertainty aware haplotype based genomic variant effect prediction
#[derive(Parser, Debug)]
#[clap(version, about)]
pub(crate) struct Predictosaurus {
    #[clap(subcommand)]
    pub(crate) command: Command,

    #[clap(short, long, global = true)]
    pub(crate) verbose: bool,
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

        /// Minimum probability for a variant to be considered in the graph
        #[clap(short, long, default_value = "0.8")]
        min_prob_present: f64,

        /// Path to the output file containing the impact graph
        #[clap(long)]
        output: PathBuf,
    },

    /// Retrieve subgraphs for individual features from the given GFF file
    Process {
        /// Path to the gff file containing the features of interest.
        #[clap(short, long)]
        features: PathBuf,

        /// Path to reference genome fasta file
        #[clap(short, long)]
        reference: PathBuf,

        /// Path to the graph file generated by the build command
        #[clap(short, long)]
        graph: PathBuf,

        /// Path to the output file containing the paths for the features given via the GFF file
        #[clap(short, long)]
        output: PathBuf,
    },
    /// Output all distinct peptides from the given features to a fastq file per given CDS in the feature file
    Peptides {
        /// Path to the gff file containing the features of interest.
        #[clap(short, long)]
        features: PathBuf,

        /// Path to reference genome fasta file
        #[clap(short, long)]
        reference: PathBuf,

        /// Path to the graph file generated by the build command
        #[clap(short, long)]
        graph: PathBuf,

        /// Interval for peptide lengths to be generated in the format `start-end`
        #[clap(short, long, default_value_t = Interval::default())]
        interval: Interval,

        /// Sample name used to retrieve allele frequencies from the graph
        #[clap(short, long)]
        sample: String,

        #[clap(short, long)]
        events: Vec<String>,

        #[clap(long)]
        min_event_prob: f64,

        #[clap(short, long)]
        background_events: Vec<String>,

        #[clap(long)]
        max_background_event_prob: f64,

        /// Path to the output directory for the fastq files
        #[clap(short, long)]
        output: PathBuf,
    },
    /// Create visualizations and output HTML, TSV, or Vega specs
    Plot {
        /// Path to the input data file
        #[clap(short, long)]
        input: PathBuf,

        /// Output format (html, tsv, vega)
        #[clap(short, long)]
        format: Format,

        /// Path to the output directory
        #[clap(short, long)]
        output: PathBuf,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct Interval {
    pub(crate) start: u32,
    pub(crate) end: u32,
}

impl Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

impl Default for Interval {
    fn default() -> Self {
        Interval { start: 8, end: 11 }
    }
}

impl Iterator for Interval {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start <= self.end {
            let current = self.start;
            self.start += 1;
            Some(current)
        } else {
            None
        }
    }
}

impl FromStr for Interval {
    type Err = String;

    fn from_str(string: &str) -> Result<Interval, Self::Err> {
        let (start, end) = string
            .split_once('-')
            .expect("Invalid interval format. Make sure to use the format `start-end`");
        let start = start.parse::<u32>().unwrap();
        let end = end.parse::<u32>().unwrap();
        if start > end {
            panic!("Invalid interval format. Make sure to use the format `start-end`");
        }
        Ok(Interval { start, end })
    }
}

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum Format {
    Html,
    Tsv,
    Vega,
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

    #[test]
    fn from_str_parses_valid_interval() {
        let input = "10-20";
        let interval = Interval::from_str(input).unwrap();
        assert_eq!(interval.start, 10);
        assert_eq!(interval.end, 20);
    }

    #[test]
    #[should_panic]
    fn from_str_fails_on_invalid_interval() {
        let input = "10:20";
        let result = Interval::from_str(input);
    }

    #[test]
    #[should_panic]
    fn from_str_fails_when_start_greater_than_end() {
        let input = "20-10";
        let result = Interval::from_str(input);
    }

    #[test]
    fn iterator_yields_all_values_in_range() {
        let mut interval = Interval { start: 1, end: 3 };
        assert_eq!(interval.next(), Some(1));
        assert_eq!(interval.next(), Some(2));
        assert_eq!(interval.next(), Some(3));
        assert_eq!(interval.next(), None);
    }

    #[test]
    fn default_interval_has_correct_start_and_end() {
        let interval = Interval::default();
        assert_eq!(interval.start, 8);
        assert_eq!(interval.end, 11);
    }

    #[test]
    fn display_formats_interval_correctly() {
        let interval = Interval { start: 5, end: 10 };
        assert_eq!(format!("{}", interval), "5-10");
    }
}
