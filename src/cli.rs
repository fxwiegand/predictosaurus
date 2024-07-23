use clap_derive::Parser;
use serde::{Deserialize, Deserializer};
use std::path::PathBuf;
use std::str::FromStr;

/// Uncertainty aware haplotype based genomic variant effect prediction
#[derive(Parser, Debug)]
#[clap(version, about)]
pub(crate) struct Predictosaurus {
    /// Path to the calls file
    #[clap(short, long)]
    pub(crate) calls: PathBuf,

    /// One or more observation files in the format `sample=observations.vcf`. Make sure the sample names match the sample names in the calls file.
    #[clap(short, long)]
    pub(crate) observations: Vec<ObservationFile>,

    /// Path to the gff file containing the features of interest.
    #[clap(short, long)]
    pub(crate) features: PathBuf,

    /// Path to reference genome fasta file
    #[clap(short, long)]
    pub(crate) reference: PathBuf,

    /// Path to the output file
    #[clap(long, default_value = ".")]
    pub(crate) output: PathBuf,
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
    use clap::Parser;

    #[test]
    fn parsing_calls_path() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
        ]);
        assert_eq!(args.calls, PathBuf::from("path/to/calls.vcf"));
    }

    #[test]
    fn parsing_single_observation_file() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
            "--observations",
            "sample1=observations1.vcf",
        ]);
        assert_eq!(args.observations.len(), 1);
        assert_eq!(args.observations[0].sample, "sample1");
        assert_eq!(
            args.observations[0].path,
            PathBuf::from("observations1.vcf")
        );
    }

    #[test]
    fn parsing_multiple_observation_files() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
            "--observations",
            "sample1=observations1.vcf",
            "--observations",
            "sample2=observations2.vcf",
        ]);
        assert_eq!(args.observations.len(), 2);
        assert_eq!(args.observations[0].sample, "sample1");
        assert_eq!(
            args.observations[0].path,
            PathBuf::from("observations1.vcf")
        );
        assert_eq!(args.observations[1].sample, "sample2");
        assert_eq!(
            args.observations[1].path,
            PathBuf::from("observations2.vcf")
        );
    }

    #[test]
    fn parsing_features_path() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
        ]);
        assert_eq!(args.features, PathBuf::from("path/to/features.gff"));
    }

    #[test]
    fn parsing_reference_path() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
        ]);
        assert_eq!(args.reference, PathBuf::from("path/to/reference.fasta"));
    }

    #[test]
    fn parsing_output_path() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
            "--output",
            "path/to/output",
        ]);
        assert_eq!(args.output, PathBuf::from("path/to/output"));
    }

    #[test]
    fn parsing_default_output_path() {
        let args = Predictosaurus::parse_from([
            "",
            "--calls",
            "path/to/calls.vcf",
            "--features",
            "path/to/features.gff",
            "--reference",
            "path/to/reference.fasta",
        ]);
        assert_eq!(args.output, PathBuf::from("."));
    }
}
