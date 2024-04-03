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

    #[test]
    fn test_observation_file_from_str() {
        let observation_file = ObservationFile::from_str("sample=observations.vcf").unwrap();
        assert_eq!(observation_file.sample, "sample");
        assert_eq!(observation_file.path, PathBuf::from("observations.vcf"));
    }
}
