use crate::graph::transcript::Transcript;
use crate::graph::EventProbs;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::dna_to_amino_acids;
use anyhow::Result;
use bio::stats::LogProb;
use itertools::Itertools;
use std::fmt::Display;
use std::path::PathBuf;

#[derive(Clone, Debug)]
pub(crate) struct Peptide {
    pub(crate) sequence: Vec<AminoAcid>,
    pub(crate) prob: EventProbs,
    pub(crate) vafs: Vec<f32>,
    pub(crate) transcript: Transcript,
}

impl Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.sequence
                .iter()
                .map(|a| a.short_abbreviation())
                .collect::<String>()
        )
    }
}

impl Peptide {
    pub(crate) fn from_rna(
        rna: Vec<u8>,
        prob: EventProbs,
        vafs: Vec<f32>,
        transcript: Transcript,
    ) -> Result<Self> {
        let sequence = dna_to_amino_acids(&rna)?;
        Ok(Peptide {
            sequence,
            prob,
            vafs,
            transcript,
        })
    }

    pub(crate) fn prob(&self, events: &Vec<String>) -> Result<LogProb> {
        let log_probs = events
            .iter()
            .map(|e| self.prob.prob(&e))
            .collect::<Result<Vec<LogProb>, _>>()?;
        Ok(LogProb::ln_sum_exp(&log_probs))
    }
}

pub(crate) fn write_peptides(peptides: Vec<Peptide>, output: &PathBuf) -> Result<()> {
    let mut fastq_writer = bio::io::fastq::Writer::to_file(output)?;
    let unique_peptides = peptides
        .iter()
        .map(|p| p.sequence.to_owned())
        .unique()
        .collect_vec();
    for peptide in unique_peptides {
        let record = peptide
            .iter()
            .map(|a| a.short_abbreviation())
            .collect::<String>();
        fastq_writer.write("", None, record.as_bytes(), b"")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::EventProbs;
    use bio::bio_types::strand::Strand;
    use bio::stats::LogProb;
    use bio::stats::Prob;
    use std::collections::HashMap;

    fn test_transcript() -> Transcript {
        Transcript::new(
            "ENSP007".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![],
        )
    }

    #[test]
    fn test_peptide_prob() {
        let probs = HashMap::from([
            ("A".to_string(), LogProb::from(Prob(0.1))),
            ("B".to_string(), LogProb::from(Prob(0.2))),
            ("C".to_string(), LogProb::from(Prob(0.3))),
        ]);
        let peptide = Peptide {
            sequence: vec![AminoAcid::Alanine, AminoAcid::Arginine],
            prob: EventProbs(probs),
            vafs: vec![],
            transcript: test_transcript(),
        };
        assert!(
            peptide
                .prob(&vec!["A".to_string(), "B".to_string()])
                .unwrap()
                < LogProb::from(Prob(0.31))
                && peptide
                    .prob(&vec!["A".to_string(), "B".to_string()])
                    .unwrap()
                    > LogProb::from(Prob(0.29)),
        );
        assert!(
            peptide
                .prob(&vec!["A".to_string(), "B".to_string(), "C".to_string()])
                .unwrap()
                < LogProb::from(Prob(0.61))
                && peptide
                    .prob(&vec!["A".to_string(), "B".to_string(), "C".to_string()])
                    .unwrap()
                    > LogProb::from(Prob(0.59)),
        );
    }

    #[test]
    fn write_peptides_creates_correct_output() {
        let peptides = vec![
            Peptide {
                sequence: vec![AminoAcid::Alanine, AminoAcid::Arginine],
                prob: EventProbs(HashMap::new()),
                vafs: vec![],
                transcript: test_transcript(),
            },
            Peptide {
                sequence: vec![AminoAcid::Alanine, AminoAcid::Arginine],
                prob: EventProbs(HashMap::new()),
                vafs: vec![],
                transcript: test_transcript(),
            },
            Peptide {
                sequence: vec![AminoAcid::Cysteine],
                prob: EventProbs(HashMap::new()),
                vafs: vec![],
                transcript: test_transcript(),
            },
        ];
        let tmp = tempfile::tempdir().unwrap();
        let output_path = tmp.path().join("output.fastq");
        write_peptides(peptides, &output_path).unwrap();
        let contents = std::fs::read_to_string(output_path).unwrap();
        assert!(contents.contains("AR"));
        assert!(contents.contains("C"));
    }
}
