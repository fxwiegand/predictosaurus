use crate::cli::Interval;
use crate::graph::node::{Node, NodeType};
use crate::graph::transcript::Transcript;
use crate::translation::distance::DistanceMetric;
use crate::translation::dna_to_amino_acids;
use crate::utils::fasta::reverse_complement;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use std::collections::HashMap;
use std::fmt::Display;

/// A protein consisting of a sequence of amino acids
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Protein {
    pub(crate) sequence: Vec<AminoAcid>,
}

impl Protein {
    /// Creates a new protein from a sequence of amino acids
    pub(crate) fn new(sequence: Vec<AminoAcid>) -> Protein {
        Protein { sequence }
    }

    /// Creates all peptides of length k from the protein
    pub(crate) fn k_peptides(&self, k: u32) -> Vec<Vec<AminoAcid>> {
        self.sequence
            .windows(k as usize)
            .map(|window| window.to_vec())
            .unique()
            .collect()
    }

    pub(crate) fn peptides(&self, interval: Interval) -> Vec<Vec<AminoAcid>> {
        interval
            .into_iter()
            .flat_map(|i| self.k_peptides(i))
            .collect()
    }

    pub(crate) fn from_transcript(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
    ) -> Result<Protein> {
        let target_ref = reference
            .get(&transcript.target)
            .ok_or_else(|| anyhow::anyhow!("Missing reference for {}", transcript.target))?;

        let cds_seq: Vec<u8> = transcript
            .cds()
            .flat_map(|cds| {
                let region = &target_ref[cds.start as usize..cds.end as usize];
                let phase = cds.phase as usize;
                match transcript.strand {
                    Strand::Reverse => reverse_complement(&region[..region.len() - phase]),
                    _ => region[phase..].to_vec(),
                }
            })
            .collect();

        Ok(Protein::new(dna_to_amino_acids(&cds_seq)?))
    }

    pub(crate) fn from_haplotype(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
        haplotype: &[Node],
    ) -> Result<Protein> {
        let target_ref = reference
            .get(&transcript.target)
            .ok_or_else(|| anyhow::anyhow!("Missing reference for {}", transcript.target))?;

        let cds_seq: Vec<u8> = transcript
            .cds()
            .flat_map(|cds| {
                let mut region = target_ref[cds.start as usize..cds.end as usize].to_vec();
                let mut offset: isize = 0;

                for node in haplotype.iter().filter(|n| {
                    n.pos >= cds.start as i64 && n.pos < cds.end as i64 && n.node_type.is_variant()
                }) {
                    let pos = (node.pos - cds.start as i64) as usize + offset as usize;
                    region.splice(
                        pos..pos + node.reference_allele.len(),
                        node.alternative_allele.bytes(),
                    );
                    offset += node.frameshift() as isize;
                }

                let phase = cds.phase as usize;
                match transcript.strand {
                    Strand::Reverse => reverse_complement(&region[..region.len() - phase]),
                    _ => region[phase..].to_vec(),
                }
            })
            .collect();

        Ok(Protein::new(dna_to_amino_acids(&cds_seq)?))
    }
}

impl Display for Protein {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.sequence
                .iter()
                .map(|amino_acid| amino_acid.abbreviation())
                .collect::<Vec<&str>>()
                .join("")
        )
    }
}

/// Amino acids and stop
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum AminoAcid {
    Alanine,
    Arginine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    GlutamicAcid,
    Glutamine,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
    Stop,
}

impl AminoAcid {
    /// Converts a codon to an amino acid or stop codon. Returns an error if the codon is invalid.
    pub(crate) fn from_codon(codon: &[u8]) -> Result<AminoAcid> {
        match codon {
            b"GCU" | b"GCC" | b"GCA" | b"GCG" => Ok(AminoAcid::Alanine),
            b"CGU" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => Ok(AminoAcid::Arginine),
            b"AAU" | b"AAC" => Ok(AminoAcid::Asparagine),
            b"GAU" | b"GAC" => Ok(AminoAcid::AsparticAcid),
            b"UGU" | b"UGC" => Ok(AminoAcid::Cysteine),
            b"GAA" | b"GAG" => Ok(AminoAcid::GlutamicAcid),
            b"CAA" | b"CAG" => Ok(AminoAcid::Glutamine),
            b"GGU" | b"GGC" | b"GGA" | b"GGG" => Ok(AminoAcid::Glycine),
            b"CAU" | b"CAC" => Ok(AminoAcid::Histidine),
            b"AUU" | b"AUC" | b"AUA" => Ok(AminoAcid::Isoleucine),
            b"UUA" | b"UUG" | b"CUU" | b"CUC" | b"CUA" | b"CUG" => Ok(AminoAcid::Leucine),
            b"AAA" | b"AAG" => Ok(AminoAcid::Lysine),
            b"AUG" => Ok(AminoAcid::Methionine),
            b"UUU" | b"UUC" => Ok(AminoAcid::Phenylalanine),
            b"CCU" | b"CCC" | b"CCA" | b"CCG" => Ok(AminoAcid::Proline),
            b"UCU" | b"UCC" | b"UCA" | b"UCG" | b"AGU" | b"AGC" => Ok(AminoAcid::Serine),
            b"ACU" | b"ACC" | b"ACA" | b"ACG" => Ok(AminoAcid::Threonine),
            b"UGG" => Ok(AminoAcid::Tryptophan),
            b"UAU" | b"UAC" => Ok(AminoAcid::Tyrosine),
            b"GUU" | b"GUC" | b"GUA" | b"GUG" => Ok(AminoAcid::Valine),
            b"UAA" | b"UAG" | b"UGA" => Ok(AminoAcid::Stop),
            _ => anyhow::bail!("Invalid codon"),
        }
    }

    /// Returns the 3 letter abbreviation of the amino acid
    fn abbreviation(&self) -> &str {
        match self {
            AminoAcid::Alanine => "Ala",
            AminoAcid::Arginine => "Arg",
            AminoAcid::Asparagine => "Asn",
            AminoAcid::AsparticAcid => "Asp",
            AminoAcid::Cysteine => "Cys",
            AminoAcid::GlutamicAcid => "Glu",
            AminoAcid::Glutamine => "Gln",
            AminoAcid::Glycine => "Gly",
            AminoAcid::Histidine => "His",
            AminoAcid::Isoleucine => "Ile",
            AminoAcid::Leucine => "Leu",
            AminoAcid::Lysine => "Lys",
            AminoAcid::Methionine => "Met",
            AminoAcid::Phenylalanine => "Phe",
            AminoAcid::Proline => "Pro",
            AminoAcid::Serine => "Ser",
            AminoAcid::Threonine => "Thr",
            AminoAcid::Tryptophan => "Trp",
            AminoAcid::Tyrosine => "Tyr",
            AminoAcid::Valine => "Val",
            AminoAcid::Stop => "Stop",
        }
    }

    pub(crate) fn short_abbreviation(&self) -> char {
        match self {
            AminoAcid::Alanine => 'A',
            AminoAcid::Arginine => 'R',
            AminoAcid::Asparagine => 'N',
            AminoAcid::AsparticAcid => 'D',
            AminoAcid::Cysteine => 'C',
            AminoAcid::GlutamicAcid => 'E',
            AminoAcid::Glutamine => 'Q',
            AminoAcid::Glycine => 'G',
            AminoAcid::Histidine => 'H',
            AminoAcid::Isoleucine => 'I',
            AminoAcid::Leucine => 'L',
            AminoAcid::Lysine => 'K',
            AminoAcid::Methionine => 'M',
            AminoAcid::Phenylalanine => 'F',
            AminoAcid::Proline => 'P',
            AminoAcid::Serine => 'S',
            AminoAcid::Threonine => 'T',
            AminoAcid::Tryptophan => 'W',
            AminoAcid::Tyrosine => 'Y',
            AminoAcid::Valine => 'V',
            AminoAcid::Stop => 'X',
        }
    }

    /// Returns true if the amino acid is a stop codon
    pub(crate) fn is_stop(&self) -> bool {
        matches!(self, AminoAcid::Stop)
    }

    /// Returns the distance between two amino acids using the specified metric
    pub(crate) fn distance(&self, other: &AminoAcid, metric: DistanceMetric) -> f64 {
        metric.compute(self, other)
    }
}

impl Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.abbreviation())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_codons() {
        assert!(matches!(
            AminoAcid::from_codon(b"GCU"),
            Ok(AminoAcid::Alanine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CGU"),
            Ok(AminoAcid::Arginine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AAU"),
            Ok(AminoAcid::Asparagine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GAU"),
            Ok(AminoAcid::AsparticAcid)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UGU"),
            Ok(AminoAcid::Cysteine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GAA"),
            Ok(AminoAcid::GlutamicAcid)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CAA"),
            Ok(AminoAcid::Glutamine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GGU"),
            Ok(AminoAcid::Glycine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CAU"),
            Ok(AminoAcid::Histidine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AUU"),
            Ok(AminoAcid::Isoleucine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UUA"),
            Ok(AminoAcid::Leucine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AAA"),
            Ok(AminoAcid::Lysine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AUG"),
            Ok(AminoAcid::Methionine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UUU"),
            Ok(AminoAcid::Phenylalanine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CCU"),
            Ok(AminoAcid::Proline)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UCU"),
            Ok(AminoAcid::Serine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"ACU"),
            Ok(AminoAcid::Threonine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UGG"),
            Ok(AminoAcid::Tryptophan)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UAU"),
            Ok(AminoAcid::Tyrosine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GUU"),
            Ok(AminoAcid::Valine)
        ));
        assert!(matches!(AminoAcid::from_codon(b"UAA"), Ok(AminoAcid::Stop)));
        assert!(matches!(AminoAcid::from_codon(b"UAG"), Ok(AminoAcid::Stop)));
        assert!(matches!(AminoAcid::from_codon(b"UGA"), Ok(AminoAcid::Stop)));
    }

    #[test]
    fn test_invalid_codons() {
        assert!(AminoAcid::from_codon(b"XXX").is_err());
        assert!(AminoAcid::from_codon(b"XYZ").is_err());
        assert!(AminoAcid::from_codon(b"UGX").is_err());
        assert!(AminoAcid::from_codon(b"ACX").is_err());
    }

    #[test]
    fn test_abbreviations() {
        assert_eq!(AminoAcid::Alanine.abbreviation(), "Ala");
        assert_eq!(AminoAcid::Arginine.abbreviation(), "Arg");
        assert_eq!(AminoAcid::Asparagine.abbreviation(), "Asn");
        assert_eq!(AminoAcid::AsparticAcid.abbreviation(), "Asp");
        assert_eq!(AminoAcid::Cysteine.abbreviation(), "Cys");
        assert_eq!(AminoAcid::GlutamicAcid.abbreviation(), "Glu");
        assert_eq!(AminoAcid::Glutamine.abbreviation(), "Gln");
        assert_eq!(AminoAcid::Glycine.abbreviation(), "Gly");
        assert_eq!(AminoAcid::Histidine.abbreviation(), "His");
        assert_eq!(AminoAcid::Isoleucine.abbreviation(), "Ile");
        assert_eq!(AminoAcid::Leucine.abbreviation(), "Leu");
        assert_eq!(AminoAcid::Lysine.abbreviation(), "Lys");
        assert_eq!(AminoAcid::Methionine.abbreviation(), "Met");
        assert_eq!(AminoAcid::Phenylalanine.abbreviation(), "Phe");
        assert_eq!(AminoAcid::Proline.abbreviation(), "Pro");
        assert_eq!(AminoAcid::Serine.abbreviation(), "Ser");
        assert_eq!(AminoAcid::Threonine.abbreviation(), "Thr");
        assert_eq!(AminoAcid::Tryptophan.abbreviation(), "Trp");
        assert_eq!(AminoAcid::Tyrosine.abbreviation(), "Tyr");
        assert_eq!(AminoAcid::Valine.abbreviation(), "Val");
        assert_eq!(AminoAcid::Stop.abbreviation(), "Stop");
    }

    #[test]
    fn test_short_abbreviations() {
        assert_eq!(AminoAcid::Alanine.short_abbreviation(), 'A');
        assert_eq!(AminoAcid::Arginine.short_abbreviation(), 'R');
        assert_eq!(AminoAcid::Asparagine.short_abbreviation(), 'N');
        assert_eq!(AminoAcid::AsparticAcid.short_abbreviation(), 'D');
        assert_eq!(AminoAcid::Cysteine.short_abbreviation(), 'C');
        assert_eq!(AminoAcid::GlutamicAcid.short_abbreviation(), 'E');
        assert_eq!(AminoAcid::Glutamine.short_abbreviation(), 'Q');
        assert_eq!(AminoAcid::Glycine.short_abbreviation(), 'G');
        assert_eq!(AminoAcid::Histidine.short_abbreviation(), 'H');
        assert_eq!(AminoAcid::Isoleucine.short_abbreviation(), 'I');
        assert_eq!(AminoAcid::Leucine.short_abbreviation(), 'L');
        assert_eq!(AminoAcid::Lysine.short_abbreviation(), 'K');
        assert_eq!(AminoAcid::Methionine.short_abbreviation(), 'M');
        assert_eq!(AminoAcid::Phenylalanine.short_abbreviation(), 'F');
        assert_eq!(AminoAcid::Proline.short_abbreviation(), 'P');
        assert_eq!(AminoAcid::Serine.short_abbreviation(), 'S');
        assert_eq!(AminoAcid::Threonine.short_abbreviation(), 'T');
        assert_eq!(AminoAcid::Tryptophan.short_abbreviation(), 'W');
        assert_eq!(AminoAcid::Tyrosine.short_abbreviation(), 'Y');
        assert_eq!(AminoAcid::Valine.short_abbreviation(), 'V');
        assert_eq!(AminoAcid::Stop.short_abbreviation(), 'X');
    }

    #[test]
    fn test_protein_to_string() {
        let protein = Protein::new(vec![
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
            AminoAcid::GlutamicAcid,
            AminoAcid::Glutamine,
            AminoAcid::Stop,
        ]);
        assert_eq!(protein.to_string(), "AlaArgAsnAspCysGluGlnStop");
    }

    #[test]
    fn k_peptides_returns_correct_peptides_of_length_3() {
        let protein = Protein::new(vec![
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
        ]);
        let peptides = protein.k_peptides(3);
        assert_eq!(
            peptides,
            vec![
                vec![
                    AminoAcid::Alanine,
                    AminoAcid::Arginine,
                    AminoAcid::Asparagine
                ],
                vec![
                    AminoAcid::Arginine,
                    AminoAcid::Asparagine,
                    AminoAcid::AsparticAcid
                ],
                vec![
                    AminoAcid::Asparagine,
                    AminoAcid::AsparticAcid,
                    AminoAcid::Cysteine
                ],
            ]
        );
    }

    #[test]
    fn k_peptides_returns_empty_when_length_is_greater_than_sequence() {
        let protein = Protein::new(vec![AminoAcid::Alanine, AminoAcid::Arginine]);
        let peptides = protein.k_peptides(3);
        assert!(peptides.is_empty());
    }

    #[test]
    fn k_peptides_handles_single_amino_acid_sequence() {
        let protein = Protein::new(vec![AminoAcid::Alanine]);
        let peptides = protein.k_peptides(1);
        assert_eq!(peptides, vec![vec![AminoAcid::Alanine]]);
    }

    #[test]
    fn k_peptides_return_no_duplicates() {
        let protein = Protein::new(vec![
            AminoAcid::Alanine,
            AminoAcid::Alanine,
            AminoAcid::Alanine,
        ]);
        let peptides = protein.k_peptides(2);
        assert_eq!(peptides, vec![vec![AminoAcid::Alanine, AminoAcid::Alanine]]);
    }

    #[test]
    fn peptides_returns_correct_peptides_for_interval() {
        let protein = Protein::new(vec![
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
        ]);
        let peptides = protein.peptides(Interval { start: 2, end: 3 });
        assert_eq!(
            peptides,
            vec![
                vec![AminoAcid::Alanine, AminoAcid::Arginine],
                vec![AminoAcid::Arginine, AminoAcid::Asparagine],
                vec![AminoAcid::Asparagine, AminoAcid::AsparticAcid],
                vec![AminoAcid::AsparticAcid, AminoAcid::Cysteine],
                vec![
                    AminoAcid::Alanine,
                    AminoAcid::Arginine,
                    AminoAcid::Asparagine
                ],
                vec![
                    AminoAcid::Arginine,
                    AminoAcid::Asparagine,
                    AminoAcid::AsparticAcid
                ],
                vec![
                    AminoAcid::Asparagine,
                    AminoAcid::AsparticAcid,
                    AminoAcid::Cysteine
                ],
            ]
        );
    }
}
