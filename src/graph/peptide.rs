use crate::graph::EventProbs;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::dna_to_amino_acids;
use anyhow::Result;

pub(crate) struct Peptide {
    pub(crate) sequence: Vec<AminoAcid>,
    pub(crate) prob: EventProbs,
}

impl Peptide {
    pub(crate) fn from_rna(rna: Vec<u8>, prob: EventProbs) -> Result<Self> {
        let sequence = dna_to_amino_acids(&rna)?;
        Ok(Peptide { sequence, prob })
    }
}
