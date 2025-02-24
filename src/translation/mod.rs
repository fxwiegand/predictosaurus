use crate::transcription;
use crate::translation::amino_acids::AminoAcid;
use anyhow::Result;
use itertools::Itertools;

pub(crate) mod amino_acids;

/// Translates a DNA sequence to a protein sequence
pub(crate) fn dna_to_amino_acids(dna: &[u8]) -> Result<Vec<AminoAcid>> {
    let rna = transcription::transcribe_dna_to_rna(dna)?;
    let amino_acids = rna
        .iter()
        .tuples::<(_, _, _)>()
        .map(|codon| AminoAcid::from_codon(&[*codon.0, *codon.1, *codon.2]))
        .take_while(|result| result.is_ok())
        .collect::<Result<Vec<AminoAcid>>>()?;
    Ok(amino_acids)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_dna_to_protein() -> Result<()> {
        assert_eq!(
            dna_to_amino_acids(b"ATGCGCGTA")?,
            vec![
                AminoAcid::Methionine,
                AminoAcid::Arginine,
                AminoAcid::Valine,
            ]
        );
        assert_eq!(
            dna_to_amino_acids(b"ATGCGCGTAAATGCGCGT")?,
            vec![
                AminoAcid::Methionine,
                AminoAcid::Arginine,
                AminoAcid::Valine,
                AminoAcid::Asparagine,
                AminoAcid::Alanine,
                AminoAcid::Arginine,
            ]
        );
        Ok(())
    }
}
