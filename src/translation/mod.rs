use crate::transcription;
use amino_acids::Protein;
use anyhow::Result;
use itertools::Itertools;

mod amino_acids;

pub(crate) fn dna_to_protein(dna: &[u8]) -> Result<Protein> {
    let rna = transcription::transcribe_dna_to_rna(dna)?;
    let protein = rna
        .iter()
        .tuples::<(_, _, _)>()
        .map(|codon| amino_acids::AminoAcid::from_codon(&[*codon.0, *codon.1, *codon.2]))
        .collect::<Result<Vec<amino_acids::AminoAcid>>>()?;
    Ok(Protein::new(protein))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_dna_to_protein() -> Result<()> {
        assert_eq!(
            dna_to_protein(b"ATGCGCGTA")?.to_string(),
            "MetArgVal".to_string()
        );
        assert_eq!(
            dna_to_protein(b"ATGCGCGTAAATGCGCGT")?.to_string(),
            "MetArgValAsnAlaArg".to_string()
        );
        Ok(())
    }

    #[test]
    fn test_dna_to_protein_with_incomplete_codon() -> Result<()> {
        assert_eq!(
            dna_to_protein(b"ATGCGCGT")?.to_string(),
            "MetArg".to_string()
        );
        Ok(())
    }
}
