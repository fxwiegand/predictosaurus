use anyhow::{anyhow, Result};

/// Transcribes a DNA sequence to an RNA sequence i.e. replaces all `T` with `U`
pub fn transcribe_dna_to_rna(dna: &[u8]) -> Result<Vec<u8>> {
    let mut rna = Vec::with_capacity(dna.len());
    for &nucleotide in dna {
        rna.push(match nucleotide {
            b'A' => b'A',
            b'T' => b'U',
            b'C' => b'C',
            b'G' => b'G',
            invalid => {
                return Err(anyhow!(
                    "Invalid nucleotide in DNA sequence: {}",
                    invalid as char
                ))
            }
        });
    }
    Ok(rna)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_transcribe_dna_to_rna() -> Result<()> {
        assert_eq!(transcribe_dna_to_rna(b"ATCG")?, b"AUCG");
        assert_eq!(transcribe_dna_to_rna(b"GATTACA")?, b"GAUUACA");
        assert_eq!(transcribe_dna_to_rna(b"CGTACG")?, b"CGUACG");
        assert_eq!(transcribe_dna_to_rna(b"TAA")?, b"UAA");
        Ok(())
    }

    #[test]
    fn test_transcribe_empty_dna() -> Result<()> {
        assert_eq!(transcribe_dna_to_rna(b"")?, Vec::<u8>::new());
        Ok(())
    }

    #[test]
    fn test_transcribe_dna_with_invalid_characters() -> Result<()> {
        match transcribe_dna_to_rna(b"\
        ") {
            Err(e) => assert_eq!(e.to_string(), "Invalid nucleotide in DNA sequence: B"),
            _ => panic!("Expected an error"),
        }
        Ok(())
    }
}
