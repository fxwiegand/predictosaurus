use bio::io::fasta;
use std::collections::HashMap;
use std::path::Path;

pub(crate) fn read_reference<P: AsRef<Path> + std::fmt::Debug>(
    fasta_path: P,
) -> HashMap<String, Vec<u8>> {
    let reader = fasta::Reader::from_file(fasta_path).expect("Error reading FASTA file");
    let mut genome_map = HashMap::new();

    for result in reader.records() {
        let record = result.expect("Error reading record");
        let id = record.id().to_string();
        let seq = record.seq().to_vec();
        genome_map.insert(id, seq);
    }

    genome_map
}

pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&nucleotide| match nucleotide {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => unreachable!("Invalid nucleotide in sequence"),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_fasta_file(contents: &str) -> NamedTempFile {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "{}", contents).unwrap();
        temp_file
    }

    #[test]
    fn reads_single_record_correctly() {
        let fasta_contents = ">seq1\nACGT";
        let temp_file = create_temp_fasta_file(fasta_contents);
        let result = read_reference(temp_file.path());
        assert_eq!(result.len(), 1);
        assert_eq!(result["seq1"], b"ACGT".to_vec());
    }

    #[test]
    fn reads_multiple_records_correctly() {
        let fasta_contents = ">seq1\nACGT\n>seq2\nTGCA";
        let temp_file = create_temp_fasta_file(fasta_contents);
        let result = read_reference(temp_file.path());
        assert_eq!(result.len(), 2);
        assert_eq!(result["seq1"], b"ACGT".to_vec());
        assert_eq!(result["seq2"], b"TGCA".to_vec());
    }

    #[test]
    #[should_panic(expected = "Error reading FASTA file")]
    fn panics_on_nonexistent_file() {
        let fasta_path = "nonexistent.fasta";
        let _ = read_reference(fasta_path);
    }

    #[test]
    fn handles_file_with_no_sequences() {
        let fasta_contents = ">seq1\n>seq2";
        let temp_file = create_temp_fasta_file(fasta_contents);
        let result = read_reference(temp_file.path());
        assert_eq!(result.len(), 2);
        assert_eq!(result["seq1"], Vec::<u8>::new());
        assert_eq!(result["seq2"], Vec::<u8>::new());
    }

    #[test]
    fn reverse_complement_of_empty_sequence() {
        let seq = b"";
        let result = reverse_complement(seq);
        assert_eq!(result, b"");
    }

    #[test]
    fn reverse_complement_of_single_nucleotide() {
        let seq = b"A";
        let result = reverse_complement(seq);
        assert_eq!(result, b"T");
    }

    #[test]
    fn reverse_complement_of_standard_sequence() {
        let seq = b"ACGT";
        let result = reverse_complement(seq);
        assert_eq!(result, b"ACGT");
    }

    #[test]
    fn reverse_complement_of_palindromic_sequence() {
        let seq = b"AGCT";
        let result = reverse_complement(seq);
        assert_eq!(result, b"AGCT");
    }

    #[test]
    #[should_panic(expected = "Invalid nucleotide in sequence")]
    fn reverse_complement_with_invalid_nucleotide() {
        let seq = b"ACGTX";
        let _ = reverse_complement(seq);
    }
}
