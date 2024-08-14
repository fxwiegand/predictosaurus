use anyhow::Result;
use rust_htslib::bcf::{Read, Reader};
use std::path::Path;

pub(crate) fn extract_event_names(bcf_file: &Path) -> Vec<String> {
    let mut event_names = Vec::new();

    // Open BCF file
    let reader = Reader::from_path(bcf_file).expect("Failed to open BCF file");

    // Iterate over header records
    for record in reader.header().header_records() {
        if let rust_htslib::bcf::HeaderRecord::Info { key: _, values, .. } = record {
            for info in values.values() {
                // Check if the INFO field starts with "PROB_"
                if let Some(stripped) = info.strip_prefix("PROB_") {
                    // Extract the event name from the INFO field key
                    let event_name = stripped.to_string();
                    event_names.push(event_name);
                }
            }
        }
    }

    event_names
}

// Get all occuring target chromosomes from a BCF file
pub(crate) fn get_targets(calls: &Path) -> Result<Vec<String>> {
    let mut targets = Vec::new();
    let reader = Reader::from_path(calls)?;

    for record in reader.header().header_records() {
        if let rust_htslib::bcf::HeaderRecord::Contig { key: _, values } = record {
            if let Some(id) = values.get("ID") {
                targets.push(id.to_string());
            }
        }
    }

    Ok(targets)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_extract_event_names_from_bcf() {
        let vcf_file_path = PathBuf::from("tests/resources/calls.bcf");
        let event_names = extract_event_names(&vcf_file_path);
        let expected_event_names = vec![
            "PRESENT".to_string(),
            "ARTIFACT".to_string(),
            "ABSENT".to_string(),
        ];
        assert_eq!(event_names, expected_event_names);
    }

    #[test]
    fn get_chromosomes_returns_correct_chromosomes() {
        let bcf_file_path = PathBuf::from("tests/resources/calls.bcf");
        let chromosomes = get_targets(&bcf_file_path).unwrap();
        let expected_chromosomes = vec!["OX512233.1".to_string()];
        assert_eq!(chromosomes, expected_chromosomes);
    }
}
