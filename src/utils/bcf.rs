use rust_htslib::bcf::{Read, Reader};
use std::path::Path;

pub(crate) fn extract_event_names(bcf_file: &Path) -> Vec<String> {
    let mut event_names = Vec::new();

    // Open BCF file
    let reader = Reader::from_path(&bcf_file).expect("Failed to open BCF file");

    // Iterate over header records
    for record in reader.header().header_records() {
        if let rust_htslib::bcf::HeaderRecord::Info { key: _, values, .. } = record {
            for info in values.values() {
                // Check if the INFO field starts with "PROB_"
                if info.starts_with("PROB_") {
                    // Extract the event name from the INFO field key
                    let event_name = info["PROB_".len()..].to_string();
                    event_names.push(event_name);
                }
            }
        }
    }

    event_names
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
}
