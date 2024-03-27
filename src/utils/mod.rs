use std::path::PathBuf;
use serde::Deserialize;
use anyhow::Result;

pub(crate) mod bcf;

#[derive(Debug, Deserialize)]
pub(crate) struct AlignmentProperties {
    #[serde(rename = "max_read_len")]
    pub(crate) max_read_length: u32,
}

impl AlignmentProperties {
    pub(crate) fn from_file(file: &PathBuf) -> Result<Self> {
        let properties = std::fs::read_to_string(file)?;
        let properties: AlignmentProperties = serde_json::from_str(&properties)?;
        Ok(properties)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deserialize_alignment_properties() {
        let alignment_properties = r#"{"max_read_len": 150, "something": "else"}"#;
        let properties: AlignmentProperties = serde_json::from_str(alignment_properties).unwrap();
        assert_eq!(properties.max_read_length, 150);
    }

    #[test]
    fn test_alignment_properties_from_file() {
        let properties = AlignmentProperties::from_file(&PathBuf::from("tests/resources/alignment-properties.json")).unwrap();
        assert_eq!(properties.max_read_length, 100);
    }
}