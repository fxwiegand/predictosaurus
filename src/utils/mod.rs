use anyhow::Result;
use std::fs;
use std::path::Path;

pub(crate) mod bcf;

pub(crate) fn create_output_dir(output_path: &Path) -> Result<()> {
    if !output_path.exists() {
        fs::create_dir_all(output_path)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_create_output_dir() {
        let output_path = PathBuf::from("/tmp/test_output_dir");
        create_output_dir(&output_path).unwrap();
        assert!(output_path.exists());
        fs::remove_dir(&output_path).unwrap();
    }
}
