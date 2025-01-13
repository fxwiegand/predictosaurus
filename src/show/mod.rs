use crate::graph::paths::{Weight, Cds};
use anyhow::Result;
use csv::Writer;
use itertools::Itertools;
use std::io::Write;
use std::path::{Path, PathBuf};
use tera::Tera;

/// Renders a Vega-Lite JSON file visualizing the given paths and writes it to the specified output directory with the feature as the file name.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the directory where the rendered JSON file will be saved.
/// * `paths` - A slice of `Weight` structs that will be serialized and passed to the template.
/// * `feature` - A string that represents the feature name, which will be used as the filename.
pub(crate) fn render_vl_paths(output_path: &PathBuf, paths: &[Weight], cds: Cds) -> Result<()> {
    let template = include_str!("../../resources/templates/paths.vl.json.tera");
    let mut context = tera::Context::new();
    context.insert("paths", paths);
    std::fs::write(
        Path::new(output_path).join(format!("{}.json", cds.name())),
        Tera::one_off(template, &context, false)?,
    )?;
    Ok(())
}

/// Renders a TSV file from the given paths and writes it to the specified output directory with the feature as the file name.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the directory where the rendered TSV file will be saved.
/// * `paths` - A slice of `Weight` structs that will be serialized and written to the TSV file.
/// * `feature` - A string that represents the feature name, which will be used as the filename.
pub(crate) fn render_tsv_paths(output_path: &PathBuf, paths: &[Weight], cds: Cds) -> Result<()> {
    let mut wtr = Writer::from_path(Path::new(output_path).join(format!("{}.tsv", cds.name())))?;
    for path in paths {
        wtr.serialize(path)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Renders an HTML file from the given paths and writes it to the specified output directory with the feature as the file name.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the directory where the rendered HTML file will be saved.
/// * `paths` - A slice of `Weight` structs that will be serialized and passed to the template.
/// * `feature` - A string that represents the feature name, which will be used as the filename.
pub(crate) fn render_html_paths(output_path: &PathBuf, paths: &[Weight], cds: Cds) -> Result<()> {
    let template = include_str!("../../resources/templates/paths.html.tera");
    let mut context = tera::Context::new();
    context.insert("paths", paths);
    std::fs::write(
        Path::new(output_path).join(format!("{}.html", cds.name())),
        Tera::one_off(template, &context, false)?,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impact::Impact;
    use serde_json::Value;

    #[test]
    fn render_vl_paths_creates_file_with_correct_content() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.into_path();
        let paths = vec![Weight {
            index: 1,
            path: Some(1),
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let cds = Cds::new("1".to_string(), "some feature".to_string(), 1, 100);
        render_vl_paths(&output_path, &paths, cds.clone()).unwrap();
        let file_path = Path::new(&output_path).join(format!("{}.json", cds.name()));
        assert!(file_path.exists());
        let file_content = std::fs::read_to_string(file_path).unwrap();
        assert!(serde_json::from_str::<Value>(&file_content).is_ok());
    }

    #[test]
    fn render_tsv_paths_creates_file_with_correct_content() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.into_path();
        let paths = vec![Weight {
            index: 1,
            path: Some(1),
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let cds = Cds::new("1".to_string(), "some feature".to_string(), 1, 100);
        render_tsv_paths(&output_path, &paths, cds.clone()).unwrap();
        let file_path = Path::new(&output_path).join(format!("{}.tsv", cds.name()));
        assert!(file_path.exists());
        let mut rdr = csv::Reader::from_path(file_path).unwrap();
        let result: Vec<Weight> = rdr.deserialize().map(|r| r.unwrap()).collect();
        assert_eq!(result, paths);
    }

    #[test]
    fn render_html_paths_creates_file_with_correct_content() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.into_path();
        let paths = vec![Weight {
            index: 1,
            path: Some(1),
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let cds = Cds::new("1".to_string(), "some feature".to_string(), 1, 100);
        render_html_paths(&output_path, &paths, cds.clone()).unwrap();
        let file_path = Path::new(&output_path).join(format!("{}.html", cds.name()));
        assert!(file_path.exists());
    }
}
