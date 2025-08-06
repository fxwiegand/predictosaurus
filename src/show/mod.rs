use crate::graph::paths::{Cds, Weight};
use anyhow::Result;
use csv::Writer;
use itertools::Itertools;
use std::io::Write;
use std::path::{Path, PathBuf};
use tera::Tera;

/// Generates a Vega-Lite JSON visualization of genomic paths and writes it to the specified output directory using the transcript name as the filename.
///
/// The function loads a Vega-Lite template, injects the provided `Weight` data, renders the template, and saves the resulting JSON to `{transcript}.json` in the given directory.
///
/// # Examples
///
/// ```
/// let output_dir = std::path::PathBuf::from("out");
/// let paths = vec![Weight { /* fields */ }];
/// render_vl_paths(&output_dir, &paths, "ENST00000367770".to_string()).unwrap();
/// ```
pub(crate) fn render_vl_paths(
    output_path: &PathBuf,
    paths: &[Weight],
    transcript: String,
) -> Result<()> {
    let template = include_str!("../../resources/templates/paths.vl.json.tera");
    let mut context = tera::Context::new();
    context.insert("paths", paths);
    std::fs::write(
        Path::new(output_path).join(format!("{transcript}.json")),
        Tera::one_off(template, &context, false)?,
    )?;
    Ok(())
}

/// Writes a slice of `Weight` structs to a TSV file named after the transcript in the specified output directory.
///
/// Each `Weight` struct is serialized as a row in the TSV file. The file is saved as `{transcript}.tsv` in the given directory.
///
/// # Examples
///
/// ```
/// # use your_crate::{render_tsv_paths, Weight};
/// # use std::path::PathBuf;
/// let output_dir = PathBuf::from("/tmp");
/// let weights = vec![Weight { /* fields */ }];
/// render_tsv_paths(&output_dir, &weights, "ENST00000367770".to_string()).unwrap();
/// // This creates "/tmp/ENST00000367770.tsv"
/// ```
pub(crate) fn render_tsv_paths(
    output_path: &PathBuf,
    paths: &[Weight],
    transcript: String,
) -> Result<()> {
    let mut wtr = Writer::from_path(Path::new(output_path).join(format!("{transcript}.tsv")))?;
    for path in paths {
        wtr.serialize(path)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Renders an HTML file visualizing the provided genomic paths and writes it to the specified output directory using the transcript name as the filename.
///
/// The function loads an HTML template, injects the given `Weight` data, and saves the rendered HTML as `{transcript}.html` in the output directory.
///
/// # Examples
///
/// ```
/// let output_dir = std::path::PathBuf::from("out");
/// let paths = vec![Weight { /* fields */ }];
/// render_html_paths(&output_dir, &paths, "ENST00000367770".to_string()).unwrap();
/// ```
pub(crate) fn render_html_paths(
    output_path: &PathBuf,
    paths: &[Weight],
    transcript: String,
) -> Result<()> {
    let template = include_str!("../../resources/templates/paths.html.tera");
    let mut context = tera::Context::new();
    context.insert("paths", paths);
    std::fs::write(
        Path::new(output_path).join(format!("{transcript}.html")),
        Tera::one_off(template, &context, false)?,
    )?;
    Ok(())
}

/// Writes a list of floating-point scores to a TSV file named after the transcript in the specified output directory.
///
/// Each score is written as a single-field record in the output file `{transcript}.tsv`.
///
/// # Arguments
///
/// * `output_path` - The directory where the TSV file will be created.
/// * `scores` - Slice of floating-point scores to write.
/// * `transcript` - The base name for the output TSV file.
///
/// # Errors
///
/// Returns an error if the file cannot be created, written to, or flushed.
///
/// # Examples
///
/// ```
/// let output_dir = std::env::temp_dir();
/// let scores = vec![0.95, 0.87, 0.76];
/// render_scores(&output_dir, &scores, "sample_transcript".to_string()).unwrap();
/// let tsv_path = output_dir.join("sample_transcript.tsv");
/// assert!(tsv_path.exists());
/// ```
pub(crate) fn render_scores(
    output_path: &PathBuf,
    scores: &[f64],
    transcript: String,
) -> Result<()> {
    let mut wtr = Writer::from_path(Path::new(output_path).join(format!("{transcript}.tsv")))?;
    for score in scores {
        wtr.write_record(&[score.to_string()])?;
    }
    wtr.flush()?;
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
        let output_path = temp_dir.keep();
        let paths = vec![Weight {
            index: 1,
            path: Some(1),
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let transcript = "some feature".to_string();
        render_vl_paths(&output_path, &paths, transcript.clone()).unwrap();
        let file_path = Path::new(&output_path).join(format!("{}.json", transcript));
        assert!(file_path.exists());
        let file_content = std::fs::read_to_string(file_path).unwrap();
        assert!(serde_json::from_str::<Value>(&file_content).is_ok());
    }

    #[test]
    fn render_tsv_paths_creates_file_with_correct_content() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.keep();
        let paths = vec![Weight {
            index: 1,
            path: Some(1),
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let transcript = "some feature".to_string();
        render_tsv_paths(&output_path, &paths, transcript.clone()).unwrap();
        let file_path = Path::new(&output_path).join(format!("{}.tsv", transcript));
        assert!(file_path.exists());
        let mut rdr = csv::Reader::from_path(file_path).unwrap();
        let result: Vec<Weight> = rdr.deserialize().map(|r| r.unwrap()).collect();
        assert_eq!(result, paths);
    }

    #[test]
    fn render_html_paths_creates_file_with_correct_content() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.keep();
        let paths = vec![Weight {
            index: 1,
            path: Some(1),
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let transcript = "some feature".to_string();
        render_html_paths(&output_path, &paths, transcript.clone()).unwrap();
        let file_path = Path::new(&output_path).join(format!("{}.html", transcript));
        assert!(file_path.exists());
    }
}
