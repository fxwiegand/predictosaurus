use crate::graph::paths::{Cds, Weight};
use crate::graph::score::HaplotypeFrequency;
use anyhow::anyhow;
use anyhow::Result;
use csv::{Writer, WriterBuilder};
use itertools::Itertools;
use rust_htslib::bcf::header::SampleSubset;
use std::collections::HashMap;
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

/// Renders a TSV file from the given paths and writes it to the specified output directory with the feature as the file name.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the directory where the rendered TSV file will be saved.
/// * `paths` - A slice of `Weight` structs that will be serialized and written to the TSV file.
/// * `feature` - A string that represents the feature name, which will be used as the filename.
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

/// Renders an HTML file from the given paths and writes it to the specified output directory with the feature as the file name.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the directory where the rendered HTML file will be saved.
/// * `paths` - A slice of `Weight` structs that will be serialized and passed to the template.
/// * `feature` - A string that represents the feature name, which will be used as the filename.
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

/// Writes given per-haplotype scores into a single TSV file with the columns:
/// `transcript`, `score`, and one column per sample . Each row corresponds to one (score, per-sample) pair for a given transcript.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the output TSV file.
/// * `scores` - A reference to a `HashMap` that holds the scores, likelihoods, and haplotypes for each transcript.
pub(crate) fn render_scores(
    output_path: &PathBuf,
    scores: &HashMap<String, Vec<(f64, HaplotypeFrequency, String)>>,
) -> Result<()> {
    let mut wtr = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(Path::new(output_path))?;

    let samples: Vec<String> = scores
        .values()
        .flat_map(|scores| {
            scores
                .iter()
                .flat_map(|(_, sample_scores, _)| sample_scores.keys().cloned())
        })
        .unique()
        .collect();

    let mut headers = vec!["transcript", "score", "haplotype"];
    headers.extend(samples.iter().map(|s| s.as_str()));
    wtr.write_record(headers)?;

    for (transcript, hap_scores) in scores {
        for (score_val, sample_scores, haplotype) in hap_scores {
            let mut row = vec![
                transcript.to_string(),
                score_val.to_string(),
                haplotype.to_string(),
            ];
            for sample in &samples {
                let val = sample_scores.get(sample).ok_or_else(|| {
                    anyhow!(
                        "No score found for sample '{}' in transcript '{}'",
                        sample,
                        transcript
                    )
                })?;
                row.push(val.to_string());
            }
            wtr.write_record(row)?;
        }
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

    #[test]
    fn test_render_scores() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.keep().join("scores.tsv");
        let scores = HashMap::from([
            (
                "transcript1".to_string(),
                vec![(0.8, HashMap::new(), "c.[100A>G;105C>T]".to_string())],
            ),
            (
                "transcript2".to_string(),
                vec![(0.6, HashMap::new(), "c.[100A>G;105C>T]".to_string())],
            ),
        ]);
        render_scores(&output_path, &scores).unwrap();
        assert!(output_path.exists());
    }
}
