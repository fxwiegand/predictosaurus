use std::io::Write;
use std::path::Path;
use anyhow::Result;
use itertools::Itertools;
use crate::graph::paths::Weight;
use tera::Tera;

/// Renders a Vega-Lite JSON file visualizing the given paths and writes it to the specified output directory with the feature as the file name.
///
/// # Arguments
///
/// * `output_path` - A string slice that holds the path of the directory where the rendered JSON file will be saved.
/// * `paths` - A vector of `Weight` structs that will be serialized and passed to the template.
/// * `feature` - A string that represents the feature name, which will be used as the filename.
pub(crate) fn render_vl_paths(output_path: &str, paths: Vec<Weight>, feature: String) -> Result<()> {
    let template = include_str!("../../resources/templates/paths.vl.json.tera");
    let mut context = tera::Context::new();
    context.insert("paths", &paths);
    std::fs::write(Path::new(output_path).join(format!("{}.json", feature)), Tera::one_off(template, &context, false)?)?;
    Ok(())
}


#[cfg(test)]
mod tests {
    use serde_json::Value;
    use crate::impact::Impact;
    use super::*;

    #[test]
    fn render_vl_paths_creates_file_with_correct_content() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().to_str().unwrap();
        let paths = vec![Weight {
            index: 1,
            vaf: 0.5,
            impact: Impact::High,
            reason: Some("Ile -> Met".to_string()),
            consequence: Some("loss".to_string()),
            sample: "sample1".to_string(),
        }];
        let feature = "test_feature".to_string();
        render_vl_paths(output_path, paths, feature.clone()).unwrap();
        let file_path = Path::new(output_path).join(format!("{}.json", feature));
        assert!(file_path.exists());
        let file_content = std::fs::read_to_string(file_path).unwrap();
        assert!(serde_json::from_str::<Value>(&file_content).is_ok());
    }
}
