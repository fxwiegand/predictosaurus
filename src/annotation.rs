use anyhow::Result;
use bio::bio_types::strand::Strand;
use genebears::{
    AnnotateOptions, AnnotatedVariant, ClientConfig, GeneBearError, GeneBears, Genome, Variant,
};
use itertools::Itertools;
use log::warn;
use serde::{Deserialize, Serialize};

use crate::graph::node::Node;
use crate::graph::transcript::Transcript;

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct Annotation {
    pub(crate) revel_score: Option<f64>,
    pub(crate) acmg_score: Option<f64>,
    pub(crate) spliceai_score: Option<f64>,
    pub(crate) alphamissense_score: Option<f64>,
}

impl Annotation {
    pub(crate) fn from_haplotype(
        haplotype: &[Node],
        transcript: &Transcript,
        genome_build: Genome,
    ) -> Result<Self> {
        let mut variants: Vec<_> = haplotype
            .iter()
            .filter(|n| n.node_type.is_variant())
            .collect();
        if transcript.strand == Strand::Reverse {
            variants.reverse();
        }
        let in_phase: Vec<&Node> = variants
            .into_iter()
            .scan(0i64, |running_fs, node| {
                let fs = node.frameshift();
                let in_phase = *running_fs % 3 == 0;
                *running_fs += fs;
                Some((node, in_phase))
            })
            .filter_map(|(node, in_phase)| in_phase.then_some(node))
            .collect();
        let variants = in_phase
            .iter()
            .map(|node| {
                Variant::new(
                    transcript.target.to_string(),
                    node.pos as u64,
                    node.reference_allele.to_string(),
                    node.alternative_allele.to_string(),
                )
            })
            .collect_vec();
        let client = GeneBears::new(ClientConfig::default())?;
        let rt = tokio::runtime::Runtime::new()?;

        let results = rt.block_on(client.annotate_variants(
            &variants,
            genome_build,
            AnnotateOptions::default(),
        ));

        let results = match results {
            Ok(r) => r,
            Err(GeneBearError::ApiServerError { message, .. }) => {
                warn!(
                    "GeneBe failed to annotate variants for {}: {}. Scores will be absent.",
                    transcript.target, message
                );
                vec![]
            }
            Err(e) => return Err(anyhow::Error::from(e)),
        };

        Ok(Self {
            revel_score: max_score(results.iter().map(|r| r.revel_score)),
            acmg_score: max_score(results.iter().map(|r| r.acmg_score)),
            spliceai_score: max_score(results.iter().map(|r| r.spliceai_max_score)),
            alphamissense_score: max_score(results.iter().map(|r| r.alphamissense_score)),
        })
    }
}

fn max_score(iter: impl Iterator<Item = Option<f64>>) -> Option<f64> {
    iter.flatten().reduce(f64::max)
}
