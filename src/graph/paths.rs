use crate::graph::{shift_phase, NodeType, VariantGraph};
use crate::impact::Impact;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use petgraph::graph::NodeIndex;
use std::cmp::max;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

pub(crate) struct Weight {
    pub(crate) index: usize,
    pub(crate) vaf: f32,
    pub(crate) impact: Impact,
    pub(crate) reason: String,
    pub(crate) consequence: String,
    pub(crate) sample: String,
}

impl HaplotypePath {
    pub(crate) fn impact(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> Result<Impact> {
        let mut impact = Impact::None;
        let ref_phase = phase;
        let mut phase = phase;
        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            let new_impact = node.impact(ref_phase, phase, reference, strand)?;
            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
            impact = max(impact, new_impact);
        }
        Ok(impact)
    }

    pub(crate) fn weights(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> Result<Vec<Weight>> {
        let mut weights = Vec::new();
        let ref_phase = phase;
        let mut phase = phase;
        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            let impact = node.impact(ref_phase, phase, reference, strand)?;
            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
            for (sample, vaf) in node.vaf.iter() {
                weights.push(Weight {
                    index: node_index.index(),
                    vaf: *vaf,
                    impact,
                    reason: node.reason(ref_phase, phase, reference, strand)?,
                    consequence: "TODO".to_string(),
                    sample: sample.clone(),
                });
            }
        }
        Ok(weights)
    }

    pub(crate) fn display(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> anyhow::Result<String> {
        let ref_phase = phase;
        let mut phase = phase;
        let mut protein = String::new();
        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            if let NodeType::Ref(_) = node.node_type {
                continue;
            }
            let ref_amino_acid = node.reference_amino_acid(ref_phase, reference, strand)?;
            let alt_amino_acid = node.variant_amino_acids(phase, reference, strand)?;
            protein.push_str(&format!(
                "{:?} -> {:?} ({:?})\n",
                ref_amino_acid,
                alt_amino_acid,
                node.impact(ref_phase, phase, reference, strand)?
            ));
            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
        }
        Ok(protein)
    }
}

mod tests {
    use crate::graph::node::{Node, NodeType};
    use crate::graph::paths::HaplotypePath;
    use crate::graph::{Edge, EventProbs, VariantGraph};
    use petgraph::{Directed, Graph};
    use std::collections::HashMap;
}
