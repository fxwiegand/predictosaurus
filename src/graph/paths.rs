use crate::graph::{shift_phase, NodeType, VariantGraph};
use crate::impact::Impact;
use petgraph::graph::NodeIndex;
use std::cmp::max;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

impl HaplotypePath {
    pub(crate) fn impact(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
    ) -> anyhow::Result<Impact> {
        let mut impact = Impact::None;
        let ref_phase = phase;
        let mut phase = phase;
        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            let new_impact = node.impact(ref_phase, phase, reference)?;
            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
            impact = max(impact, new_impact);
        }
        Ok(impact)
    }

    pub(crate) fn weight(&self, graph: &VariantGraph) -> f32 {
        self.0
            .iter()
            .filter(|n| graph.graph.node_weight(**n).unwrap().node_type.is_variant())
            .map(|n| {
                let node = graph.graph.node_weight(*n).unwrap();
                let vaf_sum: f32 = node.vaf.values().sum();
                let weight = vaf_sum / node.vaf.values().len() as f32;
                weight
            })
            .product()
    }

    pub(crate) fn display(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
    ) -> anyhow::Result<String> {
        let ref_phase = phase;
        let mut phase = phase;
        let mut protein = String::new();
        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            if let NodeType::Ref(_) = node.node_type {
                continue;
            }
            let ref_amino_acid = node.reference_amino_acid(ref_phase, reference)?;
            let alt_amino_acid = node.variant_amino_acids(phase, reference)?;
            protein.push_str(&format!(
                "{} -> {:?} ({:?})\n",
                ref_amino_acid,
                alt_amino_acid,
                node.impact(ref_phase, phase, reference)?
            ));
            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
        }
        Ok(protein)
    }
}
