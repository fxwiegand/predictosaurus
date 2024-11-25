use crate::graph::{shift_phase, NodeType, VariantGraph};
use crate::impact::Impact;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use petgraph::graph::NodeIndex;
use serde::Serialize;
use std::cmp::max;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

#[derive(Debug, Serialize, Clone)]
pub(crate) struct Weight {
    pub(crate) index: usize,
    pub(crate) path: Option<usize>,
    pub(crate) vaf: f32,
    pub(crate) impact: Impact,
    pub(crate) reason: Option<String>,
    pub(crate) consequence: Option<String>,
    pub(crate) sample: String,
}

/// Returns the maximum impact of the individual nodes on the path
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

    /// Returns the weights of the individual nodes on the path as a vector of Weights
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
                    path: None,
                    vaf: *vaf,
                    impact,
                    reason: node.reason(ref_phase, phase, reference, strand)?,
                    consequence: None, // TODO: Implement consequence
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
    use crate::graph::{Edge, EventProbs, VariantGraph};
    use crate::impact::Impact;
    use bio::bio_types::strand::Strand;
    use petgraph::{Directed, Graph};
    use std::collections::HashMap;

    fn setup_variant_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let ref_node_vaf =
            HashMap::from([("sample1".to_string(), 0.5), ("sample2".to_string(), 0.5)]);
        let alt_node_vaf =
            HashMap::from([("sample1".to_string(), 0.3), ("sample2".to_string(), 0.7)]);
        let ref_node = graph.add_node(Node {
            node_type: NodeType::Ref("".to_string()),
            vaf: ref_node_vaf,
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 0,
        });
        let alt_node = graph.add_node(Node {
            node_type: NodeType::Var("ACGTTTGTTA".to_string()),
            vaf: alt_node_vaf,
            probs: EventProbs(HashMap::new()),
            pos: 6,
            index: 1,
        });
        graph.add_edge(
            ref_node,
            alt_node,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        VariantGraph {
            graph,
            start: 0,
            end: 0,
            target: "test".to_string(),
        }
    }

    #[test]
    fn weights_returns_correct_weights() {
        let graph = setup_variant_graph();
        let path = &graph.paths()[0];
        let reference = b"ACGTTTGTTA";
        let strand = Strand::Forward;
        let weights = path.weights(&graph, 0, reference, strand).unwrap();
        assert_eq!(weights.len(), 4);
        assert_eq!(weights[0].vaf, 0.5);
        assert_eq!(weights[2].impact, Impact::Moderate);
    }
}
