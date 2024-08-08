use crate::graph::{shift_phase, NodeType, VariantGraph};
use crate::impact::Impact;
use itertools::Itertools;
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

    pub(crate) fn weights(&self, graph: &VariantGraph, sample: String) -> Vec<f32> {
        self.0
            .iter()
            .filter(|n| graph.graph.node_weight(**n).unwrap().node_type.is_variant())
            .map(|n| {
                let node = graph.graph.node_weight(*n).unwrap();
                println!("{:?}", node.vaf);
                node.vaf.get(&sample).cloned().unwrap()
            })
            .collect_vec()
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

mod tests {
    use crate::graph::node::{Node, NodeType};
    use crate::graph::paths::HaplotypePath;
    use crate::graph::{Edge, EventProbs, VariantGraph};
    use petgraph::{Directed, Graph};
    use std::collections::HashMap;

    #[test]
    fn weights_returns_correct_values_for_multiple_variant_nodes() {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let mut vaf1 = HashMap::new();
        vaf1.insert("sample".to_string(), 0.5);
        let node1 = graph.add_node(Node {
            node_type: NodeType::Var("A".to_string()),
            vaf: vaf1,
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 0,
        });

        let mut vaf2 = HashMap::new();
        vaf2.insert("sample".to_string(), 0.7);
        let node2 = graph.add_node(Node {
            node_type: NodeType::Var("A".to_string()),
            vaf: vaf2,
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 0,
        });

        let path = HaplotypePath(vec![node1, node2]);
        let variant_graph = VariantGraph {
            graph,
            start: 0,
            end: 2,
            target: "test".to_string(),
        };
        let weights = path.weights(&variant_graph, "sample".to_string());
        assert_eq!(weights, vec![0.5, 0.7]);
    }
}
