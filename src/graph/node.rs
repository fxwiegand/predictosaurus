use crate::graph::hgvs::hgvsc;
use crate::graph::transcript::Transcript;
use crate::graph::EventProbs;
use crate::transcription;
use crate::translation::amino_acids::AminoAcid;
use anyhow::anyhow;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use log::warn;
use rust_htslib::bcf::Record;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Display;
use std::str::FromStr;
use varlociraptor::variants::evidence::observations::read_observation::ProcessedReadObservation;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) enum NodeType {
    Variant,
    Reference,
}

impl FromStr for NodeType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Variant" => Ok(NodeType::Variant),
            "Reference" => Ok(NodeType::Reference),
            _ => Err(anyhow!("Invalid node type: {}", s)),
        }
    }
}

impl Display for NodeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NodeType::Variant => write!(f, "Variant"),
            NodeType::Reference => write!(f, "Reference"),
        }
    }
}

impl NodeType {
    pub(crate) fn is_variant(&self) -> bool {
        match self {
            NodeType::Variant => true,
            NodeType::Reference => false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) struct Node {
    pub(crate) node_type: NodeType,
    pub(crate) reference_allele: String,
    pub(crate) alternative_allele: String,
    pub(crate) vaf: HashMap<String, f32>,
    pub(crate) probs: EventProbs,
    pub(crate) pos: i64,
    pub(crate) index: u32,
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at position {}", self.node_type, self.pos)
    }
}

impl Node {
    pub(crate) fn from_records(
        calls_record: &Record,
        _observations_record: &HashMap<&&String, Vec<ProcessedReadObservation>>,
        event_probs: &EventProbs,
        node_type: NodeType,
        samples: &[String],
        index: u32,
    ) -> Self {
        let vafs = calls_record.format(b"AF").float().unwrap();
        let alleles = calls_record.alleles();
        let ref_allele = String::from_utf8(alleles[0].to_vec()).unwrap();
        let alt_allele = String::from_utf8(alleles[1].to_vec()).unwrap();
        Node {
            node_type,
            reference_allele: ref_allele,
            alternative_allele: alt_allele,
            vaf: samples
                .iter()
                .zip(vafs.iter())
                .map(|(s, v)| (s.to_string(), v[0]))
                .collect(),
            probs: event_probs.to_owned(),
            pos: calls_record.pos(),
            index,
        }
    }

    pub(crate) fn new(
        node_type: NodeType,
        pos: i64,
        reference_allele: String,
        alternative_allele: String,
    ) -> Self {
        Node {
            node_type,
            reference_allele,
            alternative_allele,
            vaf: Default::default(),
            probs: EventProbs(HashMap::new()),
            pos,
            index: 0,
        }
    }

    // Returns whether the node is a SNV
    pub(crate) fn is_snv(&self) -> bool {
        match &self.node_type {
            NodeType::Variant => self.alternative_allele.len() == 1,
            NodeType::Reference => false,
        }
    }

    /// Returns the frameshift caused by the variant at this node.
    pub(crate) fn frameshift(&self) -> i64 {
        match &self.node_type {
            NodeType::Variant => {
                self.alternative_allele.len() as i64 - self.reference_allele.len() as i64
            }
            NodeType::Reference => 0,
        }
    }

    /// Returns the 0-based position of the node on the reverse strand.
    pub(crate) fn position_on_reverse_strand(&self, reference_length: usize) -> i64 {
        (reference_length as i64 - 1) - self.pos
    }

    /// Returns the maximum VAF across all samples
    pub(crate) fn max_vaf(&self) -> f32 {
        *self
            .vaf
            .values()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&0.0)
    }

    pub(crate) fn hgvs_notation(&self, transcript: &Transcript) -> String {
        hgvsc(
            transcript,
            self.pos as u64,
            &self.reference_allele,
            &self.alternative_allele,
        )
        .unwrap_or_default()
    }

}

pub(crate) fn node_distance(node1: &u32, node2: &u32) -> u32 {
    node2 - node1
}

pub(crate) fn nodes_in_between(node1: &u32, node2: &u32, nodes: Vec<&u32>) -> usize {
    nodes
        .iter()
        .filter(|n| n < &&node2 && &&node1 < n)
        .filter(|n| node_distance(node1, n) != 0)
        .filter(|n| node_distance(n, node2) != 0)
        .count()
}

mod tests {
    use crate::graph::node::{node_distance, nodes_in_between, Node, NodeType};
    use crate::graph::paths::Cds;
    use crate::graph::transcript::Transcript;
    use crate::graph::{transcript, Edge, EventProbs};
    use crate::translation::amino_acids::AminoAcid;
    use bio::bio_types::strand::Strand;
    use itertools::Itertools;
    use petgraph::{Directed, Graph};
    use std::collections::HashMap;
    use std::str::FromStr;

    #[test]
    fn test_nodes_in_between() {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let vaf = HashMap::new();
        let ep = EventProbs(HashMap::new());
        let weight_1 = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 1,
            index: 1,
        };
        let node0 = graph.add_node(weight_1.clone());
        let node1 = graph.add_node(weight_1.clone());
        let weight_2 = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 4,
            index: 2,
        };
        let _node2 = graph.add_node(weight_2.clone());
        let node3 = graph.add_node(weight_2.clone());
        let weight_3 = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 5,
            index: 3,
        };
        let node4 = graph.add_node(weight_3.clone());
        let node5 = graph.add_node(weight_3.clone());

        let nodes = [node0, node1, node3, node4, node5]
            .iter()
            .map(|n| &graph.node_weight(*n).unwrap().index)
            .collect_vec();
        let nodes_in_between = nodes_in_between(
            &graph.node_weight(node0).unwrap().index,
            &graph.node_weight(node5).unwrap().index,
            nodes,
        );
        assert_eq!(nodes_in_between, 1);
    }

    #[test]
    fn test_node_distance() {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let vaf = HashMap::new();
        let ep = EventProbs(HashMap::new());
        let weight_1 = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 1,
            index: 1,
        };
        let node0 = graph.add_node(weight_1.clone());
        let node1 = graph.add_node(weight_1.clone());
        let weight_2 = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 4,
            index: 2,
        };
        let node2 = graph.add_node(weight_2.clone());
        let node3 = graph.add_node(weight_2.clone());
        let weight_3 = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 5,
            index: 3,
        };
        let node4 = graph.add_node(weight_3.clone());
        let _node5 = graph.add_node(weight_3.clone());

        let distance = node_distance(
            &graph.node_weight(node0).unwrap().index,
            &graph.node_weight(node2).unwrap().index,
        );
        assert_eq!(distance, 1);
        let distance_2 = node_distance(
            &graph.node_weight(node0).unwrap().index,
            &graph.node_weight(node1).unwrap().index,
        );
        assert_eq!(distance_2, 0);
        let distance_3 = node_distance(
            &graph.node_weight(node1).unwrap().index,
            &graph.node_weight(node4).unwrap().index,
        );
        assert_eq!(distance_3, 2);
        let distance_4 = node_distance(
            &graph.node_weight(node3).unwrap().index,
            &graph.node_weight(node4).unwrap().index,
        );
        assert_eq!(distance_4, 1);
        let distance_5 = node_distance(
            &graph.node_weight(node1).unwrap().index,
            &graph.node_weight(node3).unwrap().index,
        );
        assert_eq!(distance_5, 1);
    }

    #[test]
    fn test_frameshift() {
        let node = Node::new(NodeType::Variant, 2, "A".to_string(), "A".to_string());
        let frameshift = node.frameshift();
        assert_eq!(frameshift, 0);
        let node = Node::new(NodeType::Variant, 2, "A".to_string(), "AT".to_string());
        let frameshift = node.frameshift();
        assert_eq!(frameshift, 1);
        let node = Node::new(NodeType::Variant, 2, "A".to_string(), "".to_string());
        let frameshift = node.frameshift();
        assert_eq!(frameshift, -1);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_position() {
        let node = Node::new(NodeType::Reference, 5, "".to_string(), "".to_string());
        let reference_length = 10;
        assert_eq!(node.position_on_reverse_strand(reference_length), 4);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_middle_position() {
        let node = Node::new(NodeType::Reference, 4, "".to_string(), "".to_string());
        let reference_length = 9;
        assert_eq!(node.position_on_reverse_strand(reference_length), 4);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_start_position() {
        let node = Node::new(NodeType::Reference, 0, "".to_string(), "".to_string());
        let reference_length = 10;
        assert_eq!(node.position_on_reverse_strand(reference_length), 9);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_end_position() {
        let node = Node::new(NodeType::Reference, 10, "".to_string(), "".to_string());
        let reference_length = 11;
        assert_eq!(node.position_on_reverse_strand(reference_length), 0);
    }

    #[test]
    fn from_str_parses_variant_node_type() {
        let result = NodeType::from_str("Variant").unwrap();
        assert_eq!(result, NodeType::Variant);
    }

    #[test]
    fn from_str_parses_reference_node_type() {
        let result = NodeType::from_str("Reference").unwrap();
        assert_eq!(result, NodeType::Reference);
    }

    #[test]
    fn from_str_returns_error_for_invalid_node_type() {
        let result = NodeType::from_str("Invalid");
        assert!(result.is_err());
    }

    #[test]
    fn from_str_returns_error_for_empty_string() {
        let result = NodeType::from_str("");
        assert!(result.is_err());
    }

    #[test]
    fn test_node_display() {
        let var_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        let ref_node = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 99,
            index: 1,
        };

        assert_eq!(format!("{}", var_node), "Variant at position 42");
        assert_eq!(format!("{}", ref_node), "Reference at position 99");
    }

    #[test]
    fn test_node_max_vaf() {
        let mut vafs = HashMap::new();
        vafs.insert("A".to_string(), 0.05);
        vafs.insert("C".to_string(), 0.1);
        vafs.insert("G".to_string(), 0.02);
        vafs.insert("T".to_string(), 0.03);
        let var_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: vafs.clone(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        assert_eq!(var_node.max_vaf(), 0.1);
    }

    #[test]
    fn test_node_is_snv() {
        let snv_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        assert!(snv_node.is_snv());
        let indel_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "G".to_string(),
            alternative_allele: "CA".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        assert!(!indel_node.is_snv());
    }

    #[test]
    fn test_node_hgvs() {
        let var_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        let transcript = Transcript {
            feature: "ENST00000367770.8".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![Cds::new(0, 100, 2), Cds::new(200, 300, 1)],
        };
        assert_eq!(var_node.hgvs_notation(&transcript), "43C>A");
    }
}
