use crate::graph::EventProbs;
use crate::impact::Impact;
use crate::transcription;
use crate::translation::amino_acids::AminoAcid;
use anyhow::anyhow;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use log::warn;
use rust_htslib::bcf::Record;
use rust_htslib::htslib::WAIT_ANY;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Display;
use std::str::FromStr;
use varlociraptor::variants::evidence::observations::read_observation::ProcessedReadObservation;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) enum NodeType {
    Var(String),
    Ref(String),
}

impl Display for NodeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NodeType::Var(alt) => write!(f, "Var({alt})"),
            NodeType::Ref(_) => write!(f, "Ref"),
        }
    }
}

impl FromStr for NodeType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.starts_with("Var(") && s.ends_with(')') {
            let inner = &s[4..s.len() - 1];
            Ok(NodeType::Var(inner.to_string()))
        } else if s == "Ref" {
            Ok(NodeType::Ref("".to_string()))
        } else {
            Err(anyhow!("Invalid node type"))
        }
    }
}

impl NodeType {
    pub(crate) fn is_variant(&self) -> bool {
        match self {
            NodeType::Var(_) => true,
            NodeType::Ref(_) => false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) struct Node {
    pub(crate) node_type: NodeType,
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
        Node {
            node_type,
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

    pub(crate) fn new(node_type: NodeType, pos: i64) -> Self {
        Node {
            node_type,
            vaf: Default::default(),
            probs: EventProbs(HashMap::new()),
            pos,
            index: 0,
        }
    }

    // Returns whether the node is a SNV
    pub(crate) fn is_snv(&self) -> bool {
        match &self.node_type {
            NodeType::Var(alt_allele) => alt_allele.len() == 1,
            NodeType::Ref(_) => false,
        }
    }

    /// Returns the frameshift caused by the variant at this node.
    pub(crate) fn frameshift(&self) -> i64 {
        match &self.node_type {
            NodeType::Var(alt_allele) => alt_allele.len() as i64 - 1,
            NodeType::Ref(_) => 0,
        }
    }

    /// Returns the 0-based position of the node on the reverse strand.
    pub(crate) fn position_on_reverse_strand(&self, reference_length: usize) -> i64 {
        (reference_length as i64 - 1) - self.pos
    }

    pub(crate) fn reference_amino_acid(
        &self,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> anyhow::Result<Option<AminoAcid>> {
        let p = match strand {
            Strand::Forward => self.pos,
            Strand::Reverse => self.position_on_reverse_strand(reference.len()),
            Strand::Unknown => return Err(anyhow!("Strand is unknown")),
        };
        let start_pos = p as usize - ((p - phase as i64) % 3) as usize;
        if start_pos + 3 > reference.len() {
            return Ok(None);
        }
        let ref_codon_bases = reference[start_pos..start_pos + 3].to_vec();
        if ref_codon_bases.contains(&b'N') {
            warn!("Found nucleotide 'N' in reference sequence for Node: {self}");
            return Ok(None);
        }
        let rna_codon = transcription::transcribe_dna_to_rna(&ref_codon_bases)?;
        Ok(AminoAcid::from_codon(&rna_codon).ok())
    }

    pub(crate) fn variant_amino_acids(
        &self,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> anyhow::Result<Vec<AminoAcid>> {
        match &self.node_type {
            NodeType::Var(alt_allele) => {
                let p = match strand {
                    Strand::Forward => self.pos,
                    Strand::Reverse => self.position_on_reverse_strand(reference.len()),
                    Strand::Unknown => return Err(anyhow!("Strand is unknown")),
                };
                let alt_allele = if strand == Strand::Reverse {
                    String::from_utf8(crate::utils::fasta::reverse_complement(
                        alt_allele.as_bytes(),
                    ))?
                } else {
                    alt_allele.to_string()
                };

                let start_pos = (p - ((p - phase as i64) % 3)) as usize;
                let position_in_codon = ((p - phase as i64) % 3) as usize;
                let needed_bases = if alt_allele.is_empty() { 4 } else { 3 };
                if start_pos + needed_bases > reference.len() || p - (phase as i64) < 0 {
                    // Variant is at the end of the reference sequence and there are not enough bases to form a codon
                    return Ok(vec![]);
                }
                let ref_codon_bases = reference[start_pos..start_pos + needed_bases].to_vec();
                let alt_codon_bases = [
                    &ref_codon_bases[..position_in_codon],
                    alt_allele.as_bytes(),
                    &ref_codon_bases[position_in_codon + 1..],
                ]
                .concat();
                if alt_codon_bases.contains(&b'N') {
                    return Ok(vec![]); // Don not warn since we probably warn already for reference_amino_acid
                }
                Ok(transcription::transcribe_dna_to_rna(&alt_codon_bases)?
                    .iter()
                    .tuples::<(_, _, _)>()
                    .map(|codon| AminoAcid::from_codon(&[*codon.0, *codon.1, *codon.2]).unwrap())
                    .collect_vec())
            }
            _ => {
                unreachable!("Reference node type has no variant amino acid")
            }
        }
    }

    pub(crate) fn reason(
        &self,
        ref_phase: u8,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> anyhow::Result<Option<String>> {
        if let NodeType::Ref(_) = self.node_type {
            return Ok(None);
        }
        let ref_amino_acid = self
            .reference_amino_acid(ref_phase, reference, strand)?
            .map_or("None".to_string(), |a| a.to_string());
        let alt_amino_acids = self
            .variant_amino_acids(phase, reference, strand)?
            .iter()
            .join(", ");
        Ok(Some(format!("{ref_amino_acid} -> {alt_amino_acids}")))
    }

    pub(crate) fn impact(
        &self,
        ref_phase: u8,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> anyhow::Result<Impact> {
        match &self.node_type {
            NodeType::Var(_) => {
                let ref_amino_acid = self.reference_amino_acid(ref_phase, reference, strand)?;
                match ref_amino_acid {
                    None => Ok(Impact::None),
                    Some(ref_amino_acid) => {
                        let alt_amino_acids = self.variant_amino_acids(phase, reference, strand)?;
                        match alt_amino_acids.len() {
                            0 => Ok(Impact::Moderate),
                            1 => {
                                let alt_amino_acid = alt_amino_acids.first().unwrap();
                                let impact = match (
                                    ref_amino_acid == *alt_amino_acid,
                                    ref_amino_acid,
                                    alt_amino_acid,
                                ) {
                                    (true, _, _) => Ok(Impact::Low),
                                    (false, _, AminoAcid::Stop) => Ok(Impact::High),
                                    (false, AminoAcid::Stop, _) => Ok(Impact::High),
                                    (false, AminoAcid::Methionine, _) => Ok(Impact::High), // TODO: Check if this is always automatic start lost or can Met occur anywhere in the protein?
                                    (false, _, _) => Ok(Impact::Moderate),
                                };
                                for amino_acid in alt_amino_acids {
                                    if amino_acid == AminoAcid::Stop {
                                        return Ok(Impact::High);
                                    }
                                }
                                impact
                            }
                            _ => Ok(Impact::Moderate),
                        }
                    }
                }
            }
            NodeType::Ref(_) => Ok(Impact::None),
        }
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
    use crate::graph::{Edge, EventProbs};
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
            node_type: NodeType::Ref("".to_string()),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 1,
            index: 1,
        };
        let node0 = graph.add_node(weight_1.clone());
        let node1 = graph.add_node(weight_1.clone());
        let weight_2 = Node {
            node_type: NodeType::Ref("".to_string()),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 4,
            index: 2,
        };
        let _node2 = graph.add_node(weight_2.clone());
        let node3 = graph.add_node(weight_2.clone());
        let weight_3 = Node {
            node_type: NodeType::Ref("".to_string()),
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
            node_type: NodeType::Ref("".to_string()),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 1,
            index: 1,
        };
        let node0 = graph.add_node(weight_1.clone());
        let node1 = graph.add_node(weight_1.clone());
        let weight_2 = Node {
            node_type: NodeType::Ref("".to_string()),
            vaf: vaf.clone(),
            probs: ep.clone(),
            pos: 4,
            index: 2,
        };
        let node2 = graph.add_node(weight_2.clone());
        let node3 = graph.add_node(weight_2.clone());
        let weight_3 = Node {
            node_type: NodeType::Ref("".to_string()),
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
    fn test_variant_amino_acid_with_different_phases() {
        let node = Node::new(NodeType::Var("A".to_string()), 2);
        let reference = b"ATGCGCGTA";
        let ile = node
            .variant_amino_acids(0, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(ile, AminoAcid::Isoleucine);
        let tyr = node
            .variant_amino_acids(1, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(tyr, AminoAcid::Tyrosine);
        let thr = node
            .variant_amino_acids(2, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(thr, AminoAcid::Threonine);
    }

    #[test]
    fn test_variant_amino_acid_with_deletion() {
        let node = Node::new(NodeType::Var("".to_string()), 2);
        let reference = b"ATGCCGT";
        let ile = node
            .variant_amino_acids(0, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(ile, AminoAcid::Isoleucine);
        let ser = node
            .variant_amino_acids(1, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(ser, AminoAcid::Serine);
        let pro = node
            .variant_amino_acids(2, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(pro, AminoAcid::Proline);
    }

    #[test]
    fn test_variant_amino_acid_with_one_inserted_base() {
        let node = Node::new(NodeType::Var("GG".to_string()), 1);
        let reference = b"AGCTCT";
        let arg = node
            .variant_amino_acids(0, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(arg, AminoAcid::Arginine);
        let gly = node
            .variant_amino_acids(1, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(gly, AminoAcid::Glycine);
        let node_2 = Node::new(NodeType::Var("GG".to_string()), 3);
        let arg = node_2
            .variant_amino_acids(2, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(arg, AminoAcid::Arginine);
    }

    #[test]
    fn test_variant_amino_acid_with_two_inserted_bases() {
        let node = Node::new(NodeType::Var("TCC".to_string()), 4);
        let reference = b"ATCATCATC";
        let ile = node
            .variant_amino_acids(0, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(ile, AminoAcid::Isoleucine);
        let ser = node
            .variant_amino_acids(1, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(ser, AminoAcid::Serine);
        let his = node
            .variant_amino_acids(2, reference, Strand::Forward)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(his, AminoAcid::Histidine);
    }

    #[test]
    fn test_variant_amino_acid_with_six_inserted_bases() {
        let node = Node::new(NodeType::Var("TCCAGTT".to_string()), 4);
        let reference = b"ATCATCATC";
        let ile_gln_phe = node
            .variant_amino_acids(0, reference, Strand::Forward)
            .unwrap();
        assert_eq!(
            ile_gln_phe,
            vec![
                AminoAcid::Isoleucine,
                AminoAcid::Glutamine,
                AminoAcid::Phenylalanine
            ]
        );
        let ser_ser_ser = node
            .variant_amino_acids(1, reference, Strand::Forward)
            .unwrap();
        assert_eq!(
            ser_ser_ser,
            vec![AminoAcid::Serine, AminoAcid::Serine, AminoAcid::Serine]
        );
        let his_pro_val = node
            .variant_amino_acids(2, reference, Strand::Forward)
            .unwrap();
        assert_eq!(
            his_pro_val,
            vec![AminoAcid::Histidine, AminoAcid::Proline, AminoAcid::Valine]
        );
    }

    #[test]
    fn test_variant_amino_acid_at_feature_end() {
        let node = Node::new(NodeType::Var("A".to_string()), 9);
        let reference = b"ATGCGCGTAT";
        let no_amino_acid = node
            .variant_amino_acids(0, reference, Strand::Forward)
            .unwrap();
        assert!(no_amino_acid.is_empty());
        let node2 = Node::new(NodeType::Var("A".to_string()), 0);
        let no_amino_acid_2 = node2
            .variant_amino_acids(0, reference, Strand::Reverse)
            .unwrap();
        assert!(no_amino_acid_2.is_empty());
    }

    #[test]
    fn test_variant_amino_acid_at_feature_start() {
        let node = Node::new(NodeType::Var("A".to_string()), 0);
        let reference = b"ATGCGCGTAT";
        let no_amino_acid = node
            .variant_amino_acids(1, reference, Strand::Forward)
            .unwrap();
        assert!(no_amino_acid.is_empty());
    }

    #[test]
    fn test_variant_amino_acid_on_backward_strand() {
        let node = Node::new(NodeType::Var("A".to_string()), 2);
        let forward_reference = b"ATGCGCGTA";
        let backward_reference = &{ crate::utils::fasta::reverse_complement(forward_reference) };
        assert_eq!(backward_reference, b"TACGCGCAT");
        let tyr = node
            .variant_amino_acids(0, backward_reference, Strand::Reverse)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(tyr, AminoAcid::Tyrosine);
        let arg = node
            .variant_amino_acids(1, backward_reference, Strand::Reverse)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(arg, AminoAcid::Arginine);
        let val = node
            .variant_amino_acids(2, backward_reference, Strand::Reverse)
            .unwrap()
            .first()
            .unwrap()
            .to_owned();
        assert_eq!(val, AminoAcid::Valine);
    }

    #[test]
    fn test_reference_amino_acid_with_different_phases() {
        let node = Node::new(NodeType::Ref("".to_string()), 2);
        let reference = b"ATGCGCGTA";
        let met = node
            .reference_amino_acid(0, reference, Strand::Forward)
            .unwrap()
            .unwrap();
        assert_eq!(met, AminoAcid::Methionine);
        let cys = node
            .reference_amino_acid(1, reference, Strand::Forward)
            .unwrap()
            .unwrap();
        assert_eq!(cys, AminoAcid::Cysteine);
        let ala = node
            .reference_amino_acid(2, reference, Strand::Forward)
            .unwrap()
            .unwrap();
        assert_eq!(ala, AminoAcid::Alanine);
    }

    #[test]
    fn test_reference_amino_acid_on_backward_strand() {
        let node = Node::new(NodeType::Ref("".to_string()), 2);
        let forward_reference = b"ATGCGCGTA";
        let backward_reference = &{ crate::utils::fasta::reverse_complement(forward_reference) };
        assert_eq!(backward_reference, b"TACGCGCAT");
        let his = node
            .reference_amino_acid(0, backward_reference, Strand::Reverse)
            .unwrap()
            .unwrap();
        assert_eq!(his, AminoAcid::Histidine);
        let arg = node
            .reference_amino_acid(1, backward_reference, Strand::Reverse)
            .unwrap()
            .unwrap();
        assert_eq!(arg, AminoAcid::Arginine);
        let ala = node
            .reference_amino_acid(2, backward_reference, Strand::Reverse)
            .unwrap()
            .unwrap();
        assert_eq!(ala, AminoAcid::Alanine);
    }

    #[test]
    fn test_empty_reference_amino_acid() {
        let node = Node::new(NodeType::Ref("".to_string()), 9);
        let reference = b"ATGCGCGTAT";
        let no_amino_acid = node
            .reference_amino_acid(0, reference, Strand::Forward)
            .unwrap();
        assert!(no_amino_acid.is_none());
        let node2 = Node::new(NodeType::Ref("".to_string()), 0);
        let no_amino_acid_2 = node2
            .reference_amino_acid(0, reference, Strand::Reverse)
            .unwrap();
        assert!(no_amino_acid_2.is_none());
    }

    #[test]
    fn test_frameshift() {
        let node = Node::new(NodeType::Var("A".to_string()), 2);
        let frameshift = node.frameshift();
        assert_eq!(frameshift, 0);
        let node = Node::new(NodeType::Var("AT".to_string()), 2);
        let frameshift = node.frameshift();
        assert_eq!(frameshift, 1);
        let node = Node::new(NodeType::Var("".to_string()), 2);
        let frameshift = node.frameshift();
        assert_eq!(frameshift, -1);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_position() {
        let node = Node::new(NodeType::Ref("".to_string()), 5);
        let reference_length = 10;
        assert_eq!(node.position_on_reverse_strand(reference_length), 4);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_middle_position() {
        let node = Node::new(NodeType::Ref("".to_string()), 4);
        let reference_length = 9;
        assert_eq!(node.position_on_reverse_strand(reference_length), 4);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_start_position() {
        let node = Node::new(NodeType::Ref("".to_string()), 0);
        let reference_length = 10;
        assert_eq!(node.position_on_reverse_strand(reference_length), 9);
    }

    #[test]
    fn position_on_reverse_strand_calculates_correctly_for_end_position() {
        let node = Node::new(NodeType::Ref("".to_string()), 10);
        let reference_length = 11;
        assert_eq!(node.position_on_reverse_strand(reference_length), 0);
    }

    #[test]
    fn from_str_parses_variant_node_type() {
        let result = NodeType::from_str("Var(A)").unwrap();
        assert_eq!(result, NodeType::Var("A".to_string()));
    }

    #[test]
    fn from_str_parses_reference_node_type() {
        let result = NodeType::from_str("Ref").unwrap();
        assert_eq!(result, NodeType::Ref("".to_string()));
    }

    #[test]
    fn from_str_returns_error_for_invalid_node_type() {
        let result = NodeType::from_str("Invalid");
        assert!(result.is_err());
    }

    #[test]
    fn from_str_returns_error_for_malformed_variant_node_type() {
        let result = NodeType::from_str("Var(A");
        assert!(result.is_err());
    }

    #[test]
    fn from_str_returns_error_for_empty_string() {
        let result = NodeType::from_str("");
        assert!(result.is_err());
    }

    #[test]
    fn reason_with_valid_reference_and_variant_amino_acids() {
        let node = Node::new(NodeType::Var("A".to_string()), 2);
        let reference = b"ATGCGCGTA";
        let result = node
            .reason(0, 0, reference, Strand::Forward)
            .unwrap()
            .unwrap();
        assert_eq!(result, "Met -> Ile");
    }

    #[test]
    fn test_node_display() {
        let var_node = Node {
            node_type: NodeType::Var("A".to_string()),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        let ref_node = Node {
            node_type: NodeType::Ref("".to_string()),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 99,
            index: 1,
        };

        assert_eq!(format!("{}", var_node), "Var(A) at position 42");
        assert_eq!(format!("{}", ref_node), "Ref at position 99");
    }
}
