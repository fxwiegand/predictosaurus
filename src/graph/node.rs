use crate::graph::EventProbs;
use crate::impact::Impact;
use crate::transcription;
use crate::translation::amino_acids::AminoAcid;
use anyhow::anyhow;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use rust_htslib::bcf::Record;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use varlociraptor::variants::evidence::observations::read_observation::ProcessedReadObservation;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) enum NodeType {
    Var(String),
    Ref(String),
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

impl Node {
    pub(crate) fn from_records(
        calls_record: &Record,
        _observations_record: &HashMap<&&String, Vec<ProcessedReadObservation>>,
        tags: &Vec<String>,
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
            probs: EventProbs::from_record(calls_record, tags),
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
    ) -> anyhow::Result<AminoAcid> {
        let p = match strand {
            Strand::Forward => self.pos,
            Strand::Reverse => self.position_on_reverse_strand(reference.len()),
            Strand::Unknown => return Err(anyhow!("Strand is unknown")),
        };
        let start_pos = p as usize - ((p - phase as i64) % 3) as usize;
        let ref_codon_bases = reference[start_pos..start_pos + 3].to_vec();
        AminoAcid::from_codon(
            transcription::transcribe_dna_to_rna(ref_codon_bases.as_ref())?.as_ref(),
        )
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
                let start_pos = p as usize - ((p - phase as i64) % 3) as usize;
                let position_in_codon = (p - phase as i64) % 3;
                let needed_bases = if alt_allele.is_empty() { 4 } else { 3 };
                let ref_codon_bases = reference[start_pos..start_pos + needed_bases].to_vec();
                let alt_codon_bases = [
                    &ref_codon_bases[..position_in_codon as usize],
                    alt_allele.as_bytes(),
                    &ref_codon_bases[position_in_codon as usize + 1..],
                ]
                .concat();
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
                let alt_amino_acids = self.variant_amino_acids(phase, reference, strand)?;
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
                    (false, _, _) => Ok(Impact::Modifier),
                };
                for amino_acid in alt_amino_acids {
                    if amino_acid == AminoAcid::Stop {
                        return Ok(Impact::High);
                    }
                }
                impact
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
            .unwrap();
        assert_eq!(met, AminoAcid::Methionine);
        let cys = node
            .reference_amino_acid(1, reference, Strand::Forward)
            .unwrap();
        assert_eq!(cys, AminoAcid::Cysteine);
        let ala = node
            .reference_amino_acid(2, reference, Strand::Forward)
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
            .unwrap();
        assert_eq!(his, AminoAcid::Histidine);
        let arg = node
            .reference_amino_acid(1, backward_reference, Strand::Reverse)
            .unwrap();
        assert_eq!(arg, AminoAcid::Arginine);
        let ala = node
            .reference_amino_acid(2, backward_reference, Strand::Reverse)
            .unwrap();
        assert_eq!(ala, AminoAcid::Alanine);
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
}
