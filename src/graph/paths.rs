use crate::graph::duck::feature_graph;
use crate::graph::{shift_phase, NodeType, VariantGraph};
use crate::impact::Impact;
use crate::translation::amino_acids::{AminoAcid, Protein};
use crate::translation::dna_to_amino_acids;
use crate::utils;
use crate::utils::fasta::reverse_complement;
use anyhow::{anyhow, Result};
use bio::bio_types::strand::Strand;
use bio::io::gff;
use itertools::Itertools;
use log::info;
use petgraph::graph::NodeIndex;
use serde::{Deserialize, Serialize};
use std::cmp::max;
use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub(crate) struct Weight {
    pub(crate) index: usize,
    pub(crate) path: Option<usize>,
    pub(crate) vaf: f32,
    pub(crate) impact: Impact,
    pub(crate) reason: Option<String>,
    pub(crate) consequence: Option<String>,
    pub(crate) sample: String,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq, Hash)]
pub(crate) struct Cds {
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) phase: u8,
}

impl Cds {
    pub(crate) fn new(start: u64, end: u64, phase: u8) -> Cds {
        Cds { start, end, phase }
    }

    /// Returns true if the CDS contains the variant
    pub(crate) fn contains_variant(&self, variants: &HashMap<String, BTreeSet<i64>>) -> bool {
        for variant in variants.values() {
            if variant
                .range(self.start as i64..=self.end as i64)
                .next()
                .is_some()
            {
                return true;
            }
        }
        false
    }

    /// Length of CDS segment
    pub fn length(&self) -> usize {
        (self.end - self.start + 1) as usize
    }
}

impl HaplotypePath {
    /// Returns the maximum impact of the individual nodes on the path
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

    /// Returns the overall frameshift caused by the path
    pub(crate) fn frameshift(&self, graph: &VariantGraph) -> i64 {
        self.0
            .iter()
            .map(|node_index| graph.graph.node_weight(*node_index).unwrap().frameshift())
            .sum()
    }

    pub(crate) fn display(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> Result<String> {
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

    pub(crate) fn protein(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
        strand: Strand,
        start: usize,
        end: usize,
    ) -> Result<Vec<AminoAcid>> {
        let mut frameshift = 0;
        let mut sequence = reference[start..=end].to_vec();
        if strand == Strand::Reverse {
            sequence = reverse_complement(&sequence);
        }
        sequence = sequence[phase as usize..].to_vec();

        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            if let NodeType::Var(allele) = &node.node_type {
                let position_in_protein = match strand {
                    Strand::Forward => {
                        (node.pos - start as i64 + frameshift + phase as i64) as usize
                    }
                    Strand::Reverse => (end as i64 - node.pos + frameshift - phase as i64) as usize,
                    Strand::Unknown => return Err(anyhow::anyhow!("Strand is unknown")),
                };
                match strand {
                    Strand::Forward => {
                        sequence
                            .splice(position_in_protein..position_in_protein + 1, allele.bytes());
                    }
                    Strand::Reverse => {
                        sequence.splice(
                            position_in_protein..position_in_protein + 1,
                            String::from_utf8_lossy(
                                reverse_complement(allele.as_bytes()).as_slice(),
                            )
                            .bytes(),
                        );
                    }
                    Strand::Unknown => unreachable!(),
                }

                frameshift += node.frameshift();
            }
        }
        dna_to_amino_acids(&sequence)
    }
}

mod tests {
    use crate::graph::node::{Node, NodeType};
    use crate::graph::paths::Cds;
    use crate::graph::{Edge, EventProbs, VariantGraph};
    use crate::impact::Impact;
    use crate::translation::dna_to_amino_acids;
    use bio::bio_types::strand::Strand;
    use petgraph::{Directed, Graph};
    use std::collections::{BTreeSet, HashMap};

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
        let path = graph.paths().next().unwrap();
        let reference = b"ACGTTTGTTA";
        let strand = Strand::Forward;
        let weights = path.weights(&graph, 0, reference, strand).unwrap();
        assert_eq!(weights.len(), 4);
        assert_eq!(weights[0].vaf, 0.5);
        assert_eq!(weights[2].impact, Impact::Moderate);
    }

    fn setup_protein_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node_1 = graph.add_node(Node {
            node_type: NodeType::Var("AT".to_string()),
            vaf: HashMap::new(),
            probs: EventProbs(HashMap::new()),
            pos: 4,
            index: 0,
        });
        let node_2 = graph.add_node(Node {
            node_type: NodeType::Var("".to_string()),
            vaf: HashMap::new(),
            probs: EventProbs(HashMap::new()),
            pos: 8,
            index: 1,
        });
        graph.add_edge(
            node_1,
            node_2,
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
    fn protein_returns_correct_protein_sequence() {
        let graph = setup_protein_graph();
        let path = graph.paths().next().unwrap();
        let protein = path
            .protein(&graph, 0, b"ACGTTTGTTAG", Strand::Forward, 2, 10)
            .unwrap();
        assert_eq!(protein, dna_to_amino_acids(b"GTATTGTAG").unwrap());
    }

    #[test]
    fn protein_returns_correct_protein_sequence_for_reverse_strand() {
        let graph = setup_protein_graph();
        let path = graph.reverse_paths().next().unwrap();
        let protein = path
            .protein(&graph, 0, b"CTAACAAATGCA", Strand::Reverse, 2, 10)
            .unwrap();
        assert_eq!(protein, dna_to_amino_acids(b"GCTTTATTT").unwrap());
    }

    #[test]
    fn protein_returns_correct_protein_sequence_for_reverse_strand_with_phase_offset() {
        let graph = setup_protein_graph();
        let path = graph.reverse_paths().next().unwrap();
        let protein = path
            .protein(&graph, 1, b"AAAAAAAAAAAAAAAAAAAAAT", Strand::Reverse, 3, 21)
            .unwrap();
        assert_eq!(protein, dna_to_amino_acids(b"TTTTTTTTTTTTTTTATT").unwrap());
    }

    #[test]
    fn contains_variant_returns_true_when_variant_in_range() {
        let cds = Cds::new(10, 20, 0);
        let mut variants = HashMap::new();
        variants.insert("variant1".to_string(), BTreeSet::from([15]));
        assert!(cds.contains_variant(&variants));
    }

    #[test]
    fn contains_variant_returns_false_when_variant_out_of_range() {
        let cds = Cds::new(10, 20, 0);
        let mut variants = HashMap::new();
        variants.insert("variant1".to_string(), BTreeSet::from([25]));
        assert!(!cds.contains_variant(&variants));
    }

    #[test]
    fn contains_variant_returns_false_when_no_variants() {
        let cds = Cds::new(10, 20, 0);
        let variants: HashMap<String, BTreeSet<i64>> = HashMap::new();
        assert!(!cds.contains_variant(&variants));
    }

    #[test]
    fn contains_variant_returns_true_when_multiple_variants_in_range() {
        let cds = Cds::new(10, 20, 0);
        let mut variants = HashMap::new();
        variants.insert("variant1".to_string(), BTreeSet::from([12, 18]));
        assert!(cds.contains_variant(&variants));
    }

    #[test]
    fn contains_variant_returns_true_when_variants_on_boundary() {
        let cds = Cds::new(10, 20, 0);
        let mut variants = HashMap::new();
        variants.insert("variant1".to_string(), BTreeSet::from([10, 20]));
        assert!(cds.contains_variant(&variants));
    }
}
