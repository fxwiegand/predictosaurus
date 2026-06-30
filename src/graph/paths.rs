use crate::graph::VariantGraph;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::dna_to_amino_acids;
use crate::utils::fasta::reverse_complement;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use itertools::Itertools;
use petgraph::graph::NodeIndex;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

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
    pub(crate) fn contains_variant(
        &self,
        target: &str,
        variants: &HashMap<String, BTreeSet<i64>>,
    ) -> bool {
        if let Some(variant_positions) = variants.get(target) {
            variant_positions
                .range(self.start as i64..=self.end as i64)
                .next()
                .is_some()
        } else {
            false
        }
    }

    /// Length of CDS segment
    pub fn length(&self) -> usize {
        (self.end - self.start + 1) as usize
    }
}

impl HaplotypePath {
    /// Returns the overall frameshift caused by the path
    pub(crate) fn frameshift(&self, graph: &VariantGraph) -> i64 {
        self.0
            .iter()
            .map(|node_index| graph.graph.node_weight(*node_index).unwrap().frameshift())
            .sum()
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
            if node.node_type.is_variant() {
                let allele = node.alternative_allele.clone();
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
    
    use crate::graph::{Edge, EventProbs, VariantGraph};
    
    
    use petgraph::{Directed, Graph};
    use std::collections::HashMap;

    fn setup_protein_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node_1 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "AT".to_string(),
            vaf: HashMap::new(),
            probs: EventProbs(HashMap::new()),
            pos: 4,
            index: 0,
        });
        let node_2 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "T".to_string(),
            alternative_allele: "".to_string(),
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
        let path = &graph.paths()[0];
        let protein = path
            .protein(&graph, 0, b"ACGTTTGTTAG", Strand::Forward, 2, 10)
            .unwrap();
        assert_eq!(protein, dna_to_amino_acids(b"GTATTGTAG").unwrap());
    }

    #[test]
    fn protein_returns_correct_protein_sequence_for_reverse_strand() {
        let graph = setup_protein_graph();
        let path = &graph.reverse_paths()[0];
        let protein = path
            .protein(&graph, 0, b"CTAACAAATGCA", Strand::Reverse, 2, 10)
            .unwrap();
        assert_eq!(protein, dna_to_amino_acids(b"GCTTTATTT").unwrap());
    }

    #[test]
    fn protein_returns_correct_protein_sequence_for_reverse_strand_with_phase_offset() {
        let graph = setup_protein_graph();
        let path = &graph.reverse_paths()[0];
        let protein = path
            .protein(&graph, 1, b"AAAAAAAAAAAAAAAAAAAAAT", Strand::Reverse, 3, 21)
            .unwrap();
        assert_eq!(protein, dna_to_amino_acids(b"TTTTTTTTTTTTTTTATT").unwrap());
    }

    #[test]
    fn contains_variant_returns_true_when_variant_in_range() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([15]));
        assert!(cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_false_when_variant_out_of_range() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([25]));
        assert!(!cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_false_when_no_variants() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let variants: HashMap<String, BTreeSet<i64>> = HashMap::new();
        assert!(!cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_true_when_multiple_variants_in_range() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([12, 18]));
        assert!(cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_true_when_variants_on_boundary() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([10, 20]));
        assert!(cds.contains_variant(target, &variants));
    }
}
