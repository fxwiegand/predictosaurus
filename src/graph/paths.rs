use crate::graph::duck::feature_graph;
use crate::graph::{shift_phase, NodeType, VariantGraph};
use crate::impact::Impact;
use crate::translation::amino_acids::Protein;
use crate::translation::dna_to_protein;
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
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
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

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub(crate) struct Transcript {
    pub(crate) feature: String,
    pub(crate) target: String,
    pub(crate) strand: Strand,
    pub(crate) coding_sequences: Vec<Cds>,
}

impl Transcript {
    pub(crate) fn new(
        feature: String,
        target: String,
        strand: Strand,
        coding_sequences: Vec<Cds>,
    ) -> Transcript {
        Transcript {
            feature,
            target,
            strand,
            coding_sequences,
        }
    }

    pub(crate) fn name(&self) -> String {
        format!("{}:{}", self.feature, self.target)
    }

    pub(crate) fn start(&self) -> Result<u64> {
        self.coding_sequences
            .iter()
            .map(|cds| cds.start)
            .min()
            .ok_or_else(|| anyhow::anyhow!("No CDS found for transcript {}", self.name()))
    }

    pub(crate) fn end(&self) -> Result<u64> {
        self.coding_sequences
            .iter()
            .map(|cds| cds.end)
            .max()
            .ok_or_else(|| anyhow::anyhow!("No CDS found for transcript {}", self.name()))
    }

    /// Returns an iterator over the coding sequences of the transcript
    /// The iterator is sorted by start position in ascending order
    /// If the strand is reverse, the iterator is reversed
    pub(crate) fn cds(&self) -> Box<dyn Iterator<Item = &Cds> + '_> {
        match self.strand {
            Strand::Reverse => Box::new(
                self.coding_sequences
                    .iter()
                    .sorted_by_key(|cds| cds.start)
                    .rev(),
            ),
            _ => Box::new(self.coding_sequences.iter().sorted_by_key(|cds| cds.start)),
        }
    }

    pub(crate) fn weights(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
    ) -> Result<Vec<Vec<Weight>>> {
        let mut weights = Vec::new();
        for cds in self.cds() {
            if let Ok(graph) = feature_graph(
                graph.to_owned(),
                self.target.to_string(),
                cds.start,
                cds.end,
            ) {
                info!(
                    "Subgraph for CDS ({}-{}) of transcript {} has {} nodes",
                    cds.start,
                    cds.end,
                    self.name(),
                    graph.graph.node_count()
                );
                let paths = match self.strand {
                    Strand::Forward => Ok(graph.paths()),
                    Strand::Reverse => Ok(graph.reverse_paths()),
                    Strand::Unknown => Err(anyhow::anyhow!(
                        "Strand is unknown for transcript {}",
                        self.name()
                    )),
                }?;
                let reference_sequence = match self.strand {
                    Strand::Forward => reference.get(&self.target).unwrap(),
                    Strand::Reverse => &reverse_complement(reference.get(&self.target).unwrap()),
                    Strand::Unknown => {
                        unreachable!();
                    }
                };
                if weights.is_empty() {
                    let cds_weights = paths
                        .iter()
                        .map(|path| {
                            (
                                path.weights(&graph, cds.phase, reference_sequence, self.strand)
                                    .unwrap(),
                                path.frameshift(&graph),
                            )
                        })
                        .collect_vec();
                    for w in cds_weights {
                        weights.push(w);
                    }
                } else {
                    let mut new_weights = Vec::new();

                    for (accumulated_weights, accumulated_fs) in &weights {
                        let phase = shift_phase(cds.phase, ((*accumulated_fs + 3) % 3) as u8);
                        let cds_options = paths
                            .iter()
                            .map(|path| {
                                (
                                    path.weights(&graph, phase, reference_sequence, self.strand)
                                        .unwrap(),
                                    path.frameshift(&graph),
                                )
                            })
                            .collect_vec();

                        // For each option from the current CDS,
                        // create a new combination that appends its weights
                        // and adds its frameshift to the accumulated one.
                        // Offset weight index by max index in accumulated weights
                        for (mut new_weights_option, delta_fs) in cds_options {
                            let mut combined = accumulated_weights.clone();
                            let offset = accumulated_weights.iter().map(|w| w.index).max().unwrap();
                            for mut w in new_weights_option {
                                w.index += offset + 1;
                                combined.push(w);
                            }
                            new_weights.push((combined, *accumulated_fs + delta_fs));
                        }
                    }

                    // Replace the old combinations with the newly computed ones.
                    weights = new_weights;
                }
            }
            // print indexes of weights
            for (p, _) in &weights {
                println!("{:?}", p.iter().map(|w| w.index).collect::<Vec<usize>>());
            }
        }
        let weights = weights.iter().map(|(w, fs)| w.clone()).collect_vec();
        Ok(weights)
    }
}

pub(crate) fn transcripts(gff_file: &PathBuf) -> Result<Vec<Transcript>> {
    let mut feature_reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3)?;
    let mut transcripts = HashMap::new();
    for record in feature_reader
        .records()
        .filter_map(Result::ok)
        .filter(|record| record.feature_type() == "CDS")
    {
        let ensp = record.attributes().get("ID").ok_or_else(|| {
            anyhow::anyhow!("No ID found for CDS in sequence {}", record.seqname())
        })?;
        let target = record.seqname().to_string();
        let start = *record.start();
        let end = *record.end();
        let phase = record.phase().clone().try_into().unwrap();
        let strand = record.strand().ok_or_else(|| {
            anyhow::anyhow!("No strand found for CDS in sequence {}", record.seqname())
        })?;
        let cds = Cds::new(start, end, phase);
        let transcript = transcripts.entry(ensp.to_string()).or_insert_with(|| {
            Transcript::new(ensp.to_string(), target.clone(), strand, Vec::new())
        });
        transcript.coding_sequences.push(cds);
    }
    Ok(transcripts.into_values().collect())
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
    ) -> Result<Protein> {
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
        dna_to_protein(&sequence)
    }
}

mod tests {
    use crate::graph::node::{Node, NodeType};
    use crate::graph::paths::{Cds, Transcript};
    use crate::graph::{Edge, EventProbs, VariantGraph};
    use crate::impact::Impact;
    use crate::translation::dna_to_protein;
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
        let path = &graph.paths()[0];
        let protein = path
            .protein(&graph, 0, b"ACGTTTGTTAG", Strand::Forward, 2, 10)
            .unwrap();
        assert_eq!(protein, dna_to_protein(b"GTATTGTAG").unwrap());
    }

    #[test]
    fn protein_returns_correct_protein_sequence_for_reverse_strand() {
        let graph = setup_protein_graph();
        let path = &graph.reverse_paths()[0];
        let protein = path
            .protein(&graph, 0, b"CTAACAAATGCA", Strand::Reverse, 2, 10)
            .unwrap();
        assert_eq!(protein, dna_to_protein(b"GCTTTATTT").unwrap());
    }

    #[test]
    fn protein_returns_correct_protein_sequence_for_reverse_strand_with_phase_offset() {
        let graph = setup_protein_graph();
        let path = &graph.reverse_paths()[0];
        let protein = path
            .protein(&graph, 1, b"AAAAAAAAAAAAAAAAAAAAAT", Strand::Reverse, 3, 21)
            .unwrap();
        assert_eq!(protein, dna_to_protein(b"TTTTTTTTTTTTTTTATT").unwrap());
    }

    #[test]
    fn name_formats_transcript_correctly() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        assert_eq!(transcript.name(), "ENSP00000493376:test");
    }
}
