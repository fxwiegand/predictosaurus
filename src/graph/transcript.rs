use crate::graph::duck::feature_graph;
use crate::graph::node::NodeType;
use crate::graph::paths::{Cds, HaplotypePath, Weight};
use crate::graph::peptide::Peptide;
use crate::graph::{shift_phase, EventProbs, VariantGraph};
use crate::utils::fasta::reverse_complement;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use bio::io::gff;
use itertools::Itertools;
use log::info;
use rust_htslib::bgzf::CompressionLevel::Default;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

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

    pub(crate) fn start(&self) -> anyhow::Result<u64> {
        self.coding_sequences
            .iter()
            .map(|cds| cds.start)
            .min()
            .ok_or_else(|| anyhow::anyhow!("No CDS found for transcript {}", self.name()))
    }

    pub(crate) fn end(&self) -> anyhow::Result<u64> {
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

    fn paths(&self, graph: &VariantGraph) -> Result<Vec<HaplotypePath>> {
        match self.strand {
            Strand::Forward => Ok(graph.paths()),
            Strand::Reverse => Ok(graph.reverse_paths()),
            Strand::Unknown => Err(anyhow::anyhow!(
                "Strand is unknown for transcript {}",
                self.name()
            )),
        }
    }

    fn reference(&self, reference: &HashMap<String, Vec<u8>>) -> Result<Vec<u8>> {
        match self.strand {
            Strand::Forward => Ok(reference.get(&self.target).unwrap().to_vec()),
            Strand::Reverse => Ok(reverse_complement(reference.get(&self.target).unwrap())),
            Strand::Unknown => Err(anyhow::anyhow!(
                "Strand is unknown for transcript {}",
                self.name()
            )),
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
                let paths = self.paths(&graph)?;
                let reference_sequence = self.reference(reference)?;
                if weights.is_empty() {
                    let cds_weights = paths
                        .iter()
                        .map(|path| {
                            (
                                path.weights(&graph, cds.phase, &reference_sequence, self.strand)
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
                                    path.weights(&graph, phase, &reference_sequence, self.strand)
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
        }
        let weights = weights.iter().map(|(w, fs)| w.clone()).collect_vec();
        Ok(weights)
    }

    pub(crate) fn peptides() -> Result<Vec<Peptide>> {
        unimplemented!()
    }

    // Generate RNA sequences for the transcript, create all possible paths by combining the CDSs paths and then apply the variants to the RNA sequences
    pub(crate) fn rna(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
    ) -> Result<Vec<RnaPath>> {
        let reference = reference.get(&self.target).unwrap();
        let mut rna = Vec::new();
        for cds in self.cds() {
            let mut cds_rna = Vec::new();
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
                let mut sequence = reference[cds.start as usize..=cds.end as usize].to_vec();
                if self.strand == Strand::Reverse {
                    sequence = reverse_complement(&sequence);
                }
                sequence = sequence[cds.phase as usize..].to_vec();
                for path in self.paths(&graph)? {
                    let mut variants = HashMap::new();
                    let mut frameshift = 0;
                    let mut path_sequence = sequence.clone();
                    for node_index in path.0.iter() {
                        let node = graph.graph.node_weight(*node_index).unwrap();
                        if let NodeType::Var(allele) = &node.node_type {
                            let position_in_cds = match self.strand {
                                Strand::Forward => {
                                    (node.pos - cds.start as i64 + frameshift + cds.phase as i64)
                                        as usize
                                }
                                Strand::Reverse => {
                                    (cds.end as i64 - node.pos + frameshift - cds.phase as i64)
                                        as usize
                                }
                                Strand::Unknown => {
                                    return Err(anyhow::anyhow!("Strand is unknown"))
                                }
                            };
                            variants.insert(
                                position_in_cds,
                                (
                                    graph.graph.node_weight(*node_index).unwrap().probs.clone(),
                                    node.frameshift(),
                                ),
                            );
                            match self.strand {
                                Strand::Forward => {
                                    path_sequence.splice(
                                        position_in_cds..position_in_cds + 1,
                                        allele.bytes(),
                                    );
                                }
                                Strand::Reverse => {
                                    path_sequence.splice(
                                        position_in_cds..position_in_cds + 1,
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
                    let rna_path = RnaPath::new(path_sequence, variants);
                    cds_rna.push(rna_path);
                }
            }
            if rna.is_empty() {
                rna = cds_rna;
            } else {
                let mut new_rna = Vec::new();
                for rna_path in &rna {
                    for cds_rna_path in &cds_rna {
                        new_rna.push(rna_path.merge(cds_rna_path.clone())?);
                    }
                }
                rna = new_rna;
            }
        }
        Ok(rna)
    }
}

/// Represents an RNA sequence for a transcript and holds all variants on it
#[derive(Clone, Debug)]
pub(crate) struct RnaPath {
    pub(crate) sequence: Vec<u8>,
    pub(crate) variants: HashMap<usize, (EventProbs, i64)>,
}

impl RnaPath {
    pub(crate) fn new(sequence: Vec<u8>, variants: HashMap<usize, (EventProbs, i64)>) -> RnaPath {
        RnaPath { sequence, variants }
    }

    pub(crate) fn merge(&self, other: RnaPath) -> Result<RnaPath> {
        let mut sequence = self.sequence.clone();
        let mut variants = self.variants.clone();
        let offset = self.sequence.len();
        sequence.extend(other.sequence);
        for (position, probs) in other.variants {
            variants.insert(position + offset, probs);
        }
        Ok(RnaPath { sequence, variants })
    }
}

pub(crate) fn transcripts(gff_file: &PathBuf) -> anyhow::Result<Vec<Transcript>> {
    let mut feature_reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3)?;
    let mut transcripts = HashMap::new();
    for record in feature_reader
        .records()
        .filter_map(anyhow::Result::ok)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::duck::write_graphs;
    use crate::graph::node::Node;
    use crate::graph::Edge;
    use petgraph::{Directed, Graph};

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

    fn setup_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let alt_node_vaf_1 = HashMap::from([("s1".to_string(), 0.5)]);
        let alt_node_vaf_2 = HashMap::from([("s1".to_string(), 0.3)]);
        let alt_node_1 = graph.add_node(Node {
            node_type: NodeType::Var("".to_string()),
            vaf: alt_node_vaf_1,
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 0,
        });
        let alt_node_2 = graph.add_node(Node {
            node_type: NodeType::Var("G".to_string()),
            vaf: alt_node_vaf_2.clone(),
            probs: EventProbs(HashMap::new()),
            pos: 4,
            index: 1,
        });
        let alt_node_3 = graph.add_node(Node {
            node_type: NodeType::Var("A".to_string()),
            vaf: alt_node_vaf_2.clone(),
            probs: EventProbs(HashMap::new()),
            pos: 14,
            index: 3,
        });
        let ref_node_1 = graph.add_node(Node {
            node_type: NodeType::Ref("".to_string()),
            vaf: alt_node_vaf_2,
            probs: EventProbs(HashMap::new()),
            pos: 12,
            index: 2,
        });
        graph.add_edge(
            alt_node_1,
            alt_node_2,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graph.add_edge(
            ref_node_1,
            alt_node_3,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        VariantGraph {
            graph,
            start: 0,
            end: 15,
            target: "test".to_string(),
        }
    }

    #[test]
    fn rna_generates_correct_paths_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 10, 0), Cds::new(12, 15, 0)],
        );
        let graph = setup_graph();
        let tmp = tempfile::tempdir().unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(HashMap::from([("test".to_string(), graph)]), &graph_path).unwrap();
        let reference = HashMap::from([(
            "test".to_string(),
            vec![
                b'A', b'T', b'G', b'C', b'A', b'T', b'G', b'C', b'A', b'T', b'T', b'T', b'T', b'T',
                b'T', b'T', b'T', b'T',
            ],
        )]);
        let rna_paths = transcript.rna(&graph_path, &reference).unwrap();
        assert_eq!(rna_paths.len(), 1);
        assert_eq!(
            rna_paths[0].sequence,
            vec![
                b'A', b'G', b'C', b'G', b'T', b'G', b'C', b'A', b'T', b'T', b'T', b'T', b'A', b'T'
            ]
        );
    }
}
