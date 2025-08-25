use crate::cli::Interval;
use crate::graph::duck::{feature_graph, variants_on_graph};
use crate::graph::node::{Node, NodeType};
use crate::graph::paths::{Cds, HaplotypePath, Weight};
use crate::graph::peptide::Peptide;
use crate::graph::score::EffectScore;
use crate::graph::{shift_phase, EventProbs, VariantGraph};
use crate::translation::amino_acids::AminoAcid;
use crate::utils::fasta::reverse_complement;
use anyhow::{bail, Result};
use bio::bio_types::strand::Strand;
use bio::io::gff::{self, Phase};
use bio::stats::LogProb;
use itertools::Itertools;
use log::info;
use petgraph::adj::NodeIndex;
use rust_htslib::bgzf::CompressionLevel::Default;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

use super::node;

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
        format!("{}:{}", self.target, self.feature)
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

    // Returns the total length of all CDS of the transcript
    pub(crate) fn length(&self) -> usize {
        self.cds().map(|c| c.length()).sum()
    }

    // Returns the position of the variant in relation to the overall length of transcript
    pub(crate) fn position_in_transcript(&self, pos: usize) -> Result<usize> {
        let mut offset = 0;
        for cds in self.cds() {
            if ((cds.start as usize)..=(cds.end as usize)).contains(&pos) {
                return Ok(offset + (pos - cds.start as usize));
            }
            offset += cds.length();
        }
        bail!(
            "Position {} not found in CDS of transcript {}",
            pos,
            self.name()
        );
    }

    /// Returns the CDS containing the given genomic position, or None if not found.
    pub(crate) fn cds_for_position(&self, pos: i64) -> Option<&Cds> {
        self.cds()
            .find(|cds| (cds.start as i64..=cds.end as i64).contains(&pos))
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
        let sequence = reference.get(&self.target).ok_or_else(|| {
            anyhow::anyhow!("Reference sequence not found for target {}", self.target)
        })?;
        match self.strand {
            Strand::Forward => Ok(sequence.to_vec()),
            Strand::Reverse => Ok(reverse_complement(sequence)),
            Strand::Unknown => Err(anyhow::anyhow!(
                "Strand is unknown for transcript {}",
                self.name()
            )),
        }
    }

    fn haplotypes(&self, graph: &PathBuf) -> Result<Vec<Vec<Node>>> {
        let mut haplotypes = Vec::new();
        for cds in self.cds() {
            if let Ok(graph) = feature_graph(
                graph.to_owned(),
                self.target.to_string(),
                cds.start,
                cds.end,
            ) {
                info!(
                    "Calculating paths for {cds:?}. Graph has {} nodes",
                    graph.graph.node_count()
                );
                let paths = self
                    .paths(&graph)?
                    .iter()
                    .map(|p| {
                        p.0.iter()
                            .map(|v| graph.graph.node_weight(*v).unwrap().to_owned())
                            .collect_vec()
                    })
                    .collect_vec();
                if haplotypes.is_empty() {
                    haplotypes = paths;
                } else {
                    let mut extended = Vec::new();

                    for (existing_path, new_path) in
                        haplotypes.iter().cartesian_product(paths.iter())
                    {
                        let mut combined = existing_path.clone();
                        combined.extend(new_path.iter().cloned());
                        extended.push(combined);
                    }
                    haplotypes = extended;
                }
            }
        }
        Ok(haplotypes)
    }

    pub(crate) fn scores(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
    ) -> Result<Vec<EffectScore>> {
        let haplotypes = self.haplotypes(graph)?;
        let mut scores = Vec::with_capacity(haplotypes.len());
        for haplotype in haplotypes {
            scores.push(EffectScore::from_haplotype(reference, self, haplotype)?);
        }
        Ok(scores)
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
                info!(
                    "Found {} paths for CDS ({}-{}) of transcript {}",
                    paths.len(),
                    cds.start,
                    cds.end,
                    self.name()
                );
                if paths.is_empty() {
                    continue;
                }
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
        let weights = weights.into_iter().map(|(w, _)| w).collect_vec();
        Ok(weights)
    }

    pub(crate) fn peptides(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
        interval: Interval,
        sample: &str,
        events: &[String],
        min_event_prob: LogProb,
        background_events: &[String],
        max_background_event_prob: LogProb,
    ) -> Result<Vec<Peptide>> {
        let rnas = self.rna(graph, reference, sample)?;
        let mut peptides = Vec::new();
        for i in interval {
            for rna in &rnas {
                let p = rna.peptides(i)?;
                for peptide in p {
                    peptides.push(peptide);
                }
            }
        }
        peptides = peptides
            .iter()
            .filter(|peptide| {
                peptide.prob(events).unwrap() >= min_event_prob
                    && peptide.prob(background_events).unwrap() <= max_background_event_prob
            })
            .cloned()
            .collect_vec();
        Ok(peptides)
    }

    // Generate RNA sequences for the transcript, create all possible paths by combining the CDSs paths and then apply the variants to the RNA sequences
    pub(crate) fn rna(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
        sample: &str,
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
                                    *node.vaf.get(sample).unwrap(),
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
                    let rna_path = RnaPath::new(path_sequence, variants, self.clone());
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
    pub(crate) variants: HashMap<usize, (EventProbs, i64, f32)>,
    pub(crate) transcript: Transcript,
}

impl RnaPath {
    pub(crate) fn new(
        sequence: Vec<u8>,
        variants: HashMap<usize, (EventProbs, i64, f32)>,
        transcript: Transcript,
    ) -> RnaPath {
        RnaPath {
            sequence,
            variants,
            transcript,
        }
    }

    pub(crate) fn merge(&self, other: RnaPath) -> Result<RnaPath> {
        assert_eq!(
            self.transcript, other.transcript,
            "Transcripts must be the same to merge both paths"
        );
        let mut sequence = self.sequence.clone();
        let mut variants = self.variants.clone();
        let offset = self.sequence.len();
        sequence.extend(other.sequence);
        for (position, probs) in other.variants {
            variants.insert(position + offset, probs);
        }
        Ok(RnaPath {
            sequence,
            variants,
            transcript: self.transcript.clone(),
        })
    }

    pub(crate) fn peptides(&self, length: u32) -> Result<Vec<Peptide>> {
        let mut peptides = Vec::new();
        for (i, p) in self
            .sequence
            .windows(length as usize * 3)
            .step_by(3)
            .enumerate()
        {
            let interval = (i * 3, i * 3 + length as usize * 3);
            // We consider a variant to influence the peptide if it is in the interval or if it is before the interval and causes a frameshift not equal to 0
            let probs_variants_in_interval = self
                .variants
                .iter()
                .filter(|(pos, (_, fs, _))| {
                    (*pos >= &interval.0 && *pos < &interval.1) || (*pos < &interval.0 && fs != &0)
                })
                .map(|(_, (probs, _, af))| (probs, *af))
                .collect_vec();
            // We skip peptides that have no variants in the interval or a frameshift upstream as they are not of any interest
            if probs_variants_in_interval.is_empty() {
                continue;
            }
            let mut probs = HashMap::new();
            for (variant, _) in &probs_variants_in_interval {
                for (event, prob) in variant.0.iter() {
                    // Add event or multiply probability if it already exists. We are summing because we are in log space
                    match probs.get_mut(event) {
                        Some(p) => *p += *prob,
                        None => {
                            probs.insert(event.clone(), *prob);
                        }
                    }
                }
            }
            let afs = probs_variants_in_interval
                .iter()
                .map(|(_, af)| *af)
                .collect();
            let peptide =
                Peptide::from_rna(p.to_vec(), EventProbs(probs), afs, self.transcript.clone())?;
            peptides.push(peptide);
        }
        Ok(peptides)
    }
}

pub(crate) fn transcripts(gff_file: &PathBuf, graph: &PathBuf) -> Result<Vec<Transcript>> {
    let variants = variants_on_graph(graph)?;
    let mut feature_reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3)?;
    let mut transcripts = HashMap::new();
    let mut variant_coverage = HashMap::new();
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
        if !variant_coverage.get(ensp).copied().unwrap_or(false) {
            let has_variant = cds.contains_variant(&variants);
            if has_variant {
                variant_coverage.insert(ensp.to_string(), true);
            }
        }
        let transcript = transcripts.entry(ensp.to_string()).or_insert_with(|| {
            Transcript::new(ensp.to_string(), target.clone(), strand, Vec::new())
        });
        transcript.coding_sequences.push(cds);
    }
    Ok(transcripts
        .into_values()
        .filter(|t| matches!(variant_coverage.get(&t.feature), Some(true)))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::duck::write_graphs;
    use crate::graph::node::Node;
    use crate::graph::Edge;
    use crate::translation::amino_acids::AminoAcid;
    use bio::stats::Prob;
    use petgraph::{Directed, Graph};

    #[test]
    fn name_formats_transcript_correctly() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        assert_eq!(transcript.name(), "chr1:ENSP00000493376");
    }

    fn setup_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let alt_node_vaf_1 = HashMap::from([("s1".to_string(), 0.5)]);
        let alt_node_vaf_2 = HashMap::from([("s1".to_string(), 0.3)]);
        let event_probs = EventProbs(HashMap::from([
            ("PROB_GERMLINE".to_string(), LogProb::from(Prob(0.3))),
            ("PROB_SOMATIC".to_string(), LogProb::from(Prob(0.8))),
        ]));
        let alt_node_1 = graph.add_node(Node {
            node_type: NodeType::Var("".to_string()),
            vaf: alt_node_vaf_1,
            probs: event_probs.clone(),
            pos: 1,
            index: 0,
        });
        let alt_node_2 = graph.add_node(Node {
            node_type: NodeType::Var("G".to_string()),
            vaf: alt_node_vaf_2.clone(),
            probs: event_probs.clone(),
            pos: 4,
            index: 1,
        });
        let alt_node_3 = graph.add_node(Node {
            node_type: NodeType::Var("A".to_string()),
            vaf: alt_node_vaf_2.clone(),
            probs: event_probs.clone(),
            pos: 14,
            index: 3,
        });
        let ref_node_1 = graph.add_node(Node {
            node_type: NodeType::Ref("".to_string()),
            vaf: alt_node_vaf_2,
            probs: event_probs.clone(),
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
            target: "chr1".to_string(),
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
        let rna_paths = transcript.rna(&graph_path, &reference, "s1").unwrap();
        assert_eq!(rna_paths.len(), 4);
        assert_eq!(
            rna_paths[0].sequence,
            vec![
                b'A', b'G', b'C', b'G', b'T', b'G', b'C', b'A', b'T', b'T', b'T', b'T', b'A', b'T'
            ]
        );
    }

    #[test]
    fn peptides_generate_correctly_for_valid_input() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 10, 0)],
        );
        let graph = setup_graph();
        let tmp = tempfile::tempdir().unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(HashMap::from([("test".to_string(), graph)]), &graph_path).unwrap();
        let reference = HashMap::from([(
            "test".to_string(),
            vec![
                b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T',
            ],
        )]);
        let interval = Interval { start: 2, end: 2 };
        let events = vec!["SOMATIC".to_string()];
        let background_events = vec!["GERMLINE".to_string()];
        let peptides = transcript
            .peptides(
                &graph_path,
                &reference,
                interval,
                "s1",
                &events,
                LogProb::from(Prob(0.4)),
                &background_events,
                LogProb::from(Prob(0.1)),
            )
            .unwrap();
        let peptide_sequences = peptides
            .iter()
            .map(|p| p.sequence.clone())
            .unique()
            .collect::<Vec<Vec<AminoAcid>>>();
        let expected = vec![
            vec![AminoAcid::Phenylalanine, AminoAcid::Valine],
            vec![AminoAcid::Valine, AminoAcid::Phenylalanine],
        ];
        assert_eq!(peptide_sequences, expected);
    }

    #[test]
    fn transcripts_parses_gff_file_correctly() {
        let gff_content = "\
            ##gff-version 3
            chr1\tsource\tCDS\t1\t100\t.\t+\t0\tID=ENSP00000493376
            chr1\tsource\tCDS\t200\t300\t.\t+\t0\tID=ENSP00000493376
            chr1\tsource\tCDS\t400\t500\t.\t-\t0\tID=ENSP00000493377
        ";
        let tmp = tempfile::tempdir().unwrap();
        let gff_path = tmp.path().join("test.gff");
        std::fs::write(&gff_path, gff_content).unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(
            HashMap::from([("chr1".to_string(), setup_graph())]),
            &graph_path,
        )
        .unwrap();

        let transcripts = transcripts(&gff_path, &graph_path).unwrap();
        assert_eq!(transcripts.len(), 1);

        let transcript1 = transcripts
            .iter()
            .find(|t| t.feature == "ENSP00000493376")
            .unwrap();
        assert_eq!(transcript1.coding_sequences.len(), 2);
        assert_eq!(transcript1.strand, Strand::Forward);
    }

    #[test]
    fn reference_returns_correct_sequence_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        let reference = HashMap::from([("test".to_string(), vec![b'A', b'T', b'G', b'C'])]);
        let result = transcript.reference(&reference).unwrap();
        assert_eq!(result, vec![b'A', b'T', b'G', b'C']);
    }

    #[test]
    fn reference_returns_reverse_complement_for_reverse_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Reverse,
            vec![Cds::new(1, 10, 0)],
        );
        let reference = HashMap::from([("test".to_string(), vec![b'A', b'T', b'G', b'C'])]);
        let result = transcript.reference(&reference).unwrap();
        assert_eq!(result, vec![b'G', b'C', b'A', b'T']);
    }

    #[test]
    fn start_returns_min_cds_start() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(5, 15, 0)],
        );
        assert_eq!(transcript.start().unwrap(), 5);
    }

    #[test]
    fn end_returns_max_cds_end() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(5, 25, 0)],
        );
        assert_eq!(transcript.end().unwrap(), 25);
    }

    #[test]
    fn cds_returns_sorted_cds_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(5, 15, 0)],
        );
        let cds: Vec<&Cds> = transcript.cds().collect();
        assert_eq!(cds.len(), 2);
        assert_eq!(cds[0].start, 5);
        assert_eq!(cds[1].start, 10);
    }

    #[test]
    fn cds_returns_sorted_cds_for_reverse_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Reverse,
            vec![Cds::new(10, 20, 0), Cds::new(5, 15, 0)],
        );
        let cds: Vec<&Cds> = transcript.cds().collect();
        assert_eq!(cds.len(), 2);
        assert_eq!(cds[0].start, 10);
        assert_eq!(cds[1].start, 5);
    }

    #[test]
    fn cds_returns_empty_iterator_when_no_cds() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![],
        );
        let cds: Vec<&Cds> = transcript.cds().collect();
        assert!(cds.is_empty());
    }

    #[test]
    fn test_transcript_length() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(30, 40, 0)],
        );
        assert_eq!(transcript.length(), 22)
    }

    #[test]
    fn paths_returns_paths_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let paths = transcript.paths(&graph).unwrap();
        assert_eq!(paths, graph.paths());
    }

    #[test]
    fn paths_returns_reverse_paths_for_reverse_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Reverse,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let paths = transcript.paths(&graph).unwrap();
        assert_eq!(paths, graph.reverse_paths());
    }

    #[test]
    fn paths_returns_error_for_unknown_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "chr1".to_string(),
            Strand::Unknown,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let result = transcript.paths(&graph);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Strand is unknown for transcript chr1:ENSP00000493376"
        );
    }

    #[test]
    fn weights_returns_correct_weights_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let tmp = tempfile::tempdir().unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(HashMap::from([("test".to_string(), graph)]), &graph_path).unwrap();
        let reference = HashMap::from([("test".to_string(), vec![b'A', b'T', b'G', b'C'])]);
        let weights = transcript.weights(&graph_path, &reference).unwrap();
        assert_eq!(weights.len(), 2);
        assert_eq!(weights[0].len(), 2);
    }

    #[test]
    fn weights_returns_correct_weights_for_reverse_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Reverse,
            vec![Cds::new(1, 10, 0), Cds::new(12, 15, 0)],
        );
        let graph = setup_graph();
        let tmp = tempfile::tempdir().unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(HashMap::from([("test".to_string(), graph)]), &graph_path).unwrap();
        let reference = HashMap::from([(
            "test".to_string(),
            vec![
                b'A', b'T', b'G', b'C', b'A', b'T', b'G', b'C', b'A', b'T', b'G', b'C', b'A', b'T',
                b'G', b'C',
            ],
        )]);
        let weights = transcript.weights(&graph_path, &reference).unwrap();
        assert_eq!(weights.len(), 4);
        assert_eq!(weights[0].len(), 3);
    }

    #[test]
    fn test_position_in_transcript() {
        let transcript = Transcript {
            feature: "Test".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![
                Cds {
                    start: 100,
                    end: 104,
                    phase: 0,
                },
                Cds {
                    start: 200,
                    end: 202,
                    phase: 0,
                },
            ],
        };

        assert_eq!(transcript.position_in_transcript(100).unwrap(), 0);
        assert_eq!(transcript.position_in_transcript(104).unwrap(), 4);
        assert_eq!(transcript.position_in_transcript(200).unwrap(), 5);
        assert_eq!(transcript.position_in_transcript(202).unwrap(), 7);
        assert!(transcript.position_in_transcript(150).is_err());
    }

    #[test]
    fn test_cds_for_position_option() {
        let transcript = Transcript {
            feature: "Test".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![
                Cds {
                    start: 100,
                    end: 199,
                    phase: 0,
                },
                Cds {
                    start: 300,
                    end: 399,
                    phase: 0,
                },
            ],
        };

        assert_eq!(
            transcript
                .cds_for_position(150)
                .map(|cds| (cds.start, cds.end)),
            Some((100, 199))
        );
        assert_eq!(
            transcript
                .cds_for_position(350)
                .map(|cds| (cds.start, cds.end)),
            Some((300, 399))
        );
        assert_eq!(transcript.cds_for_position(250), None);
    }
}
