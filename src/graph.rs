use crate::cli::ObservationFile;
use crate::impact::Impact;
use crate::transcription;
use crate::translation::amino_acids::AminoAcid;
use crate::utils::bcf::extract_event_names;
use anyhow::Result;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;
use petgraph::dot::{Config, Dot};
use petgraph::graph::NodeIndex;
use petgraph::{Directed, Graph};
use rust_htslib::bcf::{Read, Reader, Record};
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use varlociraptor::calling::variants::preprocessing::read_observations;
use varlociraptor::utils::collect_variants::collect_variants;
use varlociraptor::variants::evidence::observations::read_observation::ProcessedReadObservation;

pub(crate) struct VariantGraph {
    pub(crate) graph: Graph<Node, Edge, Directed>,
    start: i64,
    end: i64,
    target: String,
}

impl VariantGraph {
    pub(crate) fn build(
        calls_file: &PathBuf,
        observation_files: &[ObservationFile],
        target: &str,
        start: i64,
        end: i64,
    ) -> Result<VariantGraph> {
        let mut calls_reader = Reader::from_path(calls_file)?;
        let header = calls_reader.header().clone();
        let mut observation_readers: HashMap<_, _> = observation_files
            .iter()
            .map(|o| (o.sample.to_string(), Reader::from_path(&o.path).unwrap()))
            .collect();
        let mut observations_records = observation_readers
            .iter_mut()
            .map(|(sample, reader)| (sample.clone(), reader.records()))
            .collect::<HashMap<_, _>>();

        let samples = calls_reader
            .header()
            .samples()
            .iter()
            .map(|s| String::from_utf8(s.to_vec()).unwrap())
            .collect_vec();
        for sample in &samples {
            assert!(
                observation_files.iter().any(|o| o.sample == *sample),
                "Sample {} in calls file not found in observation files",
                sample
            );
        }

        let event_names = extract_event_names(calls_file);

        let tags = event_names
            .iter()
            .map(|event| format!("PROB_{event}"))
            .collect();

        let mut supporting_reads = HashMap::new();

        let mut variant_graph = Graph::<Node, Edge, Directed>::new();
        let mut index = 0;
        let mut last_position = -1;

        for calls_record in calls_reader.records() {
            let mut calls_record = calls_record?;
            let position = calls_record.pos();
            let mut observations_records = observations_records
                .iter_mut()
                .map(|(sample, records)| {
                    let record = records.next().unwrap().unwrap();
                    (sample, record)
                })
                .collect::<HashMap<_, _>>();
            let _variants = collect_variants(&mut calls_record, false, None)?;
            let observations = observations_records
                .iter_mut()
                .map(|(sample, record)| {
                    let observations = read_observations(record).unwrap();
                    (sample, observations.pileup.read_observations().clone())
                })
                .collect::<HashMap<_, _>>();

            if header.rid2name(calls_record.rid().unwrap()).unwrap() != target.as_bytes()
                || position < start
                || position > end
            {
                continue;
            }

            let _fragment_ids: HashSet<_> = observations
                .iter()
                .flat_map(|(s, v)| v.iter().map(move |o| (s, o.fragment_id)))
                .collect();

            let alleles = calls_record.alleles();
            let ref_allele = String::from_utf8(alleles[0].to_vec())?;
            let alt_allele = String::from_utf8(alleles[1].to_vec())?;

            let var_node = Node::from_records(
                &calls_record,
                &observations,
                &tags,
                NodeType::Var(alt_allele),
                &samples,
                index,
            );
            let var_node_index = variant_graph.add_node(var_node);
            let mut ref_node_index = None;
            if last_position != position {
                let ref_node = Node::from_records(
                    &calls_record,
                    &observations,
                    &tags,
                    NodeType::Ref(ref_allele),
                    &samples,
                    index,
                );
                ref_node_index = Some(variant_graph.add_node(ref_node));
                index += 1;
            }

            for (sample, observations) in observations {
                for observation in observations {
                    let evidence = BayesFactor::new(observation.prob_alt, observation.prob_ref)
                        .evidence_kass_raftery();
                    match evidence {
                        KassRaftery::Strong | KassRaftery::VeryStrong => {
                            // Read supports variant
                            let entry = supporting_reads
                                .entry((sample.to_string(), observation.fragment_id))
                                .or_insert(Vec::new());
                            entry.push(var_node_index);
                        }
                        KassRaftery::None | KassRaftery::Barely | KassRaftery::Positive => {
                            // Read supports reference
                            let entry = supporting_reads
                                .entry((sample.to_string(), observation.fragment_id))
                                .or_insert(Vec::new());
                            if let Some(index) = ref_node_index {
                                entry.push(index);
                            }
                        }
                    }
                }
            }

            last_position = position;
        }

        let mut variant_graph = VariantGraph {
            graph: variant_graph,
            start,
            end: last_position,
            target: target.to_string(),
        };

        variant_graph.create_edges(&supporting_reads)?;

        Ok(variant_graph)
    }

    pub(crate) fn create_edges(
        &mut self,
        supporting_reads: &HashMap<(String, Option<u64>), Vec<NodeIndex>>,
    ) -> Result<()> {
        for ((sample, _), nodes) in supporting_reads {
            let weights = self.graph.clone();
            for node_tuple in nodes.iter().sorted().dedup().combinations(2).filter(|v| {
                node_distance(
                    &weights.node_weight(*v[0]).unwrap().index,
                    &weights.node_weight(*v[1]).unwrap().index,
                ) <= 1
                    || nodes_in_between(
                        &weights.node_weight(*v[0]).unwrap().index,
                        &weights.node_weight(*v[1]).unwrap().index,
                        nodes
                            .iter()
                            .map(|n| &weights.node_weight(*n).unwrap().index)
                            .collect_vec(),
                    ) == 0
            }) {
                let edge = self.graph.find_edge(*node_tuple[0], *node_tuple[1]);
                if let Some(edge) = edge {
                    let edge = self.graph.edge_weight_mut(edge).unwrap();
                    edge.supporting_reads
                        .entry(sample.to_string())
                        .and_modify(|v| *v += 1)
                        .or_insert(1);
                } else {
                    let edge = Edge {
                        supporting_reads: HashMap::from([(sample.to_string(), 1)]),
                    };
                    self.graph.add_edge(*node_tuple[0], *node_tuple[1], edge);
                }
            }
        }
        Ok(())
    }

    pub(crate) fn to_dot(&self) -> String {
        format!(
            "digraph {{ {:?} }}",
            Dot::with_config(&self.graph, &[Config::GraphContentOnly])
        )
    }

    pub(crate) fn to_file(&self, path: &Path) -> Result<()> {
        let path = path.join("graph.dot");
        let mut file = std::fs::File::create(path)?;
        file.write_all(self.to_dot().as_bytes())?;
        Ok(())
    }

    /// Finds all paths starting from the first two nodes in the graph.
    ///
    /// This method performs a depth-first search (DFS) starting from the first two nodes
    /// in the graph. It collects all possible paths, including those that reach leaf nodes
    /// (nodes with no outgoing edges). The paths are represented as vectors of `NodeIndex`.
    ///
    /// Returns:
    ///     A vector of vectors, where each inner vector represents a path through the graph
    ///     starting from one of the initial nodes. Each path is a sequence of `NodeIndex`
    ///     elements, representing the nodes in the order they are visited
    pub(crate) fn paths(&self) -> Vec<HaplotypePath> {
        let mut all_paths = Vec::new();
        let start_nodes = self.graph.node_indices().take(2).collect::<Vec<_>>();

        for start_node in start_nodes {
            let mut stack = vec![(start_node, vec![start_node])];

            while let Some((node, path)) = stack.pop() {
                for neighbor in self.graph.neighbors(node) {
                    let mut new_path = path.clone();
                    new_path.push(neighbor);
                    stack.push((neighbor, new_path.clone()));

                    // Check if the neighbor is a leaf node (no outgoing edges)
                    if self.graph.edges(neighbor).next().is_none() {
                        all_paths.push(HaplotypePath(new_path));
                    }
                }
            }
        }

        all_paths
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.graph.node_count() == 0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

fn shift_phase(phase: u8, frameshift: u8) -> u8 {
    match phase {
        0 => [0, 2, 1][frameshift as usize],
        1 => [1, 0, 2][frameshift as usize],
        2 => [2, 1, 0][frameshift as usize],
        _ => unreachable!(),
    }
}

impl HaplotypePath {
    pub(crate) fn impact(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
    ) -> Result<Impact> {
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

    pub(crate) fn weight(&self, graph: &VariantGraph) -> f32 {
        self.0
            .iter()
            .filter(|n| graph.graph.node_weight(**n).unwrap().node_type.is_variant())
            .map(|n| {
                let node = graph.graph.node_weight(*n).unwrap();
                let vaf_sum: f32 = node.vaf.values().sum();
                let weight = vaf_sum / node.vaf.values().len() as f32;
                weight
            })
            .product()
    }

    pub(crate) fn display(
        &self,
        graph: &VariantGraph,
        phase: u8,
        reference: &[u8],
    ) -> Result<String> {
        let ref_phase = phase;
        let mut phase = phase;
        let mut protein = String::new();
        for node_index in self.0.iter() {
            let node = graph.graph.node_weight(*node_index).unwrap();
            let ref_amino_acid = node.reference_amino_acid(ref_phase, reference)?;
            let alt_amino_acid = node.variant_amino_acid(phase, reference)?;
            protein.push_str(&format!(
                "{} -> {} ({:?})\n",
                ref_amino_acid,
                alt_amino_acid,
                node.impact(ref_phase, phase, reference)?
            ));
            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
        }
        Ok(protein)
    }
}

#[derive(Debug, Clone)]
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

#[derive(Debug, Clone)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) struct Node {
    node_type: NodeType,
    vaf: HashMap<String, f32>,
    probs: EventProbs,
    pos: i64,
    index: u32,
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

    pub(crate) fn reference_amino_acid(&self, phase: u8, reference: &[u8]) -> Result<AminoAcid> {
        let start_pos = self.pos as usize - ((self.pos - phase as i64) % 3) as usize;
        let ref_codon_bases = reference[start_pos..start_pos + 3].to_vec();
        AminoAcid::from_codon(
            transcription::transcribe_dna_to_rna(ref_codon_bases.as_ref())?.as_ref(),
        )
    }

    pub(crate) fn variant_amino_acid(&self, phase: u8, reference: &[u8]) -> Result<AminoAcid> {
        match &self.node_type {
            NodeType::Var(alt_allele) => {
                let start_pos = self.pos as usize - ((self.pos - phase as i64) % 3) as usize;
                let position_in_codon = (self.pos - phase as i64) % 3;
                let needed_bases = if alt_allele.is_empty() { 4 } else { 3 };
                let ref_codon_bases = reference[start_pos..start_pos + needed_bases].to_vec();
                let alt_codon_bases = [
                    &ref_codon_bases[..position_in_codon as usize],
                    alt_allele.as_bytes(),
                    &ref_codon_bases[position_in_codon as usize + 1..],
                ]
                .concat();
                println!("Ref: {:?}", String::from_utf8_lossy(&ref_codon_bases));
                println!("Alt: {:?}", String::from_utf8_lossy(&alt_codon_bases));
                AminoAcid::from_codon(
                    transcription::transcribe_dna_to_rna(&alt_codon_bases[..3])?.as_ref(), // TODO: How do we want to consider insertions greater than 2 that will span multiple codons?
                )
            }
            _ => {
                unreachable!("Reference node type has no variant amino acid")
            }
        }
    }

    pub(crate) fn impact(&self, ref_phase: u8, phase: u8, reference: &[u8]) -> Result<Impact> {
        match &self.node_type {
            NodeType::Var(_) => {
                let ref_amino_acid = self.reference_amino_acid(ref_phase, reference)?;
                let alt_amino_acid = self.variant_amino_acid(phase, reference)?;
                match (
                    ref_amino_acid == alt_amino_acid,
                    ref_amino_acid,
                    alt_amino_acid,
                ) {
                    (true, _, _) => Ok(Impact::Low),
                    (false, _, AminoAcid::Stop) => Ok(Impact::High),
                    (false, AminoAcid::Stop, _) => Ok(Impact::High),
                    (false, AminoAcid::Methionine, _) => Ok(Impact::High), // TODO: Check if this is always automatic start lost or can Met occur anywhere in the protein?
                    (false, _, _) => Ok(Impact::Modifier),
                }
            }
            NodeType::Ref(_) => Ok(Impact::None),
        }
    }
}

#[derive(Debug, Clone)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
struct EventProbs(HashMap<String, f32>);

impl EventProbs {
    fn from_record(record: &Record, tags: &Vec<String>) -> Self {
        let mut probs = HashMap::new();
        for tag in tags {
            let prob = record.info(tag.as_bytes()).float().unwrap().unwrap()[0];
            probs.insert(tag.to_string(), prob);
        }
        EventProbs(probs)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Edge {
    pub(crate) supporting_reads: HashMap<String, u32>,
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

// test node distance
#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::{Directed, Graph};
    use rust_htslib::bcf::{Read, Reader};
    use std::fs;

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
    fn test_event_probs_from_record() {
        let mut reader = Reader::from_path("tests/resources/calls.bcf").unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let tags = vec![
            "PROB_ABSENT".to_string(),
            "PROB_PRESENT".to_string(),
            "PROB_ARTIFACT".to_string(),
        ];
        let event_probs = EventProbs::from_record(&record, &tags);
        assert_eq!(event_probs.0.len(), 3);
        assert_eq!(event_probs.0.get("PROB_ABSENT").unwrap(), &0.036097374);
        assert_eq!(event_probs.0.get("PROB_PRESENT").unwrap(), &20.82111);
        assert_eq!(event_probs.0.get("PROB_ARTIFACT").unwrap(), &f32::INFINITY);
    }

    #[test]
    fn test_build_graph() {
        let calls_file = PathBuf::from("tests/resources/test_calls.vcf");
        let observations_file = PathBuf::from("tests/resources/test_observations.vcf");
        let observations = vec![ObservationFile {
            path: observations_file,
            sample: "sample".to_string(),
        }];
        let variant_graph = VariantGraph::build(&calls_file, &observations, "OX512233.1", 60, 85);
        assert!(variant_graph.is_ok());
    }

    #[test]
    fn test_graph_is_empty() {
        let calls_file = PathBuf::from("tests/resources/test_calls.vcf");
        let observations_file = PathBuf::from("tests/resources/test_observations.vcf");
        let observations = vec![ObservationFile {
            path: observations_file,
            sample: "sample".to_string(),
        }];
        let variant_graph =
            VariantGraph::build(&calls_file, &observations, "not actually in file", 60, 85);
        assert!(variant_graph.unwrap().is_empty());
    }

    #[test]
    fn test_graph_paths() {
        let calls_file = PathBuf::from("tests/resources/test_calls.vcf");
        let observations_file = PathBuf::from("tests/resources/test_observations.vcf");
        let observations = vec![ObservationFile {
            path: observations_file,
            sample: "sample".to_string(),
        }];
        let mut variant_graph =
            VariantGraph::build(&calls_file, &observations, "OX512233.1", 60, 85).unwrap();
        variant_graph.graph.add_edge(
            NodeIndex::new(0),
            NodeIndex::new(2),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        variant_graph.graph.add_edge(
            NodeIndex::new(2),
            NodeIndex::new(5),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        variant_graph.graph.add_edge(
            NodeIndex::new(5),
            NodeIndex::new(6),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let paths = variant_graph.paths();
        let expected_paths = vec![
            HaplotypePath(vec![
                NodeIndex::new(0),
                NodeIndex::new(2),
                NodeIndex::new(5),
                NodeIndex::new(6),
            ]),
            HaplotypePath(vec![
                NodeIndex::new(0),
                NodeIndex::new(2),
                NodeIndex::new(5),
                NodeIndex::new(7),
            ]),
            HaplotypePath(vec![
                NodeIndex::new(1),
                NodeIndex::new(3),
                NodeIndex::new(5),
                NodeIndex::new(6),
            ]),
            HaplotypePath(vec![
                NodeIndex::new(1),
                NodeIndex::new(3),
                NodeIndex::new(5),
                NodeIndex::new(7),
            ]),
        ];

        assert_eq!(paths, expected_paths);
    }

    #[test]
    fn test_variant_amino_acid_with_different_phases() {
        let node = Node::new(NodeType::Var("A".to_string()), 2);
        let reference = b"ATGCGCGTA";
        let ile = node.variant_amino_acid(0, reference).unwrap();
        assert_eq!(ile, AminoAcid::Isoleucine);
        let tyr = node.variant_amino_acid(1, reference).unwrap();
        assert_eq!(tyr, AminoAcid::Tyrosine);
        let thr = node.variant_amino_acid(2, reference).unwrap();
        assert_eq!(thr, AminoAcid::Threonine);
    }

    #[test]
    fn test_variant_amino_acid_with_deletion() {
        let node = Node::new(NodeType::Var("".to_string()), 2);
        let reference = b"ATGCCGT";
        let ile = node.variant_amino_acid(0, reference).unwrap();
        assert_eq!(ile, AminoAcid::Isoleucine);
        let ser = node.variant_amino_acid(1, reference).unwrap();
        assert_eq!(ser, AminoAcid::Serine);
        let pro = node.variant_amino_acid(2, reference).unwrap();
        assert_eq!(pro, AminoAcid::Proline);
    }

    #[test]
    fn test_variant_amino_acid_with_one_inserted_base() {
        let node = Node::new(NodeType::Var("GG".to_string()), 1);
        let reference = b"AGCTCT";
        let arg = node.variant_amino_acid(0, reference).unwrap();
        assert_eq!(arg, AminoAcid::Arginine);
        let gly = node.variant_amino_acid(1, reference).unwrap();
        assert_eq!(gly, AminoAcid::Glycine);
        let node_2 = Node::new(NodeType::Var("GG".to_string()), 3);
        let arg = node_2.variant_amino_acid(2, reference).unwrap();
        assert_eq!(arg, AminoAcid::Arginine);
    }

    #[test]
    fn test_variant_amino_acid_with_two_inserted_bases() {
        let node = Node::new(NodeType::Var("TCC".to_string()), 4);
        let reference = b"ATCATCATC";
        let ile = node.variant_amino_acid(0, reference).unwrap();
        assert_eq!(ile, AminoAcid::Isoleucine);
        let ser = node.variant_amino_acid(1, reference).unwrap();
        assert_eq!(ser, AminoAcid::Serine);
        let his = node.variant_amino_acid(2, reference).unwrap();
        assert_eq!(his, AminoAcid::Histidine);
    }

    #[test]
    fn test_reference_amino_acid_with_different_phases() {
        let node = Node::new(NodeType::Ref("".to_string()), 2);
        let reference = b"ATGCGCGTA";
        let met = node.reference_amino_acid(0, reference).unwrap();
        assert_eq!(met, AminoAcid::Methionine);
        let cys = node.reference_amino_acid(1, reference).unwrap();
        assert_eq!(cys, AminoAcid::Cysteine);
        let ala = node.reference_amino_acid(2, reference).unwrap();
        assert_eq!(ala, AminoAcid::Alanine);
    }
    #[test]
    fn impact_returns_none_for_ref_node_type() {
        let node = Node::new(NodeType::Ref("".to_string()), 0);
        let impact = node.impact(0, 0, &[]).unwrap();
        assert_eq!(impact, Impact::None);
    }

    #[test]
    fn impact_identifies_low_for_identical_ref_and_alt_amino_acids() {
        let node = Node::new(NodeType::Var("C".to_string()), 2);
        let reference = b"ATA";
        let impact = node.impact(0, 0, reference).unwrap();
        assert_eq!(impact, Impact::Low);
    }

    #[test]
    fn impact_identifies_modifier_for_different_ref_and_alt_amino_acids() {
        let node = Node::new(NodeType::Var("G".to_string()), 3);
        let reference = b"ATTTG";
        let impact = node.impact(2, 2, reference).unwrap();
        assert_eq!(impact, Impact::Modifier);
    }

    #[test]
    fn impact_identifies_high_for_early_stop() {
        let node = Node::new(NodeType::Var("A".to_string()), 6);
        let reference = b"CATATAC";
        let impact = node.impact(1, 1, reference).unwrap();
        assert_eq!(impact, Impact::High);
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

    fn setup_variant_graph_with_nodes() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node1 = graph.add_node(Node::new(NodeType::Var("A".to_string()), 1));
        let node2 = graph.add_node(Node::new(NodeType::Ref("".to_string()), 2));
        let node3 = graph.add_node(Node::new(NodeType::Var("T".to_string()), 3));
        let node4 = graph.add_node(Node::new(NodeType::Var("".to_string()), 4));
        let node5 = graph.add_node(Node::new(NodeType::Var("A".to_string()), 8));
        let node6 = graph.add_node(Node::new(NodeType::Var("TT".to_string()), 9));
        VariantGraph {
            graph,
            start: 0,
            end: 2,
            target: "test".to_string(),
        }
    }

    #[test]
    fn impact_calculates_none_for_empty_path() {
        let graph = setup_variant_graph_with_nodes();
        let path = HaplotypePath(vec![]);
        let impact = path.impact(&graph, 0, b"TTG").unwrap();
        assert_eq!(impact, Impact::None);
    }

    #[test]
    fn impact_calculates_correctly_for_single_variant_node() {
        let mut graph = setup_variant_graph_with_nodes();
        let node_index = graph.graph.node_indices().next().unwrap();
        let path = HaplotypePath(vec![node_index]);
        let impact = path.impact(&graph, 0, b"TTC").unwrap();
        assert_eq!(impact, Impact::Modifier);
    }

    #[test]
    fn impact_accumulates_over_multiple_nodes() {
        let mut graph = setup_variant_graph_with_nodes();
        let node_indices = graph.graph.node_indices().take(3).collect::<Vec<_>>();
        let path = HaplotypePath(node_indices.clone());
        let impact = path.impact(&graph, 0, b"TTCAAA").unwrap();
        assert_eq!(impact, Impact::High);
    }

    #[test]
    fn impact_handles_phase_shift_correctly() {
        let mut graph = setup_variant_graph_with_nodes();
        let node_index = graph.graph.node_indices().next().unwrap();
        let path = HaplotypePath(vec![node_index]);
        let impact = path.impact(&graph, 1, b"ATGA").unwrap();
        assert_eq!(impact, Impact::High);
    }

    #[test]
    fn impact_handles_phase_shift_caused_by_frameshift() {
        let mut graph = setup_variant_graph_with_nodes();
        let node_indices = graph
            .graph
            .node_indices()
            .skip(3)
            .take(2)
            .collect::<Vec<_>>();
        let path = HaplotypePath(node_indices.clone());
        let impact = path.impact(&graph, 0, b"GGGAAATTTAAA").unwrap();
        assert_eq!(impact, Impact::High);
    }

    #[test]
    fn impact_handles_phase_shift_caused_by_frameshift_2() {
        let mut graph = setup_variant_graph_with_nodes();
        let node_indices = graph.graph.node_indices().skip(3).take(3).collect_vec();
        let path = HaplotypePath(vec![node_indices[0], node_indices[2]]);
        let impact = path.impact(&graph, 0, b"GGGAAATTTAAC").unwrap();
        assert_eq!(impact, Impact::Modifier);
    }

    fn setup_variant_graph_with_nodes_2() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node1 = graph.add_node(Node::new(NodeType::Var("AA".to_string()), 1));
        let node2 = graph.add_node(Node::new(NodeType::Var("".to_string()), 7));
        VariantGraph {
            graph,
            start: 0,
            end: 2,
            target: "test".to_string(),
        }
    }

    #[test]
    fn impact_handles_phase_shift_caused_by_frameshift_3() {
        let mut graph = setup_variant_graph_with_nodes_2();
        let node_indices = graph.graph.node_indices().collect_vec();
        let path = HaplotypePath(node_indices.clone());
        let impact = path.impact(&graph, 0, b"ATGAAATGGAT").unwrap();
        assert_eq!(impact, Impact::High);
    }

    #[test]
    fn impact_handles_phase_shift_caused_by_big_frameshift() {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node1 = graph.add_node(Node::new(NodeType::Var("TTTTT".to_string()), 1));
        let node2 = graph.add_node(Node::new(NodeType::Var("GG".to_string()), 7));
        let graph = VariantGraph {
            graph,
            start: 0,
            end: 20,
            target: "test".to_string(),
        };
        let node_indices = graph.graph.node_indices().collect_vec();
        let path = HaplotypePath(node_indices.clone());
        let impact = path.impact(&graph, 0, b"TGTTTTTAATTT").unwrap();
        println!("{}", path.display(&graph, 0, b"TGTTTTTAATTT").unwrap());
        assert_eq!(impact, Impact::High);
    }

    #[test]
    fn is_variant_returns_true_for_variant_node() {
        let node_type = NodeType::Var("A".to_string());
        assert!(node_type.is_variant());
    }

    #[test]
    fn is_variant_returns_false_for_reference_node() {
        let node_type = NodeType::Ref("A".to_string());
        assert!(!node_type.is_variant());
    }

    #[test]
    fn display_returns_correct_protein_string_for_single_node() {
        let graph = setup_variant_graph_with_nodes();
        let path = HaplotypePath(vec![NodeIndex::new(0)]);
        let result = path.display(&graph, 0, b"ATG").unwrap();
        assert_eq!(result, "Met -> Lys (High)\n");
    }

    #[test]
    fn display_returns_correct_protein_string_for_multiple_nodes() {
        let graph = setup_variant_graph_with_nodes();
        let path = HaplotypePath(vec![NodeIndex::new(0), NodeIndex::new(2)]);
        let result = path.display(&graph, 0, b"ATGCGT").unwrap();
        assert_eq!(result, "Met -> Lys (High)\nArg -> Cys (Modifier)\n");
    }

    #[test]
    fn to_dot_generates_correct_dot_representation() {
        let graph = setup_variant_graph_with_nodes();
        let dot_output = graph.to_dot();
        assert!(dot_output.contains("digraph {"));
        assert!(dot_output.contains("}"));
    }

    #[test]
    fn to_file_writes_dot_file_correctly() {
        let graph = setup_variant_graph_with_nodes();
        let temp_dir = tempfile::tempdir().unwrap();
        let file_path = temp_dir.path().join("graph.dot");
        graph.to_file(&temp_dir.path()).unwrap();
        let written_content = fs::read_to_string(file_path).unwrap();
        assert!(written_content.contains("digraph {"));
        assert!(written_content.contains("}"));
    }
}
