use crate::cli::ObservationFile;
use crate::utils::bcf::extract_event_names;
use anyhow::Result;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;
use petgraph::dot::{Config, Dot};
use petgraph::graph::NodeIndex;
use petgraph::{Directed, Graph};
use rust_htslib::bcf::{Read, Reader, Record};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use varlociraptor::calling::variants::preprocessing::read_observations;
use varlociraptor::utils::collect_variants::collect_variants;
use varlociraptor::variants::evidence::observations::read_observation::ProcessedReadObservation;

pub(crate) struct VariantGraph(pub(crate) Graph<Node, Edge, Directed>);

impl VariantGraph {
    pub(crate) fn build(
        calls_file: &PathBuf,
        observation_files: &[ObservationFile],
    ) -> Result<VariantGraph> {
        let mut calls_reader = Reader::from_path(calls_file)?;
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
            let fragment_ids: HashSet<_> = observations
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

        let mut variant_graph = VariantGraph(variant_graph.clone());
        variant_graph.create_edges(&supporting_reads)?;

        Ok(variant_graph)
    }

    pub(crate) fn create_edges(
        &mut self,
        supporting_reads: &HashMap<(String, Option<u64>), Vec<NodeIndex>>,
    ) -> Result<()> {
        for ((sample, _), nodes) in supporting_reads {
            let weights = self.0.clone();
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
                let edge = self.0.find_edge(*node_tuple[0], *node_tuple[1]);
                if let Some(edge) = edge {
                    let edge = self.0.edge_weight_mut(edge).unwrap();
                    edge.supporting_reads
                        .entry(sample.to_string())
                        .and_modify(|v| *v += 1)
                        .or_insert(1);
                } else {
                    let edge = Edge {
                        supporting_reads: HashMap::from([(sample.to_string(), 1)]),
                    };
                    self.0.add_edge(*node_tuple[0], *node_tuple[1], edge);
                }
            }
        }
        Ok(())
    }

    pub(crate) fn to_dot(&self) -> String {
        format!(
            "digraph {{ {:?} }}",
            Dot::with_config(&self.0, &[Config::GraphContentOnly])
        )
    }

    pub(crate) fn to_file(&self, path: &Path) -> Result<()> {
        let path = path.join("graph.dot");
        let mut file = std::fs::File::create(path)?;
        file.write_all(self.to_dot().as_bytes())?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) enum NodeType {
    Var(String),
    Ref(String),
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
        let tmp = PathBuf::from("/tmp/");
        VariantGraph::build(&calls_file, &observations, &tmp).unwrap();
        let output = tmp.join("graph_0.dot");
        assert!(output.exists());
        let contents = fs::read_to_string(&output).unwrap();
        assert_eq!(contents.lines().count(), 12);
        fs::remove_file(&output).unwrap();
    }
}
