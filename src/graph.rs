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

pub(crate) struct VariantGraph(pub(crate) Graph<Node, Edge, Directed>);

impl VariantGraph {
    pub(crate) fn build(
        calls_file: &PathBuf,
        observations_file: &PathBuf,
        output_path: &Path,
    ) -> Result<()> {
        let mut calls_reader = Reader::from_path(calls_file)?;
        let mut observations_reader = Reader::from_path(observations_file)?;

        let event_names = extract_event_names(calls_file);

        let tags = event_names
            .iter()
            .map(|event| format!("PROB_{event}"))
            .collect();

        let mut supporting_reads = HashMap::new();

        let mut variant_graph = Graph::<Node, Edge, Directed>::new();
        let mut batch = 0;

        for (calls_record, observations_record) in
            calls_reader.records().zip(observations_reader.records())
        {
            let mut calls_record = calls_record?;
            let mut observations_record = observations_record?;

            let _variants = collect_variants(&mut calls_record, false, None)?;
            let observations = read_observations(&mut observations_record)?;
            let observations = observations.pileup.read_observations();
            let fragment_ids: HashSet<_> = observations.iter().map(|o| o.fragment_id).collect();

            if !fragment_ids
                .iter()
                .any(|id| supporting_reads.contains_key(id))
                && !supporting_reads.is_empty()
            {
                let mut batch_graph = VariantGraph(variant_graph.clone());
                batch_graph.create_edges(&supporting_reads)?;
                batch_graph.to_file(output_path, batch)?;
                variant_graph.clear();
                supporting_reads.clear();
                batch += 1;
            }

            let alleles = calls_record.alleles();
            let ref_allele = String::from_utf8(alleles[0].to_vec())?;
            let alt_allele = String::from_utf8(alleles[1].to_vec())?;

            let var_node = Node::from_records(
                &calls_record,
                &observations_record,
                &tags,
                NodeType::Var(alt_allele),
            );
            let var_node_index = variant_graph.add_node(var_node);

            let ref_node = Node::from_records(
                &calls_record,
                &observations_record,
                &tags,
                NodeType::Ref(ref_allele),
            );
            let ref_node_index = variant_graph.add_node(ref_node);

            for observation in observations {
                let evidence = BayesFactor::new(observation.prob_alt, observation.prob_ref)
                    .evidence_kass_raftery();
                match evidence {
                    KassRaftery::Strong | KassRaftery::VeryStrong => {
                        // Read supports variant
                        let entry = supporting_reads
                            .entry(observation.fragment_id)
                            .or_insert(Vec::new());
                        entry.push(var_node_index);
                    }
                    KassRaftery::None | KassRaftery::Barely | KassRaftery::Positive => {
                        // Read supports reference
                        let entry = supporting_reads
                            .entry(observation.fragment_id)
                            .or_insert(Vec::new());
                        entry.push(ref_node_index);
                    }
                }
            }
        }

        let mut variant_graph = VariantGraph(variant_graph.clone());
        variant_graph.create_edges(&supporting_reads)?;
        variant_graph.to_file(output_path, batch)?;

        Ok(())
    }

    pub(crate) fn create_edges(
        &mut self,
        supporting_reads: &HashMap<Option<u64>, Vec<NodeIndex>>,
    ) -> Result<()> {
        for nodes in supporting_reads.values() {
            for node_tuple in nodes
                .iter()
                .sorted()
                .dedup()
                .combinations(2)
                .filter(|v| node_distance(&v[0].index(), &v[1].index()) <= 1)
            {
                let edge = self.0.find_edge(*node_tuple[0], *node_tuple[1]);
                if let Some(edge) = edge {
                    let edge = self.0.edge_weight_mut(edge).unwrap();
                    edge.supporting_reads += 1;
                } else {
                    let edge = Edge {
                        supporting_reads: 1,
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

    pub(crate) fn to_file(&self, path: &Path, batch_number: u32) -> Result<()> {
        let path = path.join(format!("graph_{}.dot", batch_number));
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
    vaf: f32,
    probs: EventProbs,
}

impl Node {
    pub(crate) fn from_records(
        calls_record: &Record,
        _observations_record: &Record,
        tags: &Vec<String>,
        node_type: NodeType,
    ) -> Self {
        let vaf = calls_record.format(b"AF").float().unwrap()[0][0];
        Node {
            node_type,
            vaf: vaf.to_owned(),
            probs: EventProbs::from_record(calls_record, tags),
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
    pub(crate) supporting_reads: u32,
}

pub(crate) fn node_distance(node1: &usize, node2: &usize) -> usize {
    let distance = (*node1 as isize - *node2 as isize).unsigned_abs();
    if node1 % 2 == 0 {
        distance / 2
    } else {
        (distance + 1) / 2
    }
}

// test node distance
#[cfg(test)]
mod tests {
    use std::fs;
    use super::*;
    use petgraph::{Directed, Graph};
    use rust_htslib::bcf::{Read, Reader};

    #[test]
    fn test_node_distance() {
        let mut graph = Graph::<u32, Edge, Directed>::new();
        let weight = 1;
        let node0 = graph.add_node(weight.clone());
        let node1 = graph.add_node(weight.clone());
        let node2 = graph.add_node(weight.clone());
        let node3 = graph.add_node(weight.clone());
        let node4 = graph.add_node(weight.clone());

        let distance = node_distance(&node0.index(), &node2.index());
        assert_eq!(distance, 1);
        let distance_2 = node_distance(&node0.index(), &node1.index());
        assert_eq!(distance_2, 0);
        let distance_3 = node_distance(&node1.index(), &node4.index());
        assert_eq!(distance_3, 2);
        let distance_4 = node_distance(&node3.index(), &node4.index());
        assert_eq!(distance_4, 1);
        let distance_5 = node_distance(&node1.index(), &node3.index());
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
        let tmp = PathBuf::from("/tmp/");
        VariantGraph::build(&calls_file, &observations_file, &tmp).unwrap();
        let output = tmp.join("graph_0.dot");
        assert!(output.exists());
        fs::remove_file(output).unwrap();
    }
}
