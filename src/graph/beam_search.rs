use crate::graph::node::Node;
use crate::graph::paths::HaplotypePath;
use crate::graph::score::HaplotypeMetric;
use crate::graph::VariantGraph;
use itertools::Itertools;
use petgraph::graph::NodeIndex;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::f32;

#[derive(Debug, Clone)]
pub struct BeamSearchConfig {
    pub beam_width: usize,
    pub haplotype_metric: HaplotypeMetric,
}

impl Default for BeamSearchConfig {
    fn default() -> Self {
        BeamSearchConfig {
            beam_width: 100,
            haplotype_metric: HaplotypeMetric::GeometricMean,
        }
    }
}

#[derive(Clone, Debug)]
struct ScoredPath {
    path: Vec<NodeIndex>,
    sample_likelihoods: HashMap<String, f32>,
}

impl ScoredPath {
    fn max_likelihood(&self) -> f32 {
        self.sample_likelihoods
            .values()
            .cloned()
            .fold(f32::NEG_INFINITY, f32::max)
    }
}

impl PartialEq for ScoredPath {
    fn eq(&self, other: &Self) -> bool {
        self.max_likelihood() == other.max_likelihood()
    }
}

impl Eq for ScoredPath {}

impl PartialOrd for ScoredPath {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ScoredPath {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

impl ScoredPath {
    fn from_node(node: NodeIndex, graph: &VariantGraph, metric: HaplotypeMetric) -> Self {
        let nodes = vec![graph.graph.node_weight(node).unwrap().clone()];
        let sample_likelihoods = metric.calculate(&nodes);

        ScoredPath {
            path: vec![node],
            sample_likelihoods,
        }
    }

    fn extend(&self, node: NodeIndex, graph: &VariantGraph, metric: HaplotypeMetric) -> Self {
        let mut new_path = self.path.clone();
        new_path.push(node);

        let nodes: Vec<Node> = new_path
            .iter()
            .map(|&idx| graph.graph.node_weight(idx).unwrap().clone())
            .collect();

        let sample_likelihoods = metric.calculate(&nodes);

        ScoredPath {
            path: new_path,
            sample_likelihoods,
        }
    }
}

impl VariantGraph {
    pub fn beam_search_paths(&self, config: BeamSearchConfig) -> Vec<HaplotypePath> {
        let mut complete_paths = BinaryHeap::new();
        let start_nodes: Vec<NodeIndex> = self.graph.node_indices().take(2).collect();

        if start_nodes.is_empty() {
            return vec![];
        }

        let mut current_beam = BinaryHeap::new();
        for start_node in start_nodes {
            current_beam.push(ScoredPath::from_node(
                start_node,
                self,
                config.haplotype_metric,
            ));
        }

        let mut visited: HashSet<(NodeIndex, Vec<NodeIndex>)> = HashSet::new();

        while !current_beam.is_empty() {
            let mut next_beam = BinaryHeap::new();

            let mut beam_items: Vec<_> = current_beam.into_iter().collect();
            beam_items.sort_by(|a, b| b.cmp(a));
            beam_items.truncate(config.beam_width);

            for scored_path in beam_items {
                let current_node = *scored_path.path.last().unwrap();

                let mut variant_nodes: Vec<NodeIndex> = scored_path
                    .path
                    .iter()
                    .filter(|&&idx| self.graph.node_weight(idx).unwrap().node_type.is_variant())
                    .cloned()
                    .collect();
                variant_nodes.sort();

                let state = (current_node, variant_nodes);
                if visited.contains(&state) {
                    continue;
                }
                visited.insert(state);

                let neighbors: Vec<NodeIndex> = self
                    .graph
                    .neighbors(current_node)
                    .filter(|&n| n.index() > current_node.index())
                    .filter(|&n| !scored_path.path.contains(&n))
                    .collect();

                if neighbors.is_empty() {
                    let nodes: Vec<&Node> = scored_path
                        .path
                        .iter()
                        .map(|&idx| self.graph.node_weight(idx).unwrap())
                        .collect();

                    if nodes.iter().all(|n| n.pos != -1) {
                        complete_paths.push(scored_path);
                    }
                } else {
                    for neighbor in neighbors {
                        let extended = scored_path.extend(neighbor, self, config.haplotype_metric);
                        next_beam.push(extended);
                    }
                }
            }

            current_beam = next_beam;
        }

        let mut paths: Vec<(HaplotypePath, f32)> = complete_paths
            .into_iter()
            .map(|sp| (HaplotypePath(sp.path.clone()), sp.max_likelihood()))
            .collect();

        paths.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
        paths.into_iter().map(|(path, _)| path).unique().collect()
    }

    pub fn beam_search_reverse_paths(&self, config: BeamSearchConfig) -> Vec<HaplotypePath> {
        self.beam_search_paths(config)
            .iter()
            .map(|path| HaplotypePath(path.0.iter().rev().cloned().collect()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::node::NodeType;
    use crate::graph::{Edge, EventProbs};
    use petgraph::{Directed, Graph};
    use std::collections::HashMap;

    fn create_test_graph_small() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();

        let ref_node = graph.add_node(Node {
            node_type: NodeType::Reference,
            reference_allele: "A".to_string(),
            alternative_allele: "A".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 0.8)]),
            probs: EventProbs(HashMap::new()),
            pos: 0,
            index: 0,
        });

        let var1 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "T".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 0.6)]),
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 1,
        });

        let var2 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "G".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 0.5)]),
            probs: EventProbs(HashMap::new()),
            pos: 2,
            index: 2,
        });

        graph.add_edge(
            ref_node,
            var1,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graph.add_edge(
            var1,
            var2,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );

        VariantGraph {
            graph,
            start: 0,
            end: 2,
            target: "test".to_string(),
        }
    }

    fn create_test_graph_large() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();

        let mut all_nodes = vec![];
        for i in 0..30 {
            let ref_node = graph.add_node(Node {
                node_type: NodeType::Reference,
                reference_allele: "A".to_string(),
                alternative_allele: "A".to_string(),
                vaf: HashMap::from([("sample1".to_string(), 0.5)]),
                probs: EventProbs(HashMap::new()),
                pos: i as i64,
                index: i * 2,
            });

            let var_node = graph.add_node(Node {
                node_type: NodeType::Variant,
                reference_allele: "A".to_string(),
                alternative_allele: "T".to_string(),
                vaf: HashMap::from([("sample1".to_string(), if i < 15 { 0.7 } else { 0.3 })]),
                probs: EventProbs(HashMap::new()),
                pos: i as i64,
                index: i * 2 + 1,
            });

            all_nodes.push((ref_node, var_node));
        }

        for i in 0..all_nodes.len() - 1 {
            let (curr_ref, curr_var) = all_nodes[i];
            let (next_ref, next_var) = all_nodes[i + 1];
            graph.add_edge(
                curr_ref,
                next_ref,
                Edge {
                    supporting_reads: HashMap::new(),
                },
            );
            graph.add_edge(
                curr_ref,
                next_var,
                Edge {
                    supporting_reads: HashMap::new(),
                },
            );
            graph.add_edge(
                curr_var,
                next_ref,
                Edge {
                    supporting_reads: HashMap::new(),
                },
            );
            graph.add_edge(
                curr_var,
                next_var,
                Edge {
                    supporting_reads: HashMap::new(),
                },
            );
        }

        VariantGraph {
            graph,
            start: 0,
            end: 29,
            target: "test".to_string(),
        }
    }

    #[test]
    fn test_beam_search_small_graph() {
        let graph = create_test_graph_small();
        let config = BeamSearchConfig {
            beam_width: 10,
            haplotype_metric: HaplotypeMetric::GeometricMean,
        };

        let paths = graph.beam_search_paths(config);
        assert!(!paths.is_empty());
    }

    #[test]
    fn test_beam_search_large_graph() {
        let graph = create_test_graph_large();
        let config = BeamSearchConfig {
            beam_width: 100,
            haplotype_metric: HaplotypeMetric::GeometricMean,
        };

        let paths = graph.beam_search_paths(config);
        assert!(!paths.is_empty());
    }

    #[test]
    fn test_beam_search_respects_beam_width() {
        let graph = create_test_graph_large();

        let config_narrow = BeamSearchConfig {
            beam_width: 5,
            haplotype_metric: HaplotypeMetric::GeometricMean,
        };
        let paths_narrow = graph.beam_search_paths(config_narrow);

        let config_wide = BeamSearchConfig {
            beam_width: 100,
            haplotype_metric: HaplotypeMetric::GeometricMean,
        };
        let paths_wide = graph.beam_search_paths(config_wide);
        assert!(paths_narrow.len() < paths_wide.len());
    }

    #[test]
    fn test_beam_search_different_metrics() {
        let graph = create_test_graph_small();

        for metric in [
            HaplotypeMetric::Product,
            HaplotypeMetric::GeometricMean,
            HaplotypeMetric::Minimum,
        ] {
            let config = BeamSearchConfig {
                beam_width: 10,
                haplotype_metric: metric,
            };

            let paths = graph.beam_search_paths(config);
            assert!(!paths.is_empty());
        }
    }

    #[test]
    fn test_beam_search_reverse_paths() {
        let graph = create_test_graph_small();
        let config = BeamSearchConfig::default();

        let forward = graph.beam_search_paths(config.clone());
        let reverse = graph.beam_search_reverse_paths(config);

        assert_eq!(forward.len(), reverse.len());
    }

    #[test]
    fn test_beam_search_keeps_highest_likelihood_paths() {
        let mut graph = Graph::<Node, Edge, Directed>::new();

        let ref_node = graph.add_node(Node {
            node_type: NodeType::Reference,
            reference_allele: "A".to_string(),
            alternative_allele: "A".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 1.0)]),
            probs: EventProbs(HashMap::new()),
            pos: 0,
            index: 0,
        });
        let high_vaf_var = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "T".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 0.9)]),
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 1,
        });
        let low_vaf_var = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "G".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 0.1)]),
            probs: EventProbs(HashMap::new()),
            pos: 1,
            index: 2,
        });
        let final_node = graph.add_node(Node {
            node_type: NodeType::Reference,
            reference_allele: "C".to_string(),
            alternative_allele: "C".to_string(),
            vaf: HashMap::from([("sample1".to_string(), 1.0)]),
            probs: EventProbs(HashMap::new()),
            pos: 2,
            index: 3,
        });

        graph.add_edge(
            ref_node,
            high_vaf_var,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graph.add_edge(
            ref_node,
            low_vaf_var,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graph.add_edge(
            high_vaf_var,
            final_node,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graph.add_edge(
            low_vaf_var,
            final_node,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );

        let variant_graph = VariantGraph {
            graph,
            start: 0,
            end: 2,
            target: "test".to_string(),
        };

        let config = BeamSearchConfig {
            beam_width: 1,
            haplotype_metric: HaplotypeMetric::Product,
        };

        let paths = variant_graph.beam_search_paths(config);

        assert_eq!(paths.len(), 1);

        let path = &paths[0].0;
        assert_eq!(path.len(), 3);
        assert_eq!(path[0], NodeIndex::new(0));
        assert_eq!(path[1], NodeIndex::new(1));
        assert_eq!(path[2], NodeIndex::new(3));
    }
}
