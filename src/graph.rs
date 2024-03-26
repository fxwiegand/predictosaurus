use std::collections::HashMap;
use rust_htslib::bcf::Record;

#[derive(Debug, Clone)]
pub(crate) enum NodeType {
    Var(String),
    Ref(String),
}

#[derive(Debug, Clone)]
pub(crate) struct Node {
    node_type: NodeType,
    vaf: f32,
    probs: EventProbs,
}

impl Node {
    pub(crate) fn from_observations(
        record: &Record,
        tags: &Vec<String>,
        node_type: NodeType,
    ) -> Self {
        let vaf = record.info(b"AF").float().expect("No AF tag found in record").unwrap()[0];
        Node {
            node_type,
            vaf,
            probs: EventProbs::from_record(&record, &tags),
        }
    }
}

#[derive(Debug, Clone)]
struct  EventProbs(HashMap<String, f32>);

impl EventProbs {
    fn from_record(record: &Record, tags: &Vec<String>) -> Self {
        unimplemented!()
    }
}

#[derive(Debug)]
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
    use super::*;
    use petgraph::{Directed, Graph};

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
}
