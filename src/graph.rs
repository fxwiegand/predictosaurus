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
    pub(crate) fn from_records(
        calls_record: &Record,
        observations_record: &Record,
        tags: &Vec<String>,
        node_type: NodeType,
    ) -> Self {
        let vaf = calls_record.format(b"AF").float().unwrap()[0][0];
        Node {
            node_type,
            vaf: vaf.to_owned(),
            probs: EventProbs::from_record(&calls_record, &tags),
        }
    }
}

#[derive(Debug, Clone)]
struct  EventProbs(HashMap<String, f32>);

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
        let tags = vec!["PROB_ABSENT".to_string(), "PROB_PRESENT".to_string(), "PROB_ARTIFACT".to_string()];
        let event_probs = EventProbs::from_record(&record, &tags);
        assert_eq!(event_probs.0.len(), 3);
        assert_eq!(event_probs.0.get("PROB_ABSENT").unwrap(), &0.036097374);
        assert_eq!(event_probs.0.get("PROB_PRESENT").unwrap(), &20.82111);
        assert_eq!(event_probs.0.get("PROB_ARTIFACT").unwrap(), &f32::INFINITY);
    }
}
