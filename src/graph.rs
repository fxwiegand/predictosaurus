use varlociraptor::calling::variants::preprocessing::Observations;
use varlociraptor::utils::collect_variants::VariantInfo;

#[derive(Debug, Clone)]
pub(crate) enum NodeType {
    Var,
    Ref,
}

#[derive(Debug, Clone)]
pub(crate) struct Node {
    node_type: NodeType,
    variant_info: VariantInfo,
}

impl Node {
    pub(crate) fn from_observations(
        variant_info: &VariantInfo,
        _observations: &Observations,
        node_type: NodeType,
    ) -> Self {
        Node {
            node_type,
            variant_info: variant_info.to_owned(),
        }
    }
}

#[derive(Debug)]
pub(crate) struct Edge {
    pub(crate) supporting_reads: u32,
}

pub(crate) fn node_distance(node1: &usize, node2: &usize) -> usize {
    let distance = (*node1 as isize - *node2 as isize).unsigned_abs() / 2;
    if node1 % 2 == 0 {
        distance
    } else {
        distance + 1
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
        let node1 = graph.add_node(weight.clone());
        let node2 = graph.add_node(weight.clone());
        let node3 = graph.add_node(weight.clone());
        let node4 = graph.add_node(weight.clone());
        let node5 = graph.add_node(weight.clone());

        let distance = node_distance(&node1.index(), &node3.index());
        assert_eq!(distance, 1);
        let distance_2 = node_distance(&node1.index(), &node2.index());
        assert_eq!(distance_2, 0);
        let distance_3 = node_distance(&node2.index(), &node5.index());
        assert_eq!(distance_3, 2);
        let distance_4 = node_distance(&node4.index(), &node5.index());
        assert_eq!(distance_4, 1);
    }
}
