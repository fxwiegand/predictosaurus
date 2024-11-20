use crate::graph::node::{Node, NodeType};
use crate::graph::paths::Weight;
use crate::graph::{Edge, VariantGraph};
use anyhow::Result;
use duckdb::{params, Connection};
use petgraph::matrix_graph::NodeIndex;
use petgraph::Graph;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::str::FromStr;

pub(crate) fn write_graphs(
    graphs: HashMap<String, VariantGraph>,
    output_path: &Path,
) -> Result<()> {
    let path = output_path.join("graphs.duckdb");
    let db = Connection::open(&path)?;
    db.execute("CREATE TABLE graphs (target STRING PRIMARY KEY, start_position INTEGER, end_position INTEGER)", [])?;
    db.execute(
        "CREATE TABLE nodes (target STRING, node_index INTEGER PRIMARY KEY, node_type STRING, vaf STRING, probs STRING, pos INTEGER, index INTEGER)",
        [],
    )?;
    db.execute(
        "CREATE TABLE edges (target STRING, from_node INTEGER, to_node INTEGER, supporting_reads STRING)",
        [],
    )?;
    for (target, graph) in graphs {
        db.execute(
            "INSERT INTO graphs VALUES (?, ?, ?)",
            [
                target.to_string(),
                graph.start.to_string(),
                graph.end.to_string(),
            ],
        )?;
        for node_index in graph.graph.node_indices() {
            let node = graph.graph.node_weight(node_index).unwrap();
            db.execute(
                "INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?, ?)",
                [
                    &target.to_string(),
                    &node_index.index().to_string(),
                    &node.node_type.to_string(),
                    &serde_json::to_string(&node.vaf)?,
                    &serde_json::to_string(&node.probs)?,
                    &node.pos.to_string(),
                    &node.index.to_string(),
                ],
            )?;
        }
        for edge in graph.graph.edge_indices() {
            let (from_node, to_node) = graph.graph.edge_endpoints(edge).unwrap();
            let edge = graph.graph.edge_weight(edge).unwrap();
            db.execute(
                "INSERT INTO edges VALUES (?, ?, ?, ?)",
                [
                    &target.to_string(),
                    &from_node.index().to_string(),
                    &to_node.index().to_string(),
                    &serde_json::to_string(&edge.supporting_reads)?,
                ],
            )?;
        }
    }
    db.close().unwrap();
    Ok(())
}

pub(crate) fn feature_graph(
    path: PathBuf,
    target: String,
    start: u64,
    end: u64,
) -> Result<VariantGraph> {
    let db = Connection::open(path)?;
    let mut graph = Graph::<Node, Edge, petgraph::Directed>::new();
    let mut stmt = db.prepare(
    "SELECT node_index, node_type, vaf, probs, pos, index FROM nodes WHERE target = ? AND pos >= ? AND pos <= ?"
    )?;
    let nodes: Vec<(usize, String, String, String, i64, u32)> = stmt
        .query_map(
            params![target.to_string(), start as i64, end as i64],
            |row| {
                Ok((
                    row.get(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get(4)?,
                    row.get(5)?,
                ))
            },
        )?
        .map(Result::unwrap)
        .collect();
    for _ in 0..nodes.len() {
        graph.add_node(Node::new(NodeType::Ref("".to_string()), 0)); // Add placeholder nodes
    }
    for (node_index, node_type, vaf, probs, pos, index) in nodes {
        let node = Node {
            node_type: NodeType::from_str(&node_type).unwrap(),
            vaf: serde_json::from_str(&vaf).unwrap(),
            probs: serde_json::from_str(&probs).unwrap(),
            pos,
            index,
        };
        graph[NodeIndex::new(node_index)] = node;
    }
    let mut stmt = db
        .prepare("SELECT from_node, to_node, supporting_reads FROM edges WHERE target = ?")
        .unwrap();
    let edges: Vec<(usize, usize, String)> = stmt
        .query_map(params![target.to_string()], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?))
        })?
        .map(Result::unwrap)
        .collect();
    for (from_node, to_node, supporting_reads) in edges {
        if graph.node_indices().nth(from_node).is_none()
            || graph.node_indices().nth(to_node).is_none()
        {
            continue;
        }
        let edge_weight: HashMap<String, u32> = serde_json::from_str(&supporting_reads).unwrap();
        graph.add_edge(
            graph.node_indices().nth(from_node).unwrap(),
            graph.node_indices().nth(to_node).unwrap(),
            Edge {
                supporting_reads: edge_weight,
            },
        );
    }
    Ok(VariantGraph {
        graph,
        start: start as i64,
        end: end as i64,
        target: target.clone(),
    })
}

pub(crate) fn create_paths(output_path: &Path) -> Result<()> {
    let path = output_path.join("paths.duckdb");
    let db = Connection::open(&path)?;
    // TODO: Rethink the schema, we need a unique identifier for each CDS that we then use to iter batchwise via GROUP BY to generate the plots with the show command
    db.execute("CREATE TABLE path_nodes (path_index INTEGER, target String, feature STRING, node_index INTEGER, vaf FLOAT, impact STRING, reason STRING, consequence STRING, sample STRING)", [])?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn write_paths(
    outputh_path: &Path,
    paths: Vec<Vec<Weight>>,
    target: String,
    feature: String,
) -> Result<()> {
    let path = outputh_path.join("paths.duckdb");
    let db = Connection::open(&path)?;
    for (index, path) in paths.iter().enumerate() {
        for weight in path {
            db.execute(
                "INSERT INTO path_nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    index.to_string(),
                    target.to_string(),
                    feature.to_string(),
                    weight.index.to_string(),
                    weight.vaf.to_string(),
                    weight.impact.to_string(),
                    weight.reason.to_string(),
                    weight.consequence.to_string(),
                    weight.sample.to_string(),
                ],
            )?;
        }
    }
    db.close().unwrap();
    Ok(())
}

mod tests {
    use super::*;
    use crate::graph::node::{Node, NodeType};
    use crate::graph::Edge;
    use crate::impact::Impact;
    use petgraph::{Directed, Graph};

    pub(crate) fn setup_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node1 = graph.add_node(Node::new(NodeType::Var("A".to_string()), 1));
        let node2 = graph.add_node(Node::new(NodeType::Ref("".to_string()), 2));
        let node3 = graph.add_node(Node::new(NodeType::Var("T".to_string()), 3));
        let node4 = graph.add_node(Node::new(NodeType::Var("".to_string()), 4));
        let node5 = graph.add_node(Node::new(NodeType::Var("A".to_string()), 8));
        let node6 = graph.add_node(Node::new(NodeType::Var("TT".to_string()), 9));
        let edge1 = graph.add_edge(
            node1,
            node2,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let edge2 = graph.add_edge(
            node2,
            node3,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let edge3 = graph.add_edge(
            node3,
            node4,
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

    #[test]
    fn test_write_graphs() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path();
        write_graphs(graphs, output_path).unwrap();
        let db = Connection::open(output_path.join("graphs.duckdb")).unwrap();
        let mut stmt = db.prepare("SELECT target FROM graphs").unwrap();
        let targets: Vec<String> = stmt
            .query_map([], |row| row.get(0))
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(targets, vec!["graph1"]);
        let mut stmt = db
            .prepare("SELECT node_index, node_type, vaf, probs, pos, index FROM nodes")
            .unwrap();
        let nodes: Vec<(i64, String, String, String, i64, i64)> = stmt
            .query_map([], |row| {
                Ok((
                    row.get(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get(4)?,
                    row.get(5)?,
                ))
            })
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(nodes.len(), 6);
        let mut stmt = db
            .prepare("SELECT from_node, to_node, supporting_reads FROM edges")
            .unwrap();
        let edges: Vec<(i64, i64, String)> = stmt
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(edges.len(), 3);
        db.close().unwrap();
    }

    #[test]
    fn test_feature_graphs() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path();
        write_graphs(graphs, output_path).unwrap();
        let graph = feature_graph(
            output_path.join("graphs.duckdb"),
            "graph1".to_string(),
            1,
            3,
        )
        .unwrap();
        assert_eq!(graph.graph.node_count(), 3);
        assert_eq!(graph.graph.edge_count(), 2);
        assert_eq!(graph.start, 1);
        assert_eq!(graph.end, 3);
    }

    #[test]
    fn test_create_paths() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path();
        assert!(create_paths(output_path).is_ok());
        assert!(Connection::open(output_path.join("graphs.duckdb")).is_ok());
    }

    #[test]
    fn test_write_paths_creates_database_and_inserts_paths() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path();
        create_paths(output_path).unwrap();
        let paths = vec![
            vec![Weight {
                index: 1,
                vaf: 0.5,
                impact: Impact::High,
                reason: "Ile -> Met".to_string(),
                consequence: "loss".to_string(),
                sample: "sample1".to_string(),
            }],
            vec![Weight {
                index: 2,
                vaf: 0.3,
                impact: Impact::Low,
                reason: "Met -> Lys".to_string(),
                consequence: "gain".to_string(),
                sample: "sample2".to_string(),
            }],
        ];
        write_paths(
            output_path,
            paths,
            "1".to_string(),
            "some feature".to_string(),
        )
        .unwrap();
        let db = Connection::open(output_path.join("paths.duckdb")).unwrap();
        let mut stmt = db.prepare("SELECT path_index, target, feature, node_index, vaf, impact, reason, consequence, sample FROM path_nodes").unwrap();
        let path_nodes: Vec<(
            i64,
            String,
            String,
            i64,
            f64,
            String,
            String,
            String,
            String,
        )> = stmt
            .query_map([], |row| {
                Ok((
                    row.get(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get(4)?,
                    row.get(5)?,
                    row.get(6)?,
                    row.get(7)?,
                    row.get(8)?,
                ))
            })
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(path_nodes.len(), 2);
        db.close().unwrap();
    }
}
