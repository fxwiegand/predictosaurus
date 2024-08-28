use crate::graph::VariantGraph;
use duckdb::Connection;
use std::collections::HashMap;
use std::path::Path;

pub(crate) fn write_graphs(
    graphs: HashMap<String, VariantGraph>,
    output_path: &Path,
) -> anyhow::Result<()> {
    let path = output_path.join("graphs.duckdb");
    let db = Connection::open(&path)?;
    db.execute("CREATE TABLE graphs (target STRING PRIMARY KEY)", [])?;
    db.execute(
        "CREATE TABLE nodes (target STRING, node_index INTEGER PRIMARY KEY, node_type STRING, vaf STRING, probs STRING, pos INTEGER, index INTEGER)",
        [],
    )?;
    db.execute(
        "CREATE TABLE edges (target STRING, from_node INTEGER, to_node INTEGER, supporting_reads STRING)",
        [],
    )?;
    for (target, graph) in graphs {
        db.execute("INSERT INTO graphs VALUES (?)", [target.to_string()])?;
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

mod tests {
    use super::*;
    use crate::graph::Edge;
    use petgraph::graph::NodeIndex;

    #[test]
    fn test_write_graphs() {
        let mut graphs = HashMap::new();
        let mut g1 = crate::graph::tests::setup_variant_graph_with_nodes();
        g1.graph.add_edge(
            NodeIndex::new(0),
            NodeIndex::new(1),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graphs.insert("graph1".to_string(), g1);
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
        assert_eq!(edges.len(), 1);
        db.close().unwrap();
    }
}
