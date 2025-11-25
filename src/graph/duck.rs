use crate::graph::node::{Node, NodeType};
use crate::graph::paths::{Cds, Weight};
use crate::graph::score::EffectScore;
use crate::graph::score::HaplotypeLikelihoods;
use crate::graph::transcript::Transcript;
use crate::graph::{Edge, VariantGraph};
use crate::impact::Impact;
use anyhow::Result;
use duckdb::{params, Connection};
use petgraph::matrix_graph::NodeIndex;
use petgraph::Graph;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::str::FromStr;

pub(crate) fn write_graphs(graphs: HashMap<String, VariantGraph>, path: &Path) -> Result<()> {
    let mut db = Connection::open(path)?;
    let transaction = db.transaction()?;
    transaction.execute("CREATE TABLE graphs (target STRING PRIMARY KEY, start_position INTEGER, end_position INTEGER)", [])?;
    transaction.execute(
        "CREATE TABLE nodes (target STRING, node_index INTEGER, node_type STRING, reference_allele STRING, alternative_allele STRING, vaf STRING, probs STRING, pos INTEGER, index INTEGER)",
        [],
    )?;
    transaction.execute(
        "CREATE TABLE edges (target STRING, from_node INTEGER, to_node INTEGER, supporting_reads STRING)",
        [],
    )?;

    let mut insert_graph = transaction.prepare("INSERT INTO graphs VALUES (?, ?, ?)")?;
    let mut insert_node =
        transaction.prepare("INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)")?;
    let mut insert_edge = transaction.prepare("INSERT INTO edges VALUES (?, ?, ?, ?)")?;

    for (target, graph) in graphs {
        insert_graph.execute(params![target, graph.start, graph.end])?;
        for node_index in graph.graph.node_indices() {
            let node = graph.graph.node_weight(node_index).unwrap();
            insert_node.execute(params![
                target,
                node_index.index() as i64,
                node.node_type.to_string(),
                node.reference_allele,
                node.alternative_allele,
                json5::to_string(&node.vaf)?,
                json5::to_string(&node.probs)?,
                node.pos,
                node.index as i64,
            ])?;
        }
        for edge in graph.graph.edge_indices() {
            let (from_node, to_node) = graph.graph.edge_endpoints(edge).unwrap();
            let edge = graph.graph.edge_weight(edge).unwrap();
            insert_edge.execute(params![
                target,
                from_node.index() as i64,
                to_node.index() as i64,
                serde_json::to_string(&edge.supporting_reads)?,
            ])?;
        }
    }
    transaction.execute(
        "CREATE INDEX idx_nodes_target_pos ON nodes (target, pos)",
        [],
    )?;
    transaction.execute("CREATE INDEX idx_edges_target ON edges (target)", [])?;
    transaction.execute("CREATE INDEX idx_graphs_target ON graphs (target)", [])?;
    transaction.commit()?;
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
    "SELECT node_index, node_type, reference_allele, alternative_allele, vaf, probs, pos, index FROM nodes WHERE target = ? AND pos >= ? AND pos <= ?"
    )?;
    let nodes: Vec<(usize, String, String, String, String, String, i64, u32)> = stmt
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
                    row.get(6)?,
                    row.get(7)?,
                ))
            },
        )?
        .map(Result::unwrap)
        .collect();
    let max_index = *nodes
        .iter()
        .map(|(index, _, _, _, _, _, _, _)| index)
        .max()
        .unwrap_or(&0);
    for _ in 0..=max_index {
        graph.add_node(Node::new(
            NodeType::Reference,
            -1,
            "".to_string(),
            "".to_string(),
        )); // Add placeholder nodes
    }
    for (node_index, node_type, reference_allele, alternative_allele, vaf, probs, pos, index) in
        nodes
    {
        let node = Node {
            node_type: NodeType::from_str(&node_type)?,
            reference_allele,
            alternative_allele,
            vaf: json5::from_str(&vaf)?,
            probs: json5::from_str(&probs)?,
            pos,
            index,
        };
        graph[NodeIndex::new(node_index)] = node;
    }
    let mut stmt =
        db.prepare("SELECT from_node, to_node, supporting_reads FROM edges WHERE target = ?")?;
    let edges: Vec<(usize, usize, String)> = stmt
        .query_map(params![target.to_string()], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?))
        })?
        .collect::<Result<Vec<_>, _>>()?;
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
    let temp_graph = graph.clone();
    graph.retain_nodes(|_, node| temp_graph.node_weight(node).unwrap().pos != -1);
    Ok(VariantGraph {
        graph,
        start: start as i64,
        end: end as i64,
        target: target.clone(),
    })
}

/// Retrieves all variants from a given serialised graph file. Returns a HashMap of targets and all positions of variants in the graph.
pub(crate) fn variants_on_graph(path: &PathBuf) -> Result<HashMap<String, BTreeSet<i64>>> {
    let db = Connection::open(path)?;
    let mut stmt = db.prepare("SELECT target, pos FROM nodes")?;
    let variants: Vec<(String, i64)> = stmt
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .map(Result::unwrap)
        .collect();
    let mut variant_map = HashMap::new();
    for (target, pos) in variants {
        let target_variants = variant_map.entry(target).or_insert(BTreeSet::new());
        target_variants.insert(pos);
    }
    Ok(variant_map)
}

pub(crate) fn create_scores(output_path: &Path) -> Result<()> {
    let db = Connection::open(output_path)?;
    db.execute(
        "CREATE TABLE scores (transcript String, score FLOAT, likelihoods String)",
        [],
    )?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn write_scores(
    path: &Path,
    scores: Vec<(EffectScore, HaplotypeLikelihoods)>,
    transcript: Transcript,
) -> Result<()> {
    let mut db = Connection::open(path)?;
    let transaction = db.transaction()?;
    let mut stmt = transaction.prepare("INSERT INTO scores VALUES (?, ?, ?)")?;
    let transcript_name = transcript.name();
    for (score, likelihoods) in scores {
        stmt.execute(params![
            transcript_name.as_str(),
            score.score(),
            json5::to_string(&likelihoods)?
        ])?;
    }
    transaction.commit()?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn read_scores(
    path: &Path,
) -> Result<HashMap<String, Vec<(f64, HaplotypeLikelihoods)>>> {
    let db = Connection::open(path)?;
    let mut scores = HashMap::new();
    let mut stmt = db.prepare("SELECT transcript, score, likelihoods FROM scores")?;
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let transcript = row.get(0)?;
        let score = row.get(1)?;
        let likelihoods: String = row.get(2)?;
        scores
            .entry(transcript)
            .or_insert(Vec::new())
            .push((score, json5::from_str(&likelihoods)?));
    }
    db.close().unwrap();
    Ok(scores)
}

pub(crate) fn create_paths(output_path: &Path) -> Result<()> {
    let db = Connection::open(output_path)?;
    db.execute("CREATE TABLE path_nodes (path_index INTEGER, target String, feature STRING, node_index INTEGER, vaf FLOAT, impact STRING, reason STRING, consequence STRING, sample STRING)", [])?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn write_paths(
    path: &Path,
    paths: Vec<Vec<Weight>>,
    transcript: Transcript,
) -> Result<()> {
    let db = Connection::open(path)?;
    for (index, path) in paths.iter().enumerate() {
        for weight in path {
            db.execute(
                "INSERT INTO path_nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    index.to_string(),
                    transcript.target.to_string(),
                    transcript.feature.to_string(),
                    weight.index.to_string(),
                    weight.vaf.to_string(),
                    weight.impact.to_raw_string().parse()?,
                    weight.reason.as_deref().unwrap_or("None").to_string(),
                    weight.consequence.as_deref().unwrap_or("None").to_string(),
                    weight.sample.to_string(),
                ],
            )?;
        }
    }
    db.close().unwrap();
    Ok(())
}

/// Reads paths and their metadata from the database file and returns a HashMap
/// where each feature is mapped to a HashMap of path indices and their corresponding
/// vectors of Weight values.
///
/// # Arguments
/// * `path` - The path to the DuckDB file to read from.
///
/// # Returns
/// A Result containing a HashMap with the CDS ID as the key, and a HashMap of path indices
/// mapped to vectors of Weight structs for each feature.
pub(crate) fn read_paths(path: &Path) -> Result<HashMap<String, Vec<Weight>>> {
    let db = Connection::open(path)?;
    let mut stmt = db.prepare(
        "SELECT path_index, feature, node_index, vaf, impact, reason, consequence, sample
         FROM path_nodes",
    )?;

    let mut paths: HashMap<String, Vec<Weight>> = HashMap::new();
    let rows = stmt.query_map([], |row| {
        let transcript_id = row.get::<_, String>(1)?;
        let weight = Weight {
            index: row.get(2)?,
            path: Some(row.get(0)?),
            vaf: row.get(3)?,
            impact: Impact::from_str(&row.get::<_, String>(4)?).unwrap(),
            reason: row.get(5)?,
            consequence: row.get(6)?,
            sample: row.get(7)?,
        };
        Ok((transcript_id, weight))
    })?;

    for row in rows {
        let (cds, weight) = row?;
        let feature_entry = paths.entry(cds).or_default();
        feature_entry.push(weight);
    }

    Ok(paths)
}

mod tests {
    use super::*;
    use crate::graph::node::{Node, NodeType};
    use crate::graph::Edge;
    use crate::impact::Impact;
    use crate::translation::distance::DistanceMetric;
    use bio::bio_types::strand::Strand;
    use itertools::Itertools;
    use petgraph::{Directed, Graph};

    pub(crate) fn setup_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node1 = graph.add_node(Node::new(
            NodeType::Variant,
            1,
            "T".to_string(),
            "A".to_string(),
        ));
        let node2 = graph.add_node(Node::new(
            NodeType::Reference,
            2,
            "".to_string(),
            "".to_string(),
        ));
        let node3 = graph.add_node(Node::new(
            NodeType::Variant,
            3,
            "A".to_string(),
            "T".to_string(),
        ));
        let node4 = graph.add_node(Node::new(
            NodeType::Variant,
            4,
            "T".to_string(),
            "".to_string(),
        ));
        let _node5 = graph.add_node(Node::new(
            NodeType::Variant,
            8,
            "C".to_string(),
            "A".to_string(),
        ));
        let node6 = graph.add_node(Node::new(
            NodeType::Variant,
            9,
            "A".to_string(),
            "TT".to_string(),
        ));
        let _edge1 = graph.add_edge(
            node1,
            node2,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let _edge2 = graph.add_edge(
            node2,
            node3,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let _edge3 = graph.add_edge(
            node3,
            node4,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let _edge4 = graph.add_edge(
            node4,
            node6,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        VariantGraph {
            graph,
            start: 0,
            end: 10,
            target: "test".to_string(),
        }
    }

    #[test]
    fn test_write_graphs() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let db = Connection::open(output_path.as_path()).unwrap();
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
        assert_eq!(edges.len(), 4);
        db.close().unwrap();
    }

    #[test]
    fn test_write_graph_with_multiple_targets() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        graphs.insert("graph2".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let db = Connection::open(output_path.as_path()).unwrap();
        let mut stmt = db.prepare("SELECT target FROM graphs").unwrap();
        let targets: Vec<String> = stmt
            .query_map([], |row| row.get(0))
            .unwrap()
            .map(Result::unwrap)
            .sorted_by(|a: &String, b| a.cmp(b))
            .collect();
        assert_eq!(targets, vec!["graph1", "graph2"]);
    }

    #[test]
    fn test_feature_graphs() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let graph = feature_graph(output_path, "graph1".to_string(), 1, 3).unwrap();
        assert_eq!(graph.graph.node_count(), 3);
        assert_eq!(graph.graph.edge_count(), 2);
        assert_eq!(graph.start, 1);
        assert_eq!(graph.end, 3);
    }

    #[test]
    fn test_feature_graphs_2() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let graph = feature_graph(output_path, "graph1".to_string(), 3, 10).unwrap();
        assert_eq!(graph.graph.node_count(), 4);
        assert_eq!(graph.graph.edge_count(), 2);
        dbg!(&graph);
        assert_eq!(graph.start, 3);
        assert_eq!(graph.end, 10);
    }

    #[test]
    fn test_create_paths() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("paths.duckdb");
        assert!(create_paths(output_path.as_path()).is_ok());
        assert!(Connection::open(output_path.as_path()).is_ok());
    }

    #[test]
    fn test_write_paths_creates_database_and_inserts_paths() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("paths.duckdb");
        create_paths(output_path.as_path()).unwrap();
        let paths = vec![
            vec![Weight {
                index: 1,
                path: None,
                vaf: 0.5,
                impact: Impact::High,
                reason: Some("Ile -> Met".to_string()),
                consequence: Some("loss".to_string()),
                sample: "sample1".to_string(),
            }],
            vec![Weight {
                index: 2,
                path: None,
                vaf: 0.3,
                impact: Impact::Low,
                reason: Some("Met -> Lys".to_string()),
                consequence: Some("gain".to_string()),
                sample: "sample2".to_string(),
            }],
        ];
        let transcript = Transcript::new(
            "1".to_string(),
            "some feature".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 100, 0)],
        );
        write_paths(output_path.as_path(), paths, transcript).unwrap();
        let db = Connection::open(output_path.as_path()).unwrap();
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

    #[test]
    fn read_paths_returns_correct_structure() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("paths.duckdb");
        create_paths(output_path.as_path()).unwrap();
        let paths = vec![
            vec![Weight {
                index: 1,
                path: None,
                vaf: 0.5,
                impact: Impact::High,
                reason: Some("Ile -> Met".to_string()),
                consequence: Some("loss".to_string()),
                sample: "sample1".to_string(),
            }],
            vec![Weight {
                index: 2,
                path: None,
                vaf: 0.3,
                impact: Impact::Low,
                reason: Some("Met -> Lys".to_string()),
                consequence: Some("gain".to_string()),
                sample: "sample2".to_string(),
            }],
        ];
        let transcript = Transcript::new(
            "some feature".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 100, 0)],
        );
        write_paths(output_path.as_path(), paths, transcript).unwrap();
        let result = read_paths(output_path.as_path()).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("some feature"));
        let feature_paths = result.get("some feature").unwrap();
        assert_eq!(feature_paths.len(), 2);
    }

    #[test]
    fn test_variants_on_graph() {
        let mut graphs = HashMap::new();
        graphs.insert("chr1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let result = variants_on_graph(&output_path).unwrap();
        let chr1_variants = result.get("chr1").unwrap();
        assert_eq!(chr1_variants, &vec![1, 2, 3, 4, 8, 9].into_iter().collect());
    }

    #[test]
    fn test_scores() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("scores.duckdb");
        create_scores(output_path.as_path()).unwrap();
        assert!(output_path.exists());
        let transcript = Transcript::new(
            "some feature".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 100, 0)],
        );
        let effect_score = EffectScore {
            snvs: Vec::new(),
            fs_fraction: 0.5,
            stop_fraction: None,
            distance_metric: DistanceMetric::default(),
        };
        let likelihoods = HashMap::from([("A".to_string(), 0.1), ("C".to_string(), 0.2)]);
        let scores = vec![(effect_score, likelihoods)];
        write_scores(output_path.as_path(), scores, transcript).unwrap();
        let scores = read_scores(output_path.as_path()).unwrap();
        assert_eq!(scores.len(), 1);
        assert!((scores.get("chr1:some feature").unwrap()[0].0 - 0.3333333).abs() < 1e-6);
        assert!(
            (scores.get("chr1:some feature").unwrap()[0]
                .1
                .get("A")
                .unwrap()
                - 0.1)
                .abs()
                < 1e-6
        );
        assert!(
            (scores.get("chr1:some feature").unwrap()[0]
                .1
                .get("C")
                .unwrap()
                - 0.2)
                .abs()
                < 1e-6
        );
    }
}
