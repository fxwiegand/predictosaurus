use crate::cli::Predictosaurus;
use crate::graph::{Edge, Node, NodeType};
use crate::utils::bcf::extract_event_names;
use anyhow::Result;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::bayes_factors::BayesFactor;
use clap::Parser;
use petgraph::Directed;
use petgraph::Graph;
use rust_htslib::bcf::{Read, Reader};
use std::collections::HashMap;
use varlociraptor::calling::variants::preprocessing::read_observations;
use varlociraptor::utils::collect_variants::collect_variants;
use itertools::Itertools;
use petgraph::dot::{Config, Dot};

mod cli;
mod graph;
mod utils;

fn main() -> Result<()> {
    let args = Predictosaurus::parse();

    let calls_file = args.calls;
    let observations_file = args.observations;

    let mut calls_reader = Reader::from_path(&calls_file)?;
    let mut observations_reader = Reader::from_path(observations_file)?;

    let _event_names = extract_event_names(&calls_file);

    let mut supporting_reads = HashMap::new();

    let mut variant_graph = Graph::<Node, Edge, Directed>::new();

    // Idea: Iterate in batches of records that have a near position. This means whenever the next variant is more than for example 1000bp away, we can stop the batch - write out the existing graph - and start a new one.
    for (calls_record, observations_record) in
        calls_reader.records().zip(observations_reader.records())
    {
        let mut calls_record = calls_record?;
        let mut observations_record = observations_record?;

        let variants = collect_variants(&mut calls_record, false, None)?;
        let observations = read_observations(&mut observations_record)?;

        let var_node = Node::from_observations(variants.first().unwrap(), &observations, NodeType::Var);
        let var_node_index = variant_graph.add_node(var_node);

        let ref_node = Node::from_observations(variants.first().unwrap(), &observations, NodeType::Ref);
        let ref_node_index = variant_graph.add_node(ref_node);

        for observation in observations.pileup.read_observations() {
            let evidence = BayesFactor::new(observation.prob_alt, observation.prob_ref)
                .evidence_kass_raftery();
            match evidence {
                KassRaftery::Strong | KassRaftery::VeryStrong => {
                    // Read supports variant
                    let entry = supporting_reads.entry(observation.fragment_id).or_insert(Vec::new());
                    entry.push(var_node_index);
                }
                KassRaftery::None | KassRaftery::Barely | KassRaftery::Positive => {
                    // Read supports reference
                    let entry = supporting_reads.entry(observation.fragment_id).or_insert(Vec::new());
                    entry.push(ref_node_index);
                }
            }
        }
    }

    for (_, nodes) in supporting_reads {
        for node_tuple in nodes.into_iter().combinations(2) {
            let edge = variant_graph.find_edge(node_tuple[0], node_tuple[1]);
            if let Some(edge) = edge {
                let edge = variant_graph.edge_weight_mut(edge).unwrap();
                edge.supporting_reads += 1;
            } else {
                let edge = Edge {
                    supporting_reads: 1,
                };
                variant_graph.add_edge(node_tuple[0], node_tuple[1], edge);
            }
        }
    }

    println!("{:?}", Dot::with_config(&variant_graph, &[Config::EdgeNoLabel]));

    Ok(())
}
