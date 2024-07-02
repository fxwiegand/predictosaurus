use crate::graph::VariantGraph;

pub(crate) struct ImpactGraph {
    pub(crate) variant_graph: VariantGraph,
    pub(crate) exons: Vec<bio::io::gff::Record>,
}
