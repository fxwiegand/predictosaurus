use varlociraptor::calling::variants::preprocessing::Observations;
use varlociraptor::utils::collect_variants::VariantInfo;

#[derive(Debug)]
pub(crate) enum NodeType {
    Var,
    Ref,
}

#[derive(Debug)]
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
