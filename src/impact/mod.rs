use crate::graph::VariantGraph;

pub(crate) struct ImpactGraph {
    pub(crate) variant_graph: VariantGraph,
    pub(crate) exons: Vec<bio::io::gff::Record>,
}

#[derive(Debug, PartialEq, Eq, Ord, PartialOrd, Clone, Copy)]
pub(crate) enum Impact {
    None,
    Modifier,
    Low,
    Moderate,
    High,
}

#[cfg(test)]
mod impact_enum_tests {
    use super::Impact;

    #[test]
    fn test_impact_order() {
        assert!(Impact::High > Impact::Moderate);
        assert!(Impact::High > Impact::Low);
        assert!(Impact::High > Impact::Modifier);
        assert!(Impact::High > Impact::None);
        assert!(Impact::Moderate > Impact::Low);
        assert!(Impact::Moderate > Impact::Modifier);
        assert!(Impact::Moderate > Impact::None);
        assert!(Impact::Low > Impact::Modifier);
        assert!(Impact::Low > Impact::None);
        assert!(Impact::Modifier > Impact::None);
    }
}
