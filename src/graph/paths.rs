use petgraph::graph::NodeIndex;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct HaplotypePath(pub(crate) Vec<NodeIndex>);

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq, Hash)]
pub(crate) struct Cds {
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) phase: u8,
}

impl Cds {
    pub(crate) fn new(start: u64, end: u64, phase: u8) -> Cds {
        Cds { start, end, phase }
    }

    /// Returns true if the CDS contains the variant
    pub(crate) fn contains_variant(
        &self,
        target: &str,
        variants: &HashMap<String, BTreeSet<i64>>,
    ) -> bool {
        if let Some(variant_positions) = variants.get(target) {
            variant_positions
                .range(self.start as i64..=self.end as i64)
                .next()
                .is_some()
        } else {
            false
        }
    }

    /// Length of CDS segment
    pub fn length(&self) -> usize {
        (self.end - self.start + 1) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    

    #[test]
    fn contains_variant_returns_true_when_variant_in_range() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([15]));
        assert!(cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_false_when_variant_out_of_range() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([25]));
        assert!(!cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_false_when_no_variants() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let variants: HashMap<String, BTreeSet<i64>> = HashMap::new();
        assert!(!cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_true_when_multiple_variants_in_range() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([12, 18]));
        assert!(cds.contains_variant(target, &variants));
    }

    #[test]
    fn contains_variant_returns_true_when_variants_on_boundary() {
        let cds = Cds::new(10, 20, 0);
        let target = "chr1";
        let mut variants = HashMap::new();
        variants.insert(target.to_string(), BTreeSet::from([10, 20]));
        assert!(cds.contains_variant(target, &variants));
    }
}
