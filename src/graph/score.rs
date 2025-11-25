use crate::graph::node::{Node, NodeType};
use crate::graph::shift_phase;
use crate::graph::transcript::Transcript;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::amino_acids::Protein;
use crate::translation::distance::DistanceMetric;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use clap::ValueEnum;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct EffectScore {
    pub original_protein: Protein,
    pub altered_protein: Protein,
    pub distance_metric: DistanceMetric,
    realign: bool,
}

impl EffectScore {
    pub(crate) fn from_haplotype(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
        haplotype: &[Node],
        original_protein: Protein,
        distance_metric: DistanceMetric,
        realign: bool,
    ) -> Result<Self> {
        let altered_protein = Protein::from_haplotype(reference, transcript, haplotype)?;
        Ok(Self {
            original_protein,
            altered_protein,
            distance_metric,
            realign,
        })
    }

    pub fn score(&self) -> f64 {
        // Compare proteins:
        // If equal return 0
        // If no frameshift occured in the changed protein, simply compare amino acid by amino acid based on the distance metric and divide by the length of the protein
        // If a frameshift occurs, re-align the proteins and compare amino acid by amino acid based on the distance metric and divide by the length of the protein
        if self.original_protein == self.altered_protein {
            0.0
        } else if !self.realign {
            // We can assume both proteins have the same length
            self.original_protein
                .amino_acids()
                .into_iter()
                .zip(self.altered_protein.amino_acids())
                .map(|(a, b)| self.distance_metric.compute(&a, &b))
                .sum::<f64>()
                / self.original_protein.amino_acids().len() as f64
        } else {
            unimplemented!()
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AminoAcidChange {
    pub reference: Option<AminoAcid>,
    pub variants: Vec<AminoAcid>,
}

impl AminoAcidChange {
    pub fn distance(&self, metric: &DistanceMetric) -> f64 {
        match (&self.reference, self.variants.first(), self.variants.len()) {
            (Some(r), Some(v), 1) => metric.compute(r, v),
            (Some(_), Some(_), _) => 1.0, // Complex change adding more amino acids to the protein.
            _ => 0.0, // This will mostly happen when the variant is mapped to a region where the reference contains N. Therefore, we return 0.0 for now.
        }
    }

    pub fn from_node(
        node: &Node,
        phase: u8,
        reference: &[u8],
        strand: Strand,
    ) -> anyhow::Result<Self> {
        Ok(AminoAcidChange {
            reference: node.reference_amino_acid(phase, reference, strand)?,
            variants: node.variant_amino_acids(phase, reference, strand)?,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, ValueEnum)]
pub enum HaplotypeMetric {
    Product,
    GeometricMean,
    Minimum,
}

pub type HaplotypeLikelihoods = HashMap<String, f32>;

impl HaplotypeMetric {
    pub fn calculate(&self, haplotype: &[Node]) -> HaplotypeLikelihoods {
        let samples = haplotype
            .iter()
            .flat_map(|n| n.vaf.keys())
            .collect::<HashSet<_>>();
        let mut metrics = HashMap::new();
        for sample in samples {
            let vafs = haplotype
                .iter()
                .filter(|n| n.node_type.is_variant())
                .map(|n| *n.vaf.get(sample).unwrap_or(&0.0))
                .collect::<Vec<_>>();
            let result = match self {
                HaplotypeMetric::Product => vafs.iter().product(),
                HaplotypeMetric::GeometricMean => {
                    let product: f32 = vafs.iter().product();
                    product.powf(1.0 / vafs.len() as f32)
                }
                HaplotypeMetric::Minimum => vafs.iter().cloned().fold(1.0, f32::min),
            };
            metrics.insert(sample.to_string(), result);
        }
        metrics
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::EventProbs;
    use crate::translation::amino_acids::AminoAcid;
    use crate::Cds;

    #[test]
    fn test_from_node() {
        let node = Node::new(NodeType::Variant, 2, "G".to_string(), "A".to_string());
        let reference = b"ATGCGCGTA";
        let phase = 0;
        let strand = Strand::Forward;
        let change = AminoAcidChange::from_node(&node, phase, reference, strand).unwrap();
        assert_eq!(change.reference, Some(AminoAcid::Methionine));
        assert_eq!(change.variants, vec![AminoAcid::Isoleucine]);
    }

    #[test]
    fn test_amino_acid_change_distance() {
        use crate::translation::distance::DistanceMetric;

        let change = AminoAcidChange {
            reference: Some(AminoAcid::Isoleucine),
            variants: vec![AminoAcid::AsparticAcid],
        };
        let metric = DistanceMetric::Grantham;
        let expected = metric.compute(&AminoAcid::Isoleucine, &AminoAcid::AsparticAcid);
        assert!((change.distance(&metric) - expected).abs() < 1e-6);

        let change = AminoAcidChange {
            reference: None,
            variants: vec![AminoAcid::AsparticAcid],
        };
        assert!((change.distance(&metric) - 0.0).abs() < 1e-6);

        let change = AminoAcidChange {
            reference: Some(AminoAcid::Isoleucine),
            variants: vec![AminoAcid::AsparticAcid, AminoAcid::Threonine],
        };
        assert!((change.distance(&metric) - 1.0).abs() < 1e-6);

        let change = AminoAcidChange {
            reference: Some(AminoAcid::Isoleucine),
            variants: vec![],
        };
        assert!((change.distance(&metric) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_product_metric() {
        let nodes = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.5)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 0,
                index: 0,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.25)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 1,
            },
            Node {
                node_type: NodeType::Reference,
                reference_allele: "".to_string(),
                alternative_allele: "".to_string(),
                vaf: [("S1".to_string(), 0.75)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 2,
                index: 2,
            },
        ];
        let metric = HaplotypeMetric::Product;
        let result = metric.calculate(&nodes);
        assert_eq!(result.get("S1").unwrap(), &(0.5 * 0.25));
    }

    #[test]
    fn test_geometric_mean_metric() {
        let nodes = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.5)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 0,
                index: 0,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.25)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 1,
            },
        ];
        let metric = HaplotypeMetric::GeometricMean;
        let result = metric.calculate(&nodes);
        let expected = (0.5_f32 * 0.25_f32).powf(1.0 / 2.0);
        assert!((result.get("S1").unwrap() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_all_metrics_return_one_with_only_reference_nodes() {
        let nodes = vec![Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: [("S1".to_string(), 0.5)].into(),
            probs: EventProbs(HashMap::new()),
            pos: 0,
            index: 0,
        }];
        let metric = HaplotypeMetric::GeometricMean;
        let result = metric.calculate(&nodes);
        assert_eq!(result.get("S1").unwrap(), &1.0);
        let metric = HaplotypeMetric::Product;
        let result = metric.calculate(&nodes);
        assert_eq!(result.get("S1").unwrap(), &1.0);
        let metric = HaplotypeMetric::Minimum;
        let result = metric.calculate(&nodes);
        assert_eq!(result.get("S1").unwrap(), &1.0);
    }

    #[test]
    fn test_minimum_metric() {
        let nodes = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.5)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 0,
                index: 0,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.25)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 1,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.75)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 2,
                index: 2,
            },
        ];
        let metric = HaplotypeMetric::Minimum;
        let result = metric.calculate(&nodes);
        assert_eq!(*result.get("S1").unwrap(), 0.25);
    }
}
