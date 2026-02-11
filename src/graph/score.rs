use crate::graph::node::{Node, NodeType};
use crate::graph::shift_phase;
use crate::graph::transcript::Transcript;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::amino_acids::Protein;
use crate::translation::distance::DistanceMetric;
use anyhow::Result;
use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation::*;
use bio::bio_types::strand::Strand;
use clap::ValueEnum;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct EffectScore {
    pub original_protein: Protein,
    pub altered_protein: Protein,
    pub distance_metric: DistanceMetric,
    pub realign: bool,
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
            let mut total = 0.0;
            for (i, (ref_aa, alt_aa)) in self
                .original_protein
                .amino_acids()
                .into_iter()
                .zip(self.altered_protein.amino_acids())
                .enumerate()
            {
                if alt_aa.is_stop() && !ref_aa.is_stop() {
                    total += (self.altered_protein.amino_acids().len() - i - 1) as f64;
                }
                total += self.distance_metric.compute(&ref_aa, &alt_aa);
            }
            total / self.altered_protein.amino_acids().len() as f64
        } else {
            // Realign using semiglobal
            let prot_x = self.original_protein.as_bytes();
            let prot_y = self.altered_protein.as_bytes();

            let score_fn = |a: u8, b: u8| {
                // Convert used distance matrix to a negative cost with integer scaling as expected by bio
                let d = self
                    .distance_metric
                    .compute(&AminoAcid::from(a), &AminoAcid::from(b));
                (-d * 10.0).round() as i32
            };

            let gap_open = -8;
            let gap_extend = -1;

            let mut aligner =
                Aligner::with_capacity(prot_x.len(), prot_y.len(), gap_open, gap_extend, &score_fn);

            let alignment = aligner.semiglobal(&prot_x, &prot_y);

            let ops = alignment.operations;

            let mut total = 0.0;

            let mut i = alignment.xstart;
            let mut j = alignment.ystart;

            for op in ops {
                match op {
                    Match => {
                        i += 1;
                        j += 1;
                    }
                    Subst => {
                        let aa_x_opt = prot_x.get(i).map(|&b| AminoAcid::from(b));
                        let aa_y_opt = prot_y.get(j).map(|&b| AminoAcid::from(b));

                        match (aa_x_opt, aa_y_opt) {
                            (Some(aa_x), Some(aa_y)) => {
                                if aa_y.is_stop() && !aa_x.is_stop() {
                                    let remaining = (prot_x.len() - i).max(prot_y.len() - j);
                                    total += remaining as f64;
                                    break;
                                }
                                total += self.distance_metric.compute(&aa_x, &aa_y);
                            }
                            _ => {
                                total += 1.0;
                            }
                        }

                        i += 1;
                        j += 1;
                    }
                    Del | Ins => {
                        total += 1.0;

                        if let Del = op {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                    _ => {}
                }
            }
            total / self.altered_protein.amino_acids().len() as f64
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, ValueEnum)]
pub enum HaplotypeMetric {
    Product,
    GeometricMean,
    Minimum,
}

pub type HaplotypeFrequency = HashMap<String, f32>;

impl HaplotypeMetric {
    pub fn calculate(&self, haplotype: &[Node]) -> HaplotypeFrequency {
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

    #[test]
    fn test_score_equal_proteins() {
        let p1 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let p2 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
        };
        assert!(score.score().abs() < 1e-6);
    }

    #[test]
    fn test_score_different_proteins() {
        let p1 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let p2 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Valine]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
        };
        // Distance between Leucine and Valine is 0.03 / 2 since len = 2
        assert!((score.score() - 0.015).abs() < 1e-6)
    }

    #[test]
    fn test_score_different_proteins_with_stop() {
        let p1 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        let p2 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Stop,
            AminoAcid::Valine,
        ]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
        };
        assert!((score.score() - (2.0 / 3.0)).abs() < 1e-6)
    }

    #[test]
    fn test_score_different_proteins_realign() {
        let p1 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        let p2 = Protein::new(vec![AminoAcid::Leucine, AminoAcid::Valine]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: true,
        };
        assert!((score.score() - 0.5).abs() < 1e-6)
    }

    #[test]
    fn test_score_different_proteins_realign_with_substitution() {
        let p1 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
            AminoAcid::Leucine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
        ]);
        let p2 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Valine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
        ]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: true,
        };
        assert!((score.score() - 0.2).abs() < 1e-6)
    }
}
