use crate::graph::node::{Node, NodeType};
use crate::graph::shift_phase;
use crate::graph::transcript::Transcript;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::distance::DistanceMetric;
use anyhow::Result;
use bio::bio_types::strand::Strand;
use std::collections::HashMap;

/// Tuning constants for the score components
const SNV_WEIGHT: f64 = 0.1; // Weight per SNV
const FS_WEIGHT: f64 = 1.0; // Weight for frameshift fraction
const STOP_WEIGHT: f64 = 1.0; // Weight for stop-gained fraction

/// Breakdown of effects for one haplotype path
#[derive(Debug, Clone)]
pub struct EffectScore {
    /// SNVs on this haplotype
    pub snvs: Vec<AminoAcidChange>,
    /// Total fractional CDS length affected by frameshifts (0.0–1.0)
    pub fs_fraction: f64,
    /// Stop-gained penalty
    pub stop_fraction: Option<f64>,
    pub distance_metric: DistanceMetric,
}

impl EffectScore {
    /// Constructs a new `EffectScore` with no SNVs, zero frameshift and stop-gained fractions, and the default distance metric.
    ///
    /// # Examples
    ///
    /// ```
    /// let score = EffectScore::new();
    /// assert_eq!(score.snvs.len(), 0);
    /// assert_eq!(score.fs_fraction, 0.0);
    /// assert!(score.stop_fraction.is_none());
    /// ```
    pub fn new() -> Self {
        EffectScore {
            snvs: Vec::new(),
            fs_fraction: 0.0,
            stop_fraction: None,
            distance_metric: DistanceMetric::default(),
        }
    }

    /// Computes the cumulative effect score of a haplotype path on a transcript.
    ///
    /// Iterates over variant nodes in the haplotype to extract amino acid changes, frameshift events, and stop-gained variants. Calculates the fraction of the coding sequence affected by frameshifts and assigns a penalty for premature stop codons based on their position within the transcript. Returns an `EffectScore` summarizing these effects.
    ///
    /// # Parameters
    /// - `reference`: Reference genome sequences keyed by target name.
    /// - `transcript`: The transcript annotation to score against.
    /// - `haplotype`: Sequence of variant nodes representing the haplotype path.
    ///
    /// # Returns
    /// An `EffectScore` containing SNV amino acid changes, frameshift fraction, stop-gained penalty, and a default distance metric.
    ///
    /// # Errors
    /// Returns an error if amino acid extraction or transcript position calculations fail.
    ///
    /// # Examples
    ///
    /// ```
    /// let effect_score = EffectScore::from_haplotype(&reference, &transcript, haplotype)?;
    /// assert!(effect_score.fs_fraction >= 0.0 && effect_score.fs_fraction <= 1.0);
    /// ```
    pub(crate) fn from_haplotype(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
        haplotype: Vec<Node>,
    ) -> Result<Self> {
        let target_ref = reference.get(&transcript.target).unwrap();
        let mut snvs = Vec::new();
        let mut stop_penalty = None;
        let mut phase = 0;
        let mut cds = None;
        let mut frameshift_positions = Vec::new();
        for node in haplotype.iter().filter(|n| n.node_type.is_variant()) {
            if node.is_snv() {
                snvs.push(AminoAcidChange::from_node(
                    node,
                    phase,
                    target_ref,
                    transcript.strand,
                )?);
            } else if node.frameshift() != 0 {
                let pos = transcript.position_in_transcript(node.pos as usize)?;
                frameshift_positions.push((pos, node.frameshift()));
            }

            let node_cds = transcript
                .cds_for_position(node.pos)
                .expect("Found node located outside of CDS");
            if cds.as_ref() != Some(node_cds) {
                cds = Some(node_cds.to_owned());
                phase = node_cds.phase;
            }

            if stop_penalty.is_none()
                && node
                    .variant_amino_acids(phase, target_ref, transcript.strand)
                    .unwrap()
                    .contains(&AminoAcid::Stop)
            {
                // Calculate stop penalty based on position in Transcript. The earlier in the transcript, the higher the penalty
                let relative_position_fraction =
                    transcript.position_in_transcript(node.pos as usize)? as f64
                        / transcript.length() as f64;

                stop_penalty = Some(if transcript.strand == Strand::Forward {
                    1.0 - relative_position_fraction
                } else {
                    relative_position_fraction
                });
            }

            phase = shift_phase(phase, ((node.frameshift() + 3) % 3) as u8);
        }

        let (affected_length, remaining_fs, prev_pos) = frameshift_positions.iter().fold(
            (0i64, 0i64, 0usize), // (accumulated affected length, net frameshift, previous position)
            |(acc_len, net_fs, prev_pos), &(pos, fs)| {
                let affected_len = if net_fs % 3 != 0 {
                    // Calculate affected region length between prev_pos and current pos,
                    // respecting strand direction
                    let len = match transcript.strand {
                        Strand::Forward => pos as i64 - prev_pos as i64,
                        _ => prev_pos as i64 - pos as i64,
                    };
                    acc_len + len
                } else {
                    acc_len
                };
                let new_net_fs = net_fs + fs;
                (affected_len, new_net_fs, pos)
            },
        );

        // After fold, if net frameshift != 0, add the tail affected region till transcript end
        let fs_fraction = if remaining_fs % 3 != 0 {
            let tail_len = match transcript.strand {
                Strand::Forward => transcript.length() as i64 - prev_pos as i64,
                _ => prev_pos as i64 + 1,
            };
            ((affected_length + tail_len) as f64 / transcript.length() as f64)
        } else {
            affected_length as f64 / transcript.length() as f64
        };

        Ok(EffectScore {
            snvs,
            fs_fraction,
            stop_fraction: stop_penalty,
            distance_metric: DistanceMetric::default(),
        })
    }

    /// Calculates the total score for all single nucleotide variant (SNV) amino acid changes in the haplotype.
    ///
    /// Sums the distance metric values for each amino acid change recorded in the `snvs` field.
    ///
    /// # Returns
    /// The cumulative SNV score as a floating-point value.
    ///
    /// # Examples
    ///
    /// ```
    /// let score = effect_score.snv_score();
    /// assert!(score >= 0.0);
    /// ```
    fn snv_score(&self) -> f64 {
        self.snvs
            .iter()
            .map(|change| change.distance(&self.distance_metric))
            .sum()
    }

    /// Calculates the unnormalized combined effect score for a haplotype path.
    ///
    /// The score is computed as a weighted sum of the frameshift fraction, the total SNV (amino acid change) score, and the stop-gained penalty fraction, using predefined tuning constants.
    ///
    /// # Returns
    /// The raw (unnormalized) effect score as a floating-point value.
    pub fn raw(&self) -> f64 {
        FS_WEIGHT * self.fs_fraction
            + SNV_WEIGHT * self.snv_score()
            + STOP_WEIGHT * self.stop_fraction.unwrap_or(0.0)
    }

    /// Returns the normalized effect score as a value between 0 and 1.
    ///
    /// The normalized score is computed as `raw / (1 + raw)`, where `raw` is the unbounded effect score combining SNV, frameshift, and stop-gained penalties.
    ///
    /// # Examples
    ///
    /// ```
    /// let score = effect_score.normalized();
    /// assert!(score >= 0.0 && score < 1.0);
    /// ```
    pub fn normalized(&self) -> f64 {
        let r = self.raw();
        r / (1.0 + r)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AminoAcidChange {
    pub reference: Option<AminoAcid>,
    pub variants: Vec<AminoAcid>,
}

impl AminoAcidChange {
    /// Calculates the distance between the reference and variant amino acid using the provided metric.
    ///
    /// Returns the computed distance if there is exactly one variant and a reference amino acid; otherwise, returns a default penalty of 1.0.
    ///
    /// # Examples
    ///
    /// ```
    /// let change = AminoAcidChange {
    ///     reference: Some(AminoAcid::A),
    ///     variants: vec![AminoAcid::G],
    /// };
    /// let metric = DistanceMetric::default();
    /// let dist = change.distance(&metric);
    /// assert!(dist >= 0.0);
    /// ```
    pub fn distance(&self, metric: &DistanceMetric) -> f64 {
        match (&self.reference, self.variants.first()) {
            (Some(r), Some(v)) if self.variants.len() == 1 => metric.compute(r, v),
            _ => 1.0, // default penalty for complex changes
        }
    }

    /// Constructs an `AminoAcidChange` from a variant node, extracting the reference and variant amino acids based on the reading frame phase and strand orientation.
    ///
    /// Returns an error if amino acid extraction fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::{AminoAcidChange, Node, Strand};
    /// let node = /* obtain a Node representing a variant */;
    /// let reference = b"ATGGCC..."; // CDS reference sequence
    /// let phase = 0;
    /// let strand = Strand::Plus;
    /// let aa_change = AminoAcidChange::from_node(&node, phase, reference, strand).unwrap();
    /// ```
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::EventProbs;
    use crate::translation::amino_acids::AminoAcid;
    use crate::Cds;

    #[test]
    fn test_from_node() {
        let node = Node::new(NodeType::Var("A".to_string()), 2);
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
        assert!((change.distance(&metric) - 1.0).abs() < 1e-6);

        let change = AminoAcidChange {
            reference: Some(AminoAcid::Isoleucine),
            variants: vec![AminoAcid::AsparticAcid, AminoAcid::Threonine],
        };
        assert!((change.distance(&metric) - 1.0).abs() < 1e-6);

        let change = AminoAcidChange {
            reference: Some(AminoAcid::Isoleucine),
            variants: vec![],
        };
        assert!((change.distance(&metric) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_raw_score() {
        let aa_exchange = AminoAcidChange {
            reference: Some(AminoAcid::Isoleucine),
            variants: vec![AminoAcid::AsparticAcid],
        };
        let score = EffectScore {
            snvs: vec![aa_exchange], // 168/215 -> 0.78139534883
            fs_fraction: 0.4,
            stop_fraction: Some(0.2),
            distance_metric: DistanceMetric::Grantham,
        };
        println!("{}", score.snv_score());
        let expected = FS_WEIGHT * 0.4 + SNV_WEIGHT * 0.78139534883 + STOP_WEIGHT * 0.2;
        assert!((score.raw() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_normalized_score() {
        let score = EffectScore {
            snvs: vec![],
            fs_fraction: 0.1,
            stop_fraction: Some(0.3),
            distance_metric: DistanceMetric::Grantham,
        };
        let raw = score.raw();
        let expected = raw / (1.0 + raw);
        assert!((score.normalized() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_from_haplotype_with_frameshift() {
        let reference = HashMap::from([(
            "chr1".to_string(),
            b"ATGCGTACGTATGCGTACGTACGCGTACGTT".to_vec(),
        )]);

        let transcript = Transcript {
            feature: "test".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![
                Cds {
                    start: 1,
                    end: 10,
                    phase: 0,
                },
                Cds {
                    start: 21,
                    end: 30,
                    phase: 0,
                },
            ],
        };

        let haplotype = vec![
            Node::new(NodeType::Var("A".into()), 2),
            Node::new(NodeType::Var("TT".into()), 6),
            Node::new(NodeType::Var("C".into()), 22),
        ];

        let result = EffectScore::from_haplotype(&reference, &transcript, haplotype).unwrap();
        let aa_exchange = AminoAcidChange {
            reference: Some(AminoAcid::Methionine),
            variants: vec![AminoAcid::Isoleucine],
        };
        let aa_exchange_2 = AminoAcidChange {
            reference: Some(AminoAcid::Threonine),
            variants: vec![AminoAcid::Threonine],
        };

        assert_eq!(result.snvs, vec![aa_exchange, aa_exchange_2,]);
        assert!((result.fs_fraction - 0.75).abs() < 1e-6);
        assert!(result.stop_fraction.is_none());
    }
}
