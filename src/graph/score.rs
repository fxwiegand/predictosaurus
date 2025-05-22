use crate::graph::node::{Node, NodeType};
use crate::graph::shift_phase;
use crate::graph::transcript::Transcript;
use crate::translation::amino_acids::AminoAcid;
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
    /// Number of coding SNVs on this haplotype
    pub num_snvs: usize,
    /// Total fractional CDS length affected by frameshifts (0.0–1.0)
    pub fs_fraction: f64,
    /// Stop-gained penalty
    pub stop_fraction: Option<f64>,
}

impl EffectScore {
    /// Create a new empty score
    pub fn new() -> Self {
        EffectScore {
            num_snvs: 0,
            fs_fraction: 0.0,
            stop_fraction: None,
        }
    }

    pub(crate) fn from_haplotype(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
        haplotype: Vec<Node>,
    ) -> Result<Self> {
        let mut num_snvs = 0;
        let mut stop_penalty = None;
        let mut phase = 0;
        let mut cds = None;
        let mut frameshift_positions = Vec::new();
        for node in haplotype.iter().filter(|n| n.node_type.is_variant()) {
            if node.is_snv() {
                num_snvs += 1;
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
                    .variant_amino_acids(
                        phase,
                        reference.get(&transcript.target).unwrap(),
                        transcript.strand,
                    )
                    .unwrap()
                    .contains(&AminoAcid::Methionine)
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
            num_snvs,
            fs_fraction,
            stop_fraction: stop_penalty,
        })
    }

    /// Compute the raw combined score using tuning constants
    pub fn raw(&self) -> f64 {
        FS_WEIGHT * self.fs_fraction
            + SNV_WEIGHT * (self.num_snvs as f64)
            + STOP_WEIGHT * self.stop_fraction.unwrap_or(0.0)
    }

    /// Compute the normalized score in (0,1): raw/(1+raw)
    pub fn normalized(&self) -> f64 {
        let r = self.raw();
        r / (1.0 + r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::EventProbs;
    use crate::Cds;

    #[test]
    fn test_raw_score() {
        let score = EffectScore {
            num_snvs: 3,
            fs_fraction: 0.4,
            stop_fraction: Some(0.2),
        };
        let expected = FS_WEIGHT * 0.4 + SNV_WEIGHT * 3.0 + STOP_WEIGHT * 0.2;
        assert!((score.raw() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_normalized_score() {
        let score = EffectScore {
            num_snvs: 2,
            fs_fraction: 0.1,
            stop_fraction: Some(0.3),
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

        assert_eq!(result.num_snvs, 2);
        assert!((result.fs_fraction - 0.75).abs() < 1e-6);
        assert!(result.stop_fraction.is_none());
    }
}
