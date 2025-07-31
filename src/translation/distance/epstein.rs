use crate::translation::amino_acids::AminoAcid;
use std::collections::HashMap;
use std::sync::LazyLock;

const EPSTEIN_ORDER: [AminoAcid; 20] = [
    AminoAcid::Phenylalanine,
    AminoAcid::Methionine,
    AminoAcid::Leucine,
    AminoAcid::Isoleucine,
    AminoAcid::Valine,
    AminoAcid::Proline,
    AminoAcid::Tyrosine,
    AminoAcid::Tryptophan,
    AminoAcid::Cysteine,
    AminoAcid::Alanine,
    AminoAcid::Glycine,
    AminoAcid::Serine,
    AminoAcid::Threonine,
    AminoAcid::Histidine,
    AminoAcid::GlutamicAcid,
    AminoAcid::Glutamine,
    AminoAcid::AsparticAcid,
    AminoAcid::Asparagine,
    AminoAcid::Lysine,
    AminoAcid::Arginine,
];

pub static EPSTEIN_INDEX: LazyLock<HashMap<AminoAcid, usize>> = LazyLock::new(|| {
    EPSTEIN_ORDER
        .iter()
        .enumerate()
        .map(|(i, aa)| (*aa, i))
        .collect()
});

// Table containing Epstein distance matrix
pub static EPSTEIN_MATRIX: LazyLock<[[f64; 20]; 20]> = LazyLock::new(|| {
    [
        [
            0.00, 0.05, 0.08, 0.08, 0.10, 0.10, 0.21, 0.25, 0.22, 0.43, 0.53, 0.81, 0.81, 0.80,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
        ],
        [
            0.10, 0.00, 0.03, 0.03, 0.10, 0.10, 0.25, 0.32, 0.21, 0.41, 0.42, 0.80, 0.80, 0.80,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
        ],
        [
            0.15, 0.05, 0.00, 0.00, 0.03, 0.03, 0.28, 0.36, 0.20, 0.43, 0.51, 0.80, 0.80, 0.81,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.01,
        ],
        [
            0.15, 0.05, 0.00, 0.00, 0.03, 0.03, 0.28, 0.36, 0.20, 0.43, 0.51, 0.80, 0.80, 0.81,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.01,
        ],
        [
            0.20, 0.10, 0.05, 0.05, 0.00, 0.00, 0.32, 0.40, 0.20, 0.40, 0.50, 0.80, 0.80, 0.81,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.02,
        ],
        [
            0.20, 0.10, 0.05, 0.05, 0.00, 0.00, 0.32, 0.40, 0.20, 0.40, 0.50, 0.80, 0.80, 0.81,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.02,
        ],
        [
            0.20, 0.22, 0.22, 0.22, 0.24, 0.24, 0.00, 0.10, 0.13, 0.27, 0.36, 0.62, 0.61, 0.60,
            0.80, 0.80, 0.81, 0.81, 0.80, 0.80,
        ],
        [
            0.21, 0.24, 0.25, 0.25, 0.27, 0.27, 0.05, 0.00, 0.18, 0.30, 0.39, 0.63, 0.63, 0.61,
            0.81, 0.81, 0.81, 0.81, 0.81, 0.80,
        ],
        [
            0.28, 0.22, 0.21, 0.21, 0.20, 0.20, 0.25, 0.35, 0.00, 0.25, 0.31, 0.60, 0.60, 0.62,
            0.81, 0.81, 0.80, 0.80, 0.81, 0.82,
        ],
        [
            0.50, 0.45, 0.43, 0.43, 0.41, 0.41, 0.40, 0.49, 0.22, 0.00, 0.10, 0.40, 0.41, 0.47,
            0.63, 0.63, 0.62, 0.62, 0.63, 0.67,
        ],
        [
            0.61, 0.56, 0.54, 0.54, 0.52, 0.52, 0.50, 0.58, 0.34, 0.10, 0.00, 0.32, 0.34, 0.42,
            0.56, 0.56, 0.54, 0.54, 0.56, 0.61,
        ],
        [
            0.81, 0.80, 0.80, 0.80, 0.80, 0.80, 0.62, 0.63, 0.60, 0.40, 0.30, 0.00, 0.03, 0.10,
            0.21, 0.21, 0.20, 0.20, 0.21, 0.24,
        ],
        [
            0.81, 0.80, 0.80, 0.80, 0.80, 0.80, 0.61, 0.63, 0.60, 0.40, 0.31, 0.03, 0.00, 0.08,
            0.21, 0.21, 0.20, 0.20, 0.21, 0.22,
        ],
        [
            0.80, 0.80, 1.00, 1.00, 0.80, 0.80, 0.60, 0.61, 0.61, 0.42, 0.34, 0.10, 0.08, 0.00,
            0.20, 0.20, 0.21, 0.21, 0.20, 0.20,
        ],
        [
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.80, 0.81, 0.80, 0.61, 0.52, 0.22, 0.21, 0.20,
            0.00, 0.00, 0.03, 0.03, 0.00, 0.05,
        ],
        [
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.80, 0.81, 0.80, 0.61, 0.52, 0.22, 0.21, 0.20,
            0.00, 0.00, 0.03, 0.03, 0.00, 0.05,
        ],
        [
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.81, 0.81, 0.80, 0.61, 0.51, 0.21, 0.20, 0.21,
            0.03, 0.03, 0.00, 0.00, 0.03, 0.08,
        ],
        [
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.81, 0.81, 0.80, 0.61, 0.51, 0.21, 0.20, 0.21,
            0.03, 0.03, 0.00, 0.00, 0.03, 0.08,
        ],
        [
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.80, 0.81, 0.80, 0.61, 0.52, 0.22, 0.21, 0.20,
            0.00, 0.00, 0.03, 0.03, 0.00, 0.05,
        ],
        [
            1.00, 1.00, 1.00, 1.00, 1.01, 1.01, 0.80, 0.80, 0.81, 0.62, 0.53, 0.24, 0.22, 0.20,
            0.05, 0.05, 0.08, 0.08, 0.05, 0.00,
        ],
    ]
});

/// Normalized Epstein distance in [0.0, 1.0]
pub fn compute(a: &AminoAcid, b: &AminoAcid) -> f64 {
    match (EPSTEIN_INDEX.get(a), EPSTEIN_INDEX.get(b)) {
        (Some(&i), Some(&j)) => EPSTEIN_MATRIX[i][j],
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::translation::{amino_acids::AminoAcid::*, distance::DistanceMetric};

    #[test]
    fn expected_values() {
        let cases = [
            (Arginine, Leucine, 1.0),
            (Arginine, Serine, 0.24),
            (Tryptophan, Cysteine, 0.18),
            (AsparticAcid, Asparagine, 0.0),
            (Lysine, GlutamicAcid, 0.0),
            (Isoleucine, Glycine, 0.51),
            (Phenylalanine, Tyrosine, 0.21),
            (Glutamine, Histidine, 0.2),
        ];
        let d = DistanceMetric::Epstein;
        let tolerance = 0.001;
        for (a, b, expected) in cases {
            assert!((d.compute(&a, &b) - expected).abs() < tolerance);
        }
    }
}
