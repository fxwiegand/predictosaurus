use crate::translation::amino_acids::AminoAcid;
use std::collections::HashMap;
use std::sync::LazyLock;

const MIYATA_ORDER: [AminoAcid; 20] = [
    AminoAcid::Cysteine,
    AminoAcid::Proline,
    AminoAcid::Alanine,
    AminoAcid::Glycine,
    AminoAcid::Serine,
    AminoAcid::Threonine,
    AminoAcid::Glutamine,
    AminoAcid::GlutamicAcid,
    AminoAcid::Asparagine,
    AminoAcid::AsparticAcid,
    AminoAcid::Histidine,
    AminoAcid::Lysine,
    AminoAcid::Arginine,
    AminoAcid::Valine,
    AminoAcid::Leucine,
    AminoAcid::Isoleucine,
    AminoAcid::Methionine,
    AminoAcid::Phenylalanine,
    AminoAcid::Tyrosine,
    AminoAcid::Tryptophan,
];

pub static MIYATA_INDEX: LazyLock<HashMap<AminoAcid, usize>> = LazyLock::new(|| {
    MIYATA_ORDER
        .iter()
        .enumerate()
        .map(|(i, aa)| (*aa, i))
        .collect()
});

// Table containing Miyata distance matrix
const MIYATA_MATRIX: [[f64; 20]; 20] = [
    [
        0.0, 1.33, 1.39, 2.22, 2.84, 1.45, 2.48, 3.26, 2.83, 3.48, 2.56, 3.27, 3.06, 0.86, 1.65,
        1.63, 1.46, 2.24, 2.38, 3.34,
    ],
    [
        1.33, 0.0, 0.06, 0.97, 0.56, 0.87, 1.92, 2.48, 1.80, 2.40, 2.15, 2.94, 2.90, 1.79, 2.70,
        2.62, 2.36, 3.17, 3.12, 4.17,
    ],
    [
        1.39, 0.06, 0.0, 0.91, 0.51, 0.90, 1.92, 2.46, 1.78, 2.37, 2.17, 2.96, 2.92, 1.85, 2.76,
        2.69, 2.42, 3.23, 3.18, 4.23,
    ],
    [
        2.22, 0.97, 0.91, 0.0, 0.85, 1.70, 2.48, 2.78, 1.96, 2.37, 2.78, 3.54, 3.58, 2.76, 3.67,
        3.60, 3.34, 4.14, 4.08, 5.13,
    ],
    [
        2.84, 0.56, 0.51, 0.85, 0.0, 0.89, 1.65, 2.06, 1.31, 1.87, 1.94, 2.71, 2.74, 2.15, 3.04,
        2.95, 2.67, 3.45, 3.33, 4.38,
    ],
    [
        1.45, 0.87, 0.90, 1.70, 0.89, 0.0, 1.12, 1.83, 1.40, 2.05, 1.32, 2.10, 2.03, 1.42, 2.25,
        2.14, 1.86, 2.60, 2.45, 3.50,
    ],
    [
        2.48, 1.92, 1.92, 2.48, 1.65, 1.12, 0.0, 0.84, 0.99, 1.47, 0.32, 1.06, 1.13, 2.13, 2.70,
        2.57, 2.30, 2.81, 2.48, 3.42,
    ],
    [
        3.26, 2.48, 2.46, 2.78, 2.06, 1.83, 0.84, 0.0, 0.85, 0.90, 0.96, 1.14, 1.45, 2.97, 3.53,
        3.39, 3.13, 3.59, 3.22, 4.08,
    ],
    [
        2.83, 1.80, 1.78, 1.96, 1.31, 1.40, 0.99, 0.85, 0.0, 0.65, 1.29, 1.84, 2.04, 2.76, 3.49,
        3.37, 3.08, 3.70, 3.42, 4.39,
    ],
    [
        3.48, 2.40, 2.37, 2.37, 1.87, 2.05, 1.47, 0.90, 0.65, 0.0, 1.72, 2.05, 2.34, 3.40, 4.10,
        3.98, 3.69, 4.27, 3.95, 4.88,
    ],
    [
        2.56, 2.15, 2.17, 2.78, 1.94, 1.32, 0.32, 0.96, 1.29, 1.72, 0.0, 0.79, 0.82, 2.11, 2.59,
        2.45, 2.19, 2.63, 2.27, 3.16,
    ],
    [
        3.27, 2.94, 2.96, 3.54, 2.71, 2.10, 1.06, 1.14, 1.84, 2.05, 0.79, 0.0, 0.40, 2.70, 2.98,
        2.84, 2.63, 2.85, 2.42, 3.11,
    ],
    [
        3.06, 2.90, 2.92, 3.58, 2.74, 2.03, 1.13, 1.45, 2.04, 2.34, 0.82, 0.40, 0.0, 2.43, 2.62,
        2.49, 2.29, 2.47, 2.02, 2.72,
    ],
    [
        0.86, 1.79, 1.85, 2.76, 2.15, 1.42, 2.13, 2.97, 2.76, 3.40, 2.11, 2.70, 2.43, 0.0, 0.91,
        0.85, 0.62, 1.43, 1.52, 2.51,
    ],
    [
        1.65, 2.70, 2.76, 3.67, 3.04, 2.25, 2.70, 3.53, 3.49, 4.10, 2.59, 2.98, 2.62, 0.91, 0.0,
        0.14, 0.41, 0.63, 0.94, 1.73,
    ],
    [
        1.63, 2.62, 2.69, 3.60, 2.95, 2.14, 2.57, 3.39, 3.37, 3.98, 2.45, 2.84, 2.49, 0.85, 0.14,
        0.0, 0.29, 0.61, 0.86, 1.72,
    ],
    [
        1.46, 2.36, 2.42, 3.34, 2.67, 1.86, 2.30, 3.13, 3.08, 3.69, 2.19, 2.63, 2.29, 0.62, 0.41,
        0.29, 0.0, 0.82, 0.93, 1.89,
    ],
    [
        2.24, 3.17, 3.23, 4.14, 3.45, 2.60, 2.81, 3.59, 3.70, 4.27, 2.63, 2.85, 2.47, 1.43, 0.63,
        0.61, 0.82, 0.0, 0.48, 1.11,
    ],
    [
        2.38, 3.12, 3.18, 4.08, 3.33, 2.45, 2.48, 3.22, 3.42, 3.95, 2.27, 2.42, 2.02, 1.52, 0.94,
        0.86, 0.93, 0.48, 0.0, 1.06,
    ],
    [
        3.34, 4.17, 4.23, 5.13, 4.38, 3.50, 3.42, 4.08, 4.39, 4.88, 3.16, 3.11, 2.72, 2.51, 1.73,
        1.72, 1.89, 1.11, 1.06, 0.0,
    ],
];

/// Computes the normalized Miyata distance between two amino acids.
///
/// Returns a value in the range [0.0, 1.0], where 0.0 indicates identical amino acids or both are stop codons, and 1.0 indicates maximal distance or one is a stop codon and the other is not. For standard amino acids, the distance is based on the Miyata matrix and normalized by the maximum possible value.
///
/// # Examples
///
/// ```
/// use crate::AminoAcid;
/// use crate::compute;
///
/// let ala = AminoAcid::Ala;
/// let gly = AminoAcid::Gly;
/// let stop = AminoAcid::Stop;
///
/// let d1 = compute(&ala, &gly);
/// assert!(d1 >= 0.0 && d1 <= 1.0);
///
/// let d2 = compute(&stop, &stop);
/// assert_eq!(d2, 0.0);
///
/// let d3 = compute(&ala, &stop);
/// assert_eq!(d3, 1.0);
/// ```
pub fn compute(a: &AminoAcid, b: &AminoAcid) -> f64 {
    if a.is_stop() || b.is_stop() {
        return if a == b { 0.0 } else { 1.0 };
    }
    match (MIYATA_INDEX.get(a), MIYATA_INDEX.get(b)) {
        (Some(&i), Some(&j)) => MIYATA_MATRIX[i][j] / 5.13,
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::translation::{amino_acids::AminoAcid::*, distance::DistanceMetric};

    #[test]
    fn test_stop() {
        assert_eq!(compute(&Stop, &Stop), 0.0);
        assert_eq!(compute(&Stop, &Arginine), 1.0);
        assert_eq!(compute(&Arginine, &Stop), 1.0);
    }

    #[test]
    fn miyata_symmetry() {
        for i in 0..20 {
            for j in 0..20 {
                assert_eq!(MIYATA_MATRIX[i][j], MIYATA_MATRIX[j][i]);
            }
        }
    }

    #[test]
    fn miyata_expected_values() {
        let cases = [
            (Arginine, Leucine, 2.62),
            (Arginine, Serine, 2.74),
            (Tryptophan, Cysteine, 3.34),
            (AsparticAcid, Asparagine, 0.65),
            (Lysine, AsparticAcid, 2.05),
            (Lysine, GlutamicAcid, 1.14),
            (Isoleucine, Glycine, 3.60),
            (Phenylalanine, Tyrosine, 0.48),
            (GlutamicAcid, Histidine, 0.96),
            (Glutamine, Histidine, 0.32),
        ];
        let d = DistanceMetric::Miyata;
        let tolerance = 0.001;
        for (a, b, expected) in cases {
            assert!((d.compute(&a, &b) - (expected as f64 / 5.13)).abs() < tolerance);
        }
    }
}
