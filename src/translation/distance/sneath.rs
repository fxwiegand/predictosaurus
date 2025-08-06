use crate::translation::amino_acids::AminoAcid;
use std::collections::HashMap;
use std::sync::LazyLock;

const SNEATH_ORDER: [AminoAcid; 20] = [
    AminoAcid::Leucine,
    AminoAcid::Isoleucine,
    AminoAcid::Valine,
    AminoAcid::Glycine,
    AminoAcid::Alanine,
    AminoAcid::Proline,
    AminoAcid::Glutamine,
    AminoAcid::Asparagine,
    AminoAcid::Methionine,
    AminoAcid::Threonine,
    AminoAcid::Serine,
    AminoAcid::Cysteine,
    AminoAcid::GlutamicAcid,
    AminoAcid::AsparticAcid,
    AminoAcid::Lysine,
    AminoAcid::Arginine,
    AminoAcid::Tyrosine,
    AminoAcid::Phenylalanine,
    AminoAcid::Tryptophan,
    AminoAcid::Histidine,
];

pub static SNEATH_INDEX: LazyLock<HashMap<AminoAcid, usize>> = LazyLock::new(|| {
    SNEATH_ORDER
        .iter()
        .enumerate()
        .map(|(i, aa)| (*aa, i))
        .collect()
});

// Table containing Sneath distance matrix in the order
// Leu  Ile  Val  Gly  Ala  Pro  Gln  Asn  Met  Thr  Ser  Cys  Glu  Asp  Lys  Arg  Tyr  Phe  Trp  His
const SNEATH_MATRIX: [[u8; 20]; 20] = [
    [
        0, 5, 9, 24, 15, 23, 22, 20, 20, 23, 23, 24, 30, 25, 23, 33, 30, 19, 30, 25,
    ],
    [
        5, 0, 7, 25, 17, 24, 24, 23, 22, 21, 25, 26, 31, 28, 24, 34, 34, 22, 34, 28,
    ],
    [
        9, 7, 0, 19, 12, 20, 25, 23, 23, 17, 20, 21, 31, 28, 26, 36, 36, 26, 37, 31,
    ],
    [
        24, 25, 19, 0, 9, 17, 32, 26, 34, 20, 19, 21, 37, 33, 31, 43, 36, 29, 39, 34,
    ],
    [
        15, 17, 12, 9, 0, 16, 26, 25, 25, 20, 16, 13, 34, 30, 26, 37, 34, 26, 36, 29,
    ],
    [
        23, 24, 20, 17, 16, 0, 33, 31, 31, 25, 24, 25, 43, 40, 31, 43, 37, 27, 37, 36,
    ],
    [
        22, 24, 25, 32, 26, 33, 0, 10, 13, 24, 21, 22, 14, 22, 21, 23, 29, 24, 31, 27,
    ],
    [
        20, 23, 23, 26, 25, 31, 10, 0, 21, 19, 15, 19, 19, 14, 27, 31, 28, 24, 32, 24,
    ],
    [
        20, 22, 23, 34, 25, 31, 13, 21, 0, 25, 22, 17, 26, 31, 24, 28, 32, 24, 31, 30,
    ],
    [
        23, 21, 17, 20, 20, 25, 24, 19, 25, 0, 12, 19, 34, 29, 34, 38, 32, 28, 38, 34,
    ],
    [
        23, 25, 20, 19, 16, 24, 21, 15, 22, 12, 0, 13, 29, 25, 31, 37, 29, 25, 35, 28,
    ],
    [
        24, 26, 21, 21, 13, 25, 22, 19, 17, 19, 13, 0, 33, 28, 32, 36, 34, 29, 37, 31,
    ],
    [
        30, 31, 31, 37, 34, 43, 14, 19, 26, 34, 29, 33, 0, 7, 26, 31, 34, 35, 43, 27,
    ],
    [
        25, 28, 28, 33, 30, 40, 22, 14, 31, 29, 25, 28, 7, 0, 34, 39, 34, 35, 45, 35,
    ],
    [
        23, 24, 26, 31, 26, 31, 21, 27, 24, 34, 31, 32, 26, 34, 0, 14, 34, 28, 34, 27,
    ],
    [
        33, 34, 36, 43, 37, 43, 23, 31, 28, 38, 37, 36, 31, 39, 14, 0, 36, 34, 36, 31,
    ],
    [
        30, 34, 36, 36, 34, 37, 29, 28, 32, 32, 29, 34, 34, 34, 34, 36, 0, 13, 21, 23,
    ],
    [
        19, 22, 26, 29, 26, 27, 24, 24, 24, 28, 25, 29, 35, 35, 28, 34, 13, 0, 13, 18,
    ],
    [
        30, 34, 37, 39, 36, 37, 31, 32, 31, 38, 35, 37, 43, 45, 34, 36, 21, 13, 0, 25,
    ],
    [
        25, 28, 31, 34, 29, 36, 27, 24, 30, 34, 28, 31, 27, 35, 27, 31, 23, 18, 25, 0,
    ],
];

/// Computes the normalized Sneath distance between two amino acids.
///
/// Returns a value in the range [0.0, 1.0], where 0.0 indicates identical amino acids (or both are stop codons),
/// and 1.0 indicates maximal dissimilarity (or one is a stop codon and the other is not).
/// For standard amino acids, the distance is derived from the Sneath matrix and normalized by dividing by 45.0.
///
/// # Examples
///
/// ```
/// use translation::amino_acids::AminoAcid;
/// use translation::distance::sneath::compute;
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
    match (SNEATH_INDEX.get(a), SNEATH_INDEX.get(b)) {
        (Some(&i), Some(&j)) => SNEATH_MATRIX[i][j] as f64 / 45.0,
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
    fn sneath_symmetry() {
        for i in 0..20 {
            for j in 0..20 {
                assert_eq!(SNEATH_MATRIX[i][j], SNEATH_MATRIX[j][i]);
            }
        }
    }

    #[test]
    fn sneath_expected_values() {
        let cases = [
            (Arginine, Leucine, 33),
            (Arginine, Serine, 37),
            (Tryptophan, Cysteine, 37),
            (AsparticAcid, Asparagine, 14),
            (Lysine, AsparticAcid, 34),
            (Lysine, GlutamicAcid, 26),
            (Isoleucine, Glycine, 25),
            (Phenylalanine, Tyrosine, 13),
            (Glutamine, Histidine, 27),
        ];
        let d = DistanceMetric::Sneath;
        let tolerance = 0.001;
        for (a, b, expected) in cases {
            assert!((d.compute(&a, &b) - (expected as f64 / 45.0)).abs() < tolerance);
        }
    }
}
