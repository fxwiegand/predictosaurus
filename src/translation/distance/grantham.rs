use crate::translation::amino_acids::AminoAcid;
use std::collections::HashMap;
use std::sync::LazyLock;

const GRAN_ORDER: [AminoAcid; 20] = [
    AminoAcid::Arginine,
    AminoAcid::Leucine,
    AminoAcid::Proline,
    AminoAcid::Threonine,
    AminoAcid::Alanine,
    AminoAcid::Valine,
    AminoAcid::Glycine,
    AminoAcid::Isoleucine,
    AminoAcid::Phenylalanine,
    AminoAcid::Tyrosine,
    AminoAcid::Cysteine,
    AminoAcid::Histidine,
    AminoAcid::Glutamine,
    AminoAcid::Asparagine,
    AminoAcid::Lysine,
    AminoAcid::AsparticAcid,
    AminoAcid::GlutamicAcid,
    AminoAcid::Methionine,
    AminoAcid::Tryptophan,
    AminoAcid::Serine,
];

pub static GRAN_INDEX: LazyLock<HashMap<AminoAcid, usize>> = LazyLock::new(|| {
    GRAN_ORDER
        .iter()
        .enumerate()
        .map(|(i, aa)| (*aa, i))
        .collect()
});

// Table containing Grantham distance matrix in the order
// Arg  Leu  Pro  Thr  Ala  Val  Gly  Ile  Phe  Tyr  Cys  His  Gln  Asn  Lys  Asp  Glu  Met  Trp  Ser
pub static GRAN_MATRIX: LazyLock<[[u16; 20]; 20]> = LazyLock::new(|| {
    [
        [
            0, 102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43, 86, 26, 96, 54, 91, 101, 110,
        ],
        [
            102, 0, 98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113, 153, 107, 172, 138, 15, 61, 145,
        ],
        [
            103, 98, 0, 38, 27, 68, 42, 95, 114, 110, 169, 77, 76, 91, 103, 108, 93, 87, 147, 74,
        ],
        [
            71, 92, 38, 0, 58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78, 85, 65, 81, 128, 58,
        ],
        [
            112, 96, 27, 58, 0, 64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107, 84, 148, 99,
        ],
        [
            96, 32, 68, 69, 64, 0, 109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88, 124,
        ],
        [
            125, 138, 42, 59, 60, 109, 0, 135, 153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184, 56,
        ],
        [
            97, 5, 95, 89, 94, 29, 135, 0, 21, 33, 198, 94, 109, 149, 102, 168, 134, 10, 61, 142,
        ],
        [
            97, 22, 114, 103, 113, 50, 153, 21, 0, 22, 205, 100, 116, 158, 102, 177, 140, 28, 40,
            155,
        ],
        [
            77, 36, 110, 92, 112, 55, 147, 33, 22, 0, 194, 83, 99, 143, 85, 160, 122, 36, 37, 144,
        ],
        [
            180, 198, 169, 149, 195, 192, 159, 198, 205, 194, 0, 174, 154, 139, 202, 154, 170, 196,
            215, 112,
        ],
        [
            29, 99, 77, 47, 86, 84, 98, 94, 100, 83, 174, 0, 24, 68, 32, 81, 40, 87, 115, 89,
        ],
        [
            43, 113, 76, 42, 91, 96, 87, 109, 116, 99, 154, 24, 0, 46, 53, 61, 29, 101, 130, 68,
        ],
        [
            86, 153, 91, 65, 111, 133, 80, 149, 158, 143, 139, 68, 46, 0, 94, 23, 42, 142, 174, 46,
        ],
        [
            26, 107, 103, 78, 106, 97, 127, 102, 102, 85, 202, 32, 53, 94, 0, 101, 56, 95, 110, 121,
        ],
        [
            96, 172, 108, 85, 126, 152, 94, 168, 177, 160, 154, 81, 61, 23, 101, 0, 45, 160, 181,
            65,
        ],
        [
            54, 138, 93, 65, 107, 121, 98, 134, 140, 122, 170, 40, 29, 42, 56, 45, 0, 126, 152, 80,
        ],
        [
            91, 15, 87, 81, 84, 21, 127, 10, 28, 36, 196, 87, 101, 142, 95, 160, 126, 0, 67, 135,
        ],
        [
            101, 61, 147, 128, 148, 88, 184, 61, 40, 37, 215, 115, 130, 174, 110, 181, 152, 67, 0,
            177,
        ],
        [
            110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65, 80, 135, 177, 0,
        ],
    ]
});

/// Normalized Grantham distance in [0.0, 1.0]
pub fn compute(a: &AminoAcid, b: &AminoAcid) -> f64 {
    match (GRAN_INDEX.get(a), GRAN_INDEX.get(b)) {
        (Some(&i), Some(&j)) => GRAN_MATRIX[i][j] as f64 / 215.0,
        _ => 1.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::translation::{amino_acids::AminoAcid::*, distance::DistanceMetric};

    #[test]
    fn symmetry() {
        for i in 0..20 {
            for j in 0..20 {
                assert_eq!(GRAN_MATRIX[i][j], GRAN_MATRIX[j][i]);
            }
        }
    }

    #[test]
    fn expected_values() {
        let cases = [
            (Arginine, Leucine, 102),
            (Arginine, Serine, 110),
            (Tryptophan, Cysteine, 215),
            (AsparticAcid, Asparagine, 23),
            (Lysine, GlutamicAcid, 56),
            (Isoleucine, Glycine, 135),
            (Phenylalanine, Tyrosine, 22),
            (Glutamine, Histidine, 24),
        ];
        let d = DistanceMetric::Grantham;
        let tolerance = 0.001;
        for (a, b, expected) in cases {
            assert!((d.compute(&a, &b) - (expected as f64 / 215.0)).abs() < tolerance);
        }
    }
}
