pub mod epstein;
pub mod grantham;
pub mod miyata;
pub mod sneath;
use crate::translation::amino_acids::AminoAcid;

#[derive(Clone, Copy, Default, Debug)]
pub enum DistanceMetric {
    #[default]
    Grantham,
    Sneath,
    Epstein,
    Miyata,
}

impl DistanceMetric {
    /// Computes the distance between two amino acids using the selected metric.
    ///
    /// The distance metric is determined by the `DistanceMetric` variant. The result is a floating-point value representing the computed distance according to the chosen metric.
    ///
    /// # Examples
    ///
    /// ```
    /// use translation::distance::{DistanceMetric, AminoAcid};
    ///
    /// let metric = DistanceMetric::Grantham;
    /// let a = AminoAcid::Alanine;
    /// let b = AminoAcid::Valine;
    /// let distance = metric.compute(&a, &b);
    /// assert!(distance >= 0.0);
    /// ```
    pub fn compute(&self, a: &AminoAcid, b: &AminoAcid) -> f64 {
        match self {
            DistanceMetric::Grantham => grantham::compute(a, b),
            DistanceMetric::Sneath => sneath::compute(a, b),
            DistanceMetric::Epstein => epstein::compute(a, b),
            DistanceMetric::Miyata => miyata::compute(a, b),
        }
    }
}
