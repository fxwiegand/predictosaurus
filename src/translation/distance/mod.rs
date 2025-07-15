pub mod epstein;
pub mod grantham;
pub mod miyata;
pub mod sneath;
use crate::translation::amino_acids::AminoAcid;

#[derive(Clone, Copy)]
pub enum DistanceMetric {
    Grantham,
    Sneath,
    Epstein,
    Miyata,
}

impl DistanceMetric {
    pub fn compute(&self, a: &AminoAcid, b: &AminoAcid) -> f64 {
        match self {
            DistanceMetric::Grantham => grantham::compute(a, b),
            DistanceMetric::Sneath => sneath::compute(a, b),
            DistanceMetric::Epstein => epstein::compute(a, b),
            DistanceMetric::Miyata => miyata::compute(a, b),
        }
    }
}
