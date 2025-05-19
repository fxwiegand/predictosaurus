/// Tuning constants for the score components
const SNP_WEIGHT: f64 = 0.1; // Weight per SNP
const FS_WEIGHT: f64 = 1.0; // Weight for frameshift fraction
const STOP_WEIGHT: f64 = 1.0; // Weight for stop-gained fraction

/// Breakdown of effects for one haplotype path
#[derive(Debug, Clone)]
pub struct EffectScore {
    /// Number of coding SNPs on this haplotype
    pub num_snps: usize,
    /// Total fractional CDS length affected by frameshifts (0.0–1.0)
    pub fs_fraction: f64,
    /// Stop-gained penalty
    pub stop_fraction: f64,
}

impl EffectScore {
    /// Create a new empty score
    pub fn new() -> Self {
        EffectScore {
            num_snps: 0,
            fs_fraction: 0.0,
            stop_fraction: 0.0,
        }
    }

    /// Compute the raw combined score using tuning constants
    pub fn raw(&self) -> f64 {
        FS_WEIGHT * self.fs_fraction
            + SNP_WEIGHT * (self.num_snps as f64)
            + STOP_WEIGHT * self.stop_fraction
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

    #[test]
    fn test_raw_score() {
        let score = EffectScore {
            num_snps: 3,
            fs_fraction: 0.4,
            stop_fraction: 0.2,
        };
        let expected = FS_WEIGHT * 0.4 + SNP_WEIGHT * 3.0 + STOP_WEIGHT * 0.2;
        assert!((score.raw() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_normalized_score() {
        let score = EffectScore {
            num_snps: 2,
            fs_fraction: 0.1,
            stop_fraction: 0.3,
        };
        let raw = score.raw();
        let expected = raw / (1.0 + raw);
        assert!((score.normalized() - expected).abs() < 1e-6);
    }
}
