use crate::graph::transcript::Transcript;
use crate::utils::fasta::reverse_complement;
use bio::bio_types::strand::Strand;

pub(crate) fn hgvsc(
    transcript: &Transcript,
    position: u64,
    reference_allele: &str,
    alternative_allele: &str,
) -> Option<String> {
    let forward_pos = transcript.position_in_transcript(position as usize).ok()? + 1; // 1-based
    let cdna_pos = if transcript.strand == Strand::Reverse {
        transcript.length() - forward_pos + 1
    } else {
        forward_pos
    };

    let (ref_allele, alt_allele) = if transcript.strand == Strand::Reverse {
        (
            String::from_utf8(reverse_complement(reference_allele.as_bytes())).unwrap(),
            String::from_utf8(reverse_complement(alternative_allele.as_bytes())).unwrap(),
        )
    } else {
        (reference_allele.to_string(), alternative_allele.to_string())
    };

    let notation = match (ref_allele.is_empty(), alt_allele.is_empty()) {
        (true, false) => format!("{cdna_pos}_{}ins{alt_allele}", cdna_pos + 1),
        (false, true) => match ref_allele.len() {
            1 => format!("{cdna_pos}del"),
            n => format!("{cdna_pos}_{}del", cdna_pos + n - 1),
        },
        (false, false) => match ref_allele.len() {
            1 => format!("{cdna_pos}{ref_allele}>{alt_allele}"),
            n => format!("{cdna_pos}_{}delins{alt_allele}", cdna_pos + n - 1),
        },
        (true, true) => return None,
    };

    Some(notation)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{graph::transcript::Transcript, Cds};
    use bio::bio_types::strand::Strand;

    fn transcript(cdss: Vec<Cds>, strand: Strand) -> Transcript {
        Transcript {
            feature: "ENST1".into(),
            target: "chr1".into(),
            strand,
            coding_sequences: cdss,
        }
    }

    // Plus strand
    #[test]
    fn snv_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                114,
                "C",
                "G"
            ),
            Some("15C>G".into())
        );
    }
    #[test]
    fn snv_first_base_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                100,
                "A",
                "T"
            ),
            Some("1A>T".into())
        );
    }
    #[test]
    fn del_single_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "A",
                ""
            ),
            Some("6del".into())
        );
    }
    #[test]
    fn del_multi_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "ATG",
                ""
            ),
            Some("6_8del".into())
        );
    }
    #[test]
    fn ins_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "",
                "ATTCAC"
            ),
            Some("6_7insATTCAC".into())
        );
    }

    #[test]
    fn multi_exon_forward() {
        let t = transcript(
            vec![Cds::new(100, 109, 0), Cds::new(200, 299, 0)],
            Strand::Forward,
        );
        assert_eq!(hgvsc(&t, 200, "C", "T"), Some("11C>T".into()));
        assert_eq!(hgvsc(&t, 204, "G", "A"), Some("15G>A".into()));
    }

    // Minus strand
    #[test]
    fn snv_reverse_allele_flip() {
        // offset from 3' end: 200-114=86 => cDNA 87; C>G => G>C
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Reverse),
                114,
                "C",
                "G"
            ),
            Some("87G>C".into())
        );
    }
    #[test]
    fn snv_last_base_reverse() {
        // pos 200 = cDNA 1 on minus; A>T => T>A
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Reverse),
                200,
                "A",
                "T"
            ),
            Some("1T>A".into())
        );
    }
    #[test]
    fn ins_reverse() {
        // ATG rev_comp => CAT
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Reverse),
                195,
                "",
                "ATG"
            ),
            Some("6_7insCAT".into())
        );
    }

    // Edge cases
    #[test]
    fn outside_cds() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                50,
                "A",
                "T"
            ),
            None
        );
    }
    #[test]
    fn empty_alleles() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "",
                ""
            ),
            None
        );
    }

    #[test]
    fn rc_via_fasta_util() {
        assert_eq!(
            String::from_utf8(reverse_complement(b"ATCG")).unwrap(),
            "CGAT"
        );
    }

    #[test]
    fn delins_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "ATG",
                "TC"
            ),
            Some("6_8delinsTC".into())
        );
    }
}
