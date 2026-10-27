use anyhow::{Result, bail};

use crate::{Component, ILLUMINA_QUALITY_OFFSET, ResolvedGeometry, geometry::ReadMate};

/// An extraction struct for the UMI of a read
#[derive(Clone, Copy)]
pub struct UmiMapper {
    pos: usize,
    len: usize,
    mate: ReadMate,
    /// Minimum Phred quality required at every UMI base (0 disables the check)
    min_qual: u8,
}
impl UmiMapper {
    pub fn new(pos: usize, len: usize, mate: ReadMate, min_qual: u8) -> Self {
        Self {
            pos,
            len,
            mate,
            min_qual,
        }
    }

    pub fn resolve(geometry: &ResolvedGeometry, min_qual: u8) -> Result<Self> {
        let Some(region) = geometry.get(Component::Umi) else {
            bail!("geometry missing [umi]")
        };
        let Some(len) = region.length else {
            bail!("length missing [umi]")
        };
        Ok(Self::new(region.offset, len, region.mate, min_qual))
    }

    /// Returns the mate of the read
    #[inline]
    pub fn mate(&self) -> ReadMate {
        self.mate
    }

    /// Determines the linear range of the UMI extraction
    #[inline]
    pub fn range(&self) -> std::ops::Range<usize> {
        self.pos..self.pos + self.len
    }

    /// Extracts the UMI from the given sequence
    #[inline]
    pub fn extract_umi<'a>(&self, seq: &'a [u8]) -> Option<&'a [u8]> {
        seq.get(self.range())
    }

    /// Extracts the UMI from the given sequence as a 2-bit encoded u64 (bitnuc)
    #[inline]
    pub fn extract_2bit_umi(&self, seq: &[u8]) -> Option<Result<u64>> {
        self.extract_umi(seq)
            .map(|umi| bitnuc::as_2bit(umi).map_err(anyhow::Error::new))
    }

    /// Validates that all quality scores are above the required threshold
    #[inline]
    pub fn passes_quality_threshold(&self, qual: &[u8]) -> bool {
        if self.min_qual == 0 {
            return true;
        }
        qual.get(self.range()) // pull the range
            .is_none_or(|sub_qual| {
                // iter over all quality scores and ensure all are above the threshold
                sub_qual
                    .iter()
                    .all(|q| q.saturating_sub(ILLUMINA_QUALITY_OFFSET) >= self.min_qual)
            }) // default to true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Phred+33 encode a slice of quality scores
    fn phred(scores: &[u8]) -> Vec<u8> {
        scores.iter().map(|q| q + ILLUMINA_QUALITY_OFFSET).collect()
    }

    #[test]
    fn test_quality_threshold() {
        let umi = UmiMapper::new(2, 4, ReadMate::R1, 20);
        // bases outside the UMI range (0..2, 6..) are ignored
        assert!(umi.passes_quality_threshold(&phred(&[0, 0, 20, 30, 40, 20, 0])));
        assert!(!umi.passes_quality_threshold(&phred(&[40, 40, 20, 19, 40, 40, 40])));
    }

    #[test]
    fn test_quality_threshold_disabled() {
        let umi = UmiMapper::new(0, 4, ReadMate::R1, 0);
        assert!(umi.passes_quality_threshold(&phred(&[0, 0, 0, 0])));
    }

    #[test]
    fn test_quality_threshold_missing_or_invalid() {
        let umi = UmiMapper::new(0, 4, ReadMate::R1, crate::UMI_MIN_QUALITY);
        // missing / truncated qualities default to passing
        assert!(umi.passes_quality_threshold(&[]));
        assert!(umi.passes_quality_threshold(&phred(&[40, 40])));
        // bytes below the Phred+33 offset must fail, not wrap around
        assert!(!umi.passes_quality_threshold(&[0, 0, 0, 0]));
    }
}
