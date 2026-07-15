//! Molecular fingerprints and Tanimoto similarity.

use openbabel_sys::ffi;

use crate::mol::Molecule;

/// A molecular fingerprint: a bit vector packed into 32-bit words.
///
/// ```no_run
/// use openbabel::{Molecule, Fingerprint};
/// let a = Molecule::parse("c1ccccc1", "smi").unwrap();
/// let b = Molecule::parse("c1ccccc1C", "smi").unwrap();
/// let fa = Fingerprint::compute(&a, "FP2").unwrap();
/// let fb = Fingerprint::compute(&b, "FP2").unwrap();
/// let sim = fa.tanimoto(&fb); // 0.0..=1.0
/// assert!(sim > 0.0);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fingerprint {
    words: Vec<u32>,
}

impl Fingerprint {
    /// Compute the fingerprint of `mol` using the named plugin.
    ///
    /// Common ids are `"FP2"` (path-based, the default), `"FP3"`, `"FP4"`, and
    /// `"MACCS"`. Returns `None` if OpenBabel has no such fingerprint plugin.
    pub fn compute(mol: &Molecule, kind: &str) -> Option<Fingerprint> {
        let words = crate::with_ob(|| ffi::fingerprint(mol.as_inner(), kind));
        if words.is_empty() {
            None
        } else {
            Some(Fingerprint { words })
        }
    }

    /// Tanimoto similarity coefficient with another fingerprint, in `0.0..=1.0`.
    pub fn tanimoto(&self, other: &Fingerprint) -> f64 {
        // OBFingerprint::Tanimoto is a pure static function over the vectors,
        // but we still serialize to keep all FFI access uniformly guarded.
        crate::with_ob(|| ffi::tanimoto(&self.words, &other.words))
    }

    /// The raw fingerprint words.
    pub fn as_words(&self) -> &[u32] {
        &self.words
    }

    /// Number of set bits (popcount) across the fingerprint.
    pub fn count_ones(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }
}
