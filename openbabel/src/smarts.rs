//! SMARTS substructure matching.

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::error::Error;
use crate::mol::Molecule;

/// A compiled SMARTS query pattern.
///
/// ```no_run
/// use openbabel::{Molecule, SmartsPattern};
/// let mol = Molecule::parse("c1ccccc1O", "smi").unwrap(); // phenol
/// let oh = SmartsPattern::new("[OX2H]").unwrap();          // hydroxyl
/// assert!(oh.matches(&mol));
/// ```
pub struct SmartsPattern {
    inner: UniquePtr<ffi::Smarts>,
}

impl SmartsPattern {
    /// Compile a SMARTS pattern.
    ///
    /// Returns [`Error::InvalidSmarts`] if the pattern is malformed.
    pub fn new(pattern: &str) -> Result<Self, Error> {
        let inner = crate::with_ob(|| ffi::smarts_new(pattern));
        if inner.is_null() {
            Err(Error::InvalidSmarts {
                pattern: pattern.to_string(),
            })
        } else {
            Ok(SmartsPattern { inner })
        }
    }

    /// Number of atoms in the pattern (the length of each match).
    pub fn atom_count(&self) -> u32 {
        crate::with_ob(|| ffi::smarts_atom_count(self.as_inner()))
    }

    /// Whether the pattern matches `mol` at least once.
    pub fn matches(&self, mol: &Molecule) -> bool {
        crate::with_ob(|| ffi::smarts_matches(self.as_inner(), mol.as_inner()))
    }

    /// All unique matches, each as a `Vec` of 0-based atom indices into `mol`.
    pub fn match_indices(&self, mol: &Molecule) -> Vec<Vec<u32>> {
        crate::with_ob(|| {
            let width = ffi::smarts_atom_count(self.as_inner()) as usize;
            if width == 0 {
                return Vec::new();
            }
            let flat = ffi::smarts_match_atoms(self.as_inner(), mol.as_inner());
            flat.chunks(width)
                .map(|chunk| chunk.iter().map(|&i| i.saturating_sub(1)).collect())
                .collect()
        })
    }

    /// Number of unique matches of the pattern in `mol`.
    pub fn num_matches(&self, mol: &Molecule) -> usize {
        crate::with_ob(|| {
            let width = ffi::smarts_atom_count(self.as_inner()) as usize;
            if width == 0 {
                return 0;
            }
            ffi::smarts_match_atoms(self.as_inner(), mol.as_inner()).len() / width
        })
    }

    fn as_inner(&self) -> &ffi::Smarts {
        self.inner.as_ref().expect("SmartsPattern is never null")
    }
}
