//! A borrowed view of a ring within a [`Molecule`](crate::Molecule).
//!
//! Rings come from OpenBabel's SSSR (Smallest Set of Smallest Rings), perceived
//! on demand. Access them via [`Molecule::rings`](crate::Molecule::rings) /
//! [`ring`](crate::Molecule::ring); the count is
//! [`Molecule::num_rings`](crate::Molecule::num_rings).

use openbabel_sys::ffi;

/// A single ring of the SSSR, borrowed from its parent molecule.
#[derive(Clone, Copy)]
pub struct Ring<'mol> {
    mol: &'mol ffi::Molecule,
    /// 0-based ring index within the SSSR.
    idx: u32,
}

impl<'mol> Ring<'mol> {
    pub(crate) fn new(mol: &'mol ffi::Molecule, idx: u32) -> Self {
        Ring { mol, idx }
    }

    /// 0-based index of this ring within the SSSR.
    pub fn index(&self) -> u32 {
        self.idx
    }

    /// Number of atoms in the ring (its size, e.g. 6 for benzene).
    pub fn size(&self) -> u32 {
        crate::with_ob(|| ffi::ring_size(self.mol, self.idx))
    }

    /// 0-based indices of the atoms forming the ring, in cycle order.
    pub fn atom_indices(&self) -> Vec<u32> {
        crate::with_ob(|| ffi::ring_atom_indices(self.mol, self.idx))
    }

    /// Whether every atom in the ring is aromatic.
    pub fn is_aromatic(&self) -> bool {
        crate::with_ob(|| ffi::ring_is_aromatic(self.mol, self.idx))
    }
}

impl std::fmt::Debug for Ring<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Ring")
            .field("index", &self.index())
            .field("size", &self.size())
            .field("is_aromatic", &self.is_aromatic())
            .finish()
    }
}
