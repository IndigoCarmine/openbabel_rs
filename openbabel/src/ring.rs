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

    /// The ring's type name as classified by OpenBabel's ring typer (e.g.
    /// `"benzene"`). Empty if ring-type perception has not run for the molecule.
    pub fn ring_type(&self) -> String {
        crate::with_ob(|| ffi::ring_type(self.mol, self.idx))
    }

    /// The 0-based index of the ring's root atom — the heteroatom OpenBabel
    /// anchors the ring on (O for furan, N for pyrrole, …). `None` for an
    /// all-carbon ring, which has no distinguished root.
    pub fn root_atom(&self) -> Option<u32> {
        let one_based = crate::with_ob(|| ffi::ring_root_atom(self.mol, self.idx));
        (one_based != 0).then(|| one_based - 1)
    }

    /// Whether `atom` is a member of this ring.
    pub fn contains_atom(&self, atom: &crate::Atom<'mol>) -> bool {
        crate::with_ob(|| ffi::ring_contains_atom(self.mol, self.idx, atom.ob_idx))
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
