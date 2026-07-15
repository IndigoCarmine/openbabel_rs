//! A borrowed view of a bond within a [`Molecule`](crate::Molecule).

use openbabel_sys::ffi;

/// A single bond, borrowed from its parent molecule.
#[derive(Clone, Copy)]
pub struct Bond<'mol> {
    pub(crate) mol: &'mol ffi::Molecule,
    /// OpenBabel's 0-based bond index.
    pub(crate) ob_idx: u32,
}

impl<'mol> Bond<'mol> {
    /// 0-based index of this bond within the molecule.
    pub fn index(&self) -> u32 {
        self.ob_idx
    }

    /// 0-based index of the atom at the start of this bond.
    pub fn begin_atom_index(&self) -> u32 {
        // OpenBabel returns 1-based atom indices; convert to 0-based.
        ffi::bond_begin_idx(self.mol, self.ob_idx).saturating_sub(1)
    }

    /// 0-based index of the atom at the end of this bond.
    pub fn end_atom_index(&self) -> u32 {
        ffi::bond_end_idx(self.mol, self.ob_idx).saturating_sub(1)
    }

    /// Bond order (1 = single, 2 = double, 3 = triple; 5 denotes aromatic in
    /// OpenBabel's convention).
    pub fn order(&self) -> u32 {
        ffi::bond_order(self.mol, self.ob_idx)
    }

    /// Whether this bond is part of an aromatic system.
    pub fn is_aromatic(&self) -> bool {
        ffi::bond_is_aromatic(self.mol, self.ob_idx)
    }

    /// Whether this bond is a member of any ring.
    pub fn is_in_ring(&self) -> bool {
        ffi::bond_is_in_ring(self.mol, self.ob_idx)
    }
}

impl std::fmt::Debug for Bond<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Bond")
            .field("index", &self.index())
            .field("begin", &self.begin_atom_index())
            .field("end", &self.end_atom_index())
            .field("order", &self.order())
            .finish()
    }
}
