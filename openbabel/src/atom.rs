//! A borrowed view of an atom within a [`Molecule`](crate::Molecule).

use openbabel_sys::ffi;

/// A single atom, borrowed from its parent molecule.
///
/// An `Atom` is a lightweight handle (a reference to the molecule plus an
/// index); it cannot outlive the molecule it came from.
#[derive(Clone, Copy)]
pub struct Atom<'mol> {
    pub(crate) mol: &'mol ffi::Molecule,
    /// OpenBabel's 1-based atom index.
    pub(crate) ob_idx: u32,
}

impl<'mol> Atom<'mol> {
    /// 0-based index of this atom within the molecule.
    pub fn index(&self) -> u32 {
        self.ob_idx - 1
    }

    /// Atomic number (e.g. 6 for carbon, 8 for oxygen).
    pub fn atomic_number(&self) -> u32 {
        crate::with_ob(|| ffi::atom_atomic_num(self.mol, self.ob_idx))
    }

    /// Cartesian coordinates `(x, y, z)`. Zero unless the molecule has a
    /// conformer / 3D structure.
    pub fn coords(&self) -> (f64, f64, f64) {
        crate::with_ob(|| {
            (
                ffi::atom_x(self.mol, self.ob_idx),
                ffi::atom_y(self.mol, self.ob_idx),
                ffi::atom_z(self.mol, self.ob_idx),
            )
        })
    }

    /// Formal charge on this atom.
    pub fn formal_charge(&self) -> i32 {
        crate::with_ob(|| ffi::atom_formal_charge(self.mol, self.ob_idx))
    }

    /// Partial (fractional) atomic charge.
    ///
    /// Reflects the last model assigned via
    /// [`Molecule::compute_charges`](crate::Molecule::compute_charges); if none
    /// was, OpenBabel computes Gasteiger charges on demand.
    pub fn partial_charge(&self) -> f64 {
        crate::with_ob(|| ffi::atom_partial_charge(self.mol, self.ob_idx))
    }

    /// Whether this atom is part of an aromatic system.
    pub fn is_aromatic(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_aromatic(self.mol, self.ob_idx))
    }

    /// Whether this atom is a member of any ring.
    pub fn is_in_ring(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_in_ring(self.mol, self.ob_idx))
    }

    /// Number of explicit connections (bonds) to this atom.
    pub fn degree(&self) -> u32 {
        crate::with_ob(|| ffi::atom_degree(self.mol, self.ob_idx))
    }

    /// Total valence, counting implicit hydrogens.
    pub fn total_valence(&self) -> u32 {
        crate::with_ob(|| ffi::atom_total_valence(self.mol, self.ob_idx))
    }

    /// Number of implicit (not explicitly present) hydrogens on this atom.
    pub fn implicit_hydrogens(&self) -> u32 {
        crate::with_ob(|| ffi::atom_implicit_h_count(self.mol, self.ob_idx))
    }

    /// Hybridization: 1 = sp, 2 = sp², 3 = sp³, … (0 if unassigned).
    pub fn hybridization(&self) -> u32 {
        crate::with_ob(|| ffi::atom_hybridization(self.mol, self.ob_idx))
    }

    /// Whether this atom can donate a hydrogen bond.
    pub fn is_hbond_donor(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_hbond_donor(self.mol, self.ob_idx))
    }

    /// Whether this atom can accept a hydrogen bond.
    pub fn is_hbond_acceptor(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_hbond_acceptor(self.mol, self.ob_idx))
    }
}

impl std::fmt::Debug for Atom<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Atom")
            .field("index", &self.index())
            .field("atomic_number", &self.atomic_number())
            .finish()
    }
}
