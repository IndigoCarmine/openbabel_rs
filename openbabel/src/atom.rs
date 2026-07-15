//! A borrowed view of an atom within a [`Molecule`](crate::Molecule).

use openbabel_sys::ffi;

/// The winding of a tetrahedral stereocenter: the sense (clockwise or
/// anticlockwise) in which its neighbours are arranged when viewed from a
/// reference direction. This is OpenBabel's internal descriptor, not a CIP
/// `R`/`S` label.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Winding {
    /// The neighbours run clockwise when viewed from the reference direction.
    Clockwise,
    /// The neighbours run anticlockwise when viewed from the reference direction.
    AntiClockwise,
}

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

    /// Whether this atom is a tetrahedral stereocenter.
    pub fn is_tetrahedral_stereo(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_tetrahedral_stereo(self.mol, self.ob_idx))
    }

    /// The [`Winding`] of this atom's tetrahedral stereocenter, or `None` if it
    /// is not a specified stereocenter.
    pub fn stereo_winding(&self) -> Option<Winding> {
        match crate::with_ob(|| ffi::atom_tetrahedral_winding(self.mol, self.ob_idx)) {
            1 => Some(Winding::Clockwise),
            2 => Some(Winding::AntiClockwise),
            _ => None,
        }
    }

    /// OpenBabel's internal atom type string (e.g. `"C3"`, `"O2"`, `"Nam"`).
    pub fn type_name(&self) -> String {
        crate::with_ob(|| ffi::atom_type(self.mol, self.ob_idx))
    }

    /// Isotope number, or `0` for the element's default isotopic mix.
    pub fn isotope(&self) -> u32 {
        crate::with_ob(|| ffi::atom_isotope(self.mol, self.ob_idx))
    }

    /// Standard atomic mass of this atom (accounts for its isotope, if set).
    pub fn atomic_mass(&self) -> f64 {
        crate::with_ob(|| ffi::atom_atomic_mass(self.mol, self.ob_idx))
    }

    /// Exact (isotopic) mass of this atom.
    pub fn exact_mass(&self) -> f64 {
        crate::with_ob(|| ffi::atom_exact_mass(self.mol, self.ob_idx))
    }

    /// Spin multiplicity (0 if unset; 2 = radical, 3 = carbene/triplet, …).
    pub fn spin_multiplicity(&self) -> i32 {
        crate::with_ob(|| ffi::atom_spin_multiplicity(self.mol, self.ob_idx))
    }

    /// Number of heavy-atom (non-hydrogen) neighbours.
    pub fn heavy_degree(&self) -> u32 {
        crate::with_ob(|| ffi::atom_heavy_degree(self.mol, self.ob_idx))
    }

    /// Number of heteroatom (non-C, non-H) neighbours.
    pub fn hetero_degree(&self) -> u32 {
        crate::with_ob(|| ffi::atom_hetero_degree(self.mol, self.ob_idx))
    }

    /// Whether this atom is a chiral center.
    pub fn is_chiral(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_chiral(self.mol, self.ob_idx))
    }

    /// Whether this atom is a heteroatom (not carbon or hydrogen).
    pub fn is_heteroatom(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_heteroatom(self.mol, self.ob_idx))
    }

    /// Whether this atom is a metal.
    pub fn is_metal(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_metal(self.mol, self.ob_idx))
    }

    /// Whether this is a hydrogen bonded to N, O, P, or S (a polar hydrogen).
    pub fn is_polar_hydrogen(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_polar_hydrogen(self.mol, self.ob_idx))
    }

    /// Number of rings this atom belongs to.
    pub fn ring_count(&self) -> u32 {
        crate::with_ob(|| ffi::atom_member_of_ring_count(self.mol, self.ob_idx))
    }

    /// Size of the smallest ring containing this atom, or `0` if it is in none.
    pub fn smallest_ring_size(&self) -> u32 {
        crate::with_ob(|| ffi::atom_member_of_ring_size(self.mol, self.ob_idx))
    }

    /// Whether this atom is a member of a ring of exactly `size` atoms.
    pub fn is_in_ring_size(&self, size: u32) -> bool {
        crate::with_ob(|| ffi::atom_is_in_ring_size(self.mol, self.ob_idx, size))
    }

    /// The [`Residue`](crate::Residue) this atom belongs to, or `None` if the
    /// molecule carries no residue information (typical for small molecules
    /// parsed from SMILES).
    pub fn residue(&self) -> Option<crate::Residue<'mol>> {
        let idx = crate::with_ob(|| ffi::atom_residue_index(self.mol, self.ob_idx));
        if idx < 0 {
            None
        } else {
            Some(crate::residue::Residue::new(self.mol, idx as u32))
        }
    }

    /// This atom's PDB atom name within its residue (e.g. `" CA "`), or an
    /// empty string if it has no residue.
    pub fn residue_atom_id(&self) -> String {
        crate::with_ob(|| ffi::atom_residue_atom_id(self.mol, self.ob_idx))
    }

    /// Whether this atom is a `HETATM` (heteroatom record) in its residue.
    pub fn is_hetatm(&self) -> bool {
        crate::with_ob(|| ffi::atom_is_hetatm(self.mol, self.ob_idx))
    }

    /// This atom's PDB serial number within its residue, or `0` if it has none.
    pub fn serial_number(&self) -> u32 {
        crate::with_ob(|| ffi::atom_serial_number(self.mol, self.ob_idx))
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
