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
        crate::with_ob(|| ffi::bond_begin_idx(self.mol, self.ob_idx)).saturating_sub(1)
    }

    /// 0-based index of the atom at the end of this bond.
    pub fn end_atom_index(&self) -> u32 {
        crate::with_ob(|| ffi::bond_end_idx(self.mol, self.ob_idx)).saturating_sub(1)
    }

    /// Bond order (1 = single, 2 = double, 3 = triple; 5 denotes aromatic in
    /// OpenBabel's convention).
    pub fn order(&self) -> u32 {
        crate::with_ob(|| ffi::bond_order(self.mol, self.ob_idx))
    }

    /// Whether this bond is part of an aromatic system.
    pub fn is_aromatic(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_aromatic(self.mol, self.ob_idx))
    }

    /// Whether this bond is a member of any ring.
    pub fn is_in_ring(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_in_ring(self.mol, self.ob_idx))
    }

    /// Whether this bond carries cis/trans (double-bond) stereochemistry.
    pub fn is_cistrans_stereo(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_cistrans_stereo(self.mol, self.ob_idx))
    }

    /// Geometric length of the bond (needs coordinates; `0` without them).
    pub fn length(&self) -> f64 {
        crate::with_ob(|| ffi::bond_length(self.mol, self.ob_idx))
    }

    /// Equilibrium bond length from OpenBabel's tabulated covalent radii.
    pub fn equilibrium_length(&self) -> f64 {
        crate::with_ob(|| ffi::bond_equilibrium_length(self.mol, self.ob_idx))
    }

    /// Whether this bond is freely rotatable (single, acyclic, non-terminal).
    pub fn is_rotor(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_rotor(self.mol, self.ob_idx))
    }

    /// Whether this bond is the C–N bond of an amide.
    pub fn is_amide(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_amide(self.mol, self.ob_idx))
    }

    /// Whether this bond is the C–O single bond of an ester.
    pub fn is_ester(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_ester(self.mol, self.ob_idx))
    }

    /// Whether this bond is a carbonyl (C=O) bond.
    pub fn is_carbonyl(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_carbonyl(self.mol, self.ob_idx))
    }

    /// Whether this bond is a ring-closure bond (as written in SMILES).
    pub fn is_closure(&self) -> bool {
        crate::with_ob(|| ffi::bond_is_closure(self.mol, self.ob_idx))
    }

    /// The atom at the other end of this bond from `atom`, or `None` if `atom`
    /// is not one of this bond's endpoints.
    pub fn other_atom(&self, atom: &crate::Atom<'mol>) -> Option<crate::Atom<'mol>> {
        let idx = crate::with_ob(|| ffi::bond_other_atom(self.mol, self.ob_idx, atom.index()));
        if idx < 0 {
            None
        } else {
            Some(crate::atom::Atom {
                mol: self.mol,
                ob_idx: idx as u32 + 1,
            })
        }
    }

    /// A string annotation previously attached under `key` (see
    /// [`BondMut::set_data`]), or `None` if this bond has none.
    pub fn data(&self, key: &str) -> Option<String> {
        let mut ok = false;
        let value = crate::with_ob(|| ffi::bond_get_data(self.mol, self.ob_idx, key, &mut ok));
        ok.then_some(value)
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

/// A mutable handle to a bond, for editing it.
///
/// Obtained from [`Molecule::bond_mut`](crate::Molecule::bond_mut); it borrows
/// the molecule mutably, so only one exists at a time.
pub struct BondMut<'m> {
    mol: &'m mut crate::Molecule,
    /// OpenBabel's 0-based bond index.
    ob_idx: u32,
}

impl<'m> BondMut<'m> {
    pub(crate) fn new(mol: &'m mut crate::Molecule, ob_idx: u32) -> Self {
        BondMut { mol, ob_idx }
    }

    /// 0-based index of this bond within the molecule.
    pub fn index(&self) -> u32 {
        self.ob_idx
    }

    /// Set the bond order (1 = single, 2 = double, 3 = triple; 5 = aromatic in
    /// OpenBabel's convention).
    pub fn set_order(&mut self, order: u32) -> &mut Self {
        crate::with_ob(|| ffi::bond_set_order(self.mol.as_inner_pin_mut(), self.ob_idx, order));
        self
    }

    /// Move the end atom so the bond has geometric length `length` (Å), keeping
    /// the begin atom fixed. Returns `false` if the bond is invalid.
    pub fn set_length(&mut self, length: f64) -> bool {
        crate::with_ob(|| ffi::bond_set_length(self.mol.as_inner_pin_mut(), self.ob_idx, length))
    }

    /// Attach (or replace) an arbitrary string annotation under `key`, readable
    /// later with [`Bond::data`].
    pub fn set_data(&mut self, key: &str, value: &str) -> &mut Self {
        crate::with_ob(|| ffi::bond_set_data(self.mol.as_inner_pin_mut(), self.ob_idx, key, value));
        self
    }
}
