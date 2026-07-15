//! A borrowed view of a residue within a [`Molecule`](crate::Molecule).
//!
//! Residues carry the biopolymer / PDB substructure grouping (chain, residue
//! name and number, per-atom PDB names) that OpenBabel attaches when it reads
//! formats like PDB, PDBQT, or mmCIF. Small molecules parsed from SMILES have
//! no residues, so [`Molecule::num_residues`](crate::Molecule::num_residues)
//! is `0`.

use openbabel_sys::ffi;

use crate::Atom;

/// A single residue, borrowed from its parent molecule.
///
/// A `Residue` is a lightweight handle (a reference to the molecule plus a
/// 0-based residue index); it cannot outlive the molecule it came from.
#[derive(Clone, Copy)]
pub struct Residue<'mol> {
    mol: &'mol ffi::Molecule,
    /// OpenBabel's 0-based residue index.
    idx: u32,
}

impl<'mol> Residue<'mol> {
    pub(crate) fn new(mol: &'mol ffi::Molecule, idx: u32) -> Self {
        Residue { mol, idx }
    }

    /// 0-based index of this residue within the molecule.
    pub fn index(&self) -> u32 {
        self.idx
    }

    /// Residue name, e.g. `"GLY"`, `"ALA"`, `"HOH"`.
    pub fn name(&self) -> String {
        crate::with_ob(|| ffi::residue_name(self.mol, self.idx))
    }

    /// Residue sequence number (PDB `resSeq`).
    pub fn number(&self) -> i32 {
        crate::with_ob(|| ffi::residue_number(self.mol, self.idx))
    }

    /// Residue sequence number as text; may carry a trailing insertion code.
    pub fn number_string(&self) -> String {
        crate::with_ob(|| ffi::residue_number_string(self.mol, self.idx))
    }

    /// Chain identifier (e.g. `"A"`), or an empty string if unset.
    pub fn chain(&self) -> String {
        crate::with_ob(|| ffi::residue_chain(self.mol, self.idx))
    }

    /// Insertion code, or an empty string if unset.
    pub fn insertion_code(&self) -> String {
        crate::with_ob(|| ffi::residue_insertion_code(self.mol, self.idx))
    }

    /// Number of atoms in this residue.
    pub fn num_atoms(&self) -> u32 {
        crate::with_ob(|| ffi::residue_num_atoms(self.mol, self.idx))
    }

    /// Number of heavy (non-hydrogen) atoms in this residue.
    pub fn num_heavy_atoms(&self) -> u32 {
        crate::with_ob(|| ffi::residue_num_heavy_atoms(self.mol, self.idx))
    }

    /// The atoms belonging to this residue.
    pub fn atoms(&self) -> Vec<Atom<'mol>> {
        let indices = crate::with_ob(|| ffi::residue_atom_indices(self.mol, self.idx));
        indices
            .into_iter()
            .map(|i| Atom {
                mol: self.mol,
                ob_idx: i + 1, // shim returns 0-based; Atom uses 1-based ob_idx
            })
            .collect()
    }
}

impl std::fmt::Debug for Residue<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Residue")
            .field("index", &self.index())
            .field("name", &self.name())
            .field("number", &self.number())
            .field("chain", &self.chain())
            .finish()
    }
}
