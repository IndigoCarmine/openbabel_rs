//! The [`Molecule`] type — a safe, owning wrapper around OpenBabel's `OBMol`.

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::atom::Atom;
use crate::bond::Bond;
use crate::error::Error;
use crate::with_ob;

/// A molecule.
///
/// Construct one by parsing (`Molecule::parse`) or building up from empty
/// (`Molecule::new`). Dropping it frees the underlying OpenBabel object.
pub struct Molecule {
    inner: UniquePtr<ffi::Molecule>,
}

impl Molecule {
    /// Create an empty molecule.
    pub fn new() -> Self {
        Molecule {
            inner: with_ob(ffi::mol_new),
        }
    }

    /// Parse `data` using the OpenBabel format id `format`
    /// (e.g. `"smi"`, `"mol"`, `"sdf"`, `"pdb"`, `"inchi"`).
    ///
    /// Returns [`Error::Parse`] if the data cannot be read in that format.
    pub fn parse(data: &str, format: &str) -> Result<Self, Error> {
        let inner = with_ob(|| ffi::mol_read(format, data));
        if inner.is_null() {
            Err(Error::Parse {
                format: format.to_string(),
            })
        } else {
            Ok(Molecule { inner })
        }
    }

    /// Serialize this molecule to `format`, returning the text.
    ///
    /// Returns [`Error::UnknownFormat`] if OpenBabel doesn't know `format`.
    pub fn write(&self, format: &str) -> Result<String, Error> {
        let mut ok = true;
        let out = with_ob(|| ffi::mol_write(self.as_inner(), format, &mut ok));
        if ok {
            Ok(out)
        } else {
            Err(Error::UnknownFormat {
                format: format.to_string(),
            })
        }
    }

    /// Molecular formula in Hill order, e.g. `"C2H6O"` (counts implicit H).
    pub fn formula(&self) -> String {
        with_ob(|| ffi::mol_formula(self.as_inner()))
    }

    /// Standard molar mass in g/mol (counts implicit H).
    pub fn molar_mass(&self) -> f64 {
        with_ob(|| ffi::mol_mol_wt(self.as_inner()))
    }

    /// Monoisotopic exact mass (counts implicit H).
    pub fn exact_mass(&self) -> f64 {
        with_ob(|| ffi::mol_exact_mass(self.as_inner()))
    }

    /// Net formal charge of the molecule.
    pub fn total_charge(&self) -> i32 {
        with_ob(|| ffi::mol_total_charge(self.as_inner()))
    }

    /// Number of (explicit) atoms.
    pub fn num_atoms(&self) -> u32 {
        with_ob(|| ffi::mol_num_atoms(self.as_inner()))
    }

    /// Number of bonds.
    pub fn num_bonds(&self) -> u32 {
        with_ob(|| ffi::mol_num_bonds(self.as_inner()))
    }

    /// The molecule's title (often a name or identifier; may be empty).
    pub fn title(&self) -> String {
        with_ob(|| ffi::mol_title(self.as_inner()))
    }

    /// Set the molecule's title.
    pub fn set_title(&mut self, title: &str) {
        with_ob(|| ffi::mol_set_title(self.inner.pin_mut(), title));
    }

    /// Make implicit hydrogens explicit (adds H atoms to the graph).
    pub fn add_hydrogens(&mut self) {
        with_ob(|| ffi::mol_add_hydrogens(self.inner.pin_mut()));
    }

    /// Remove explicit hydrogens (they become implicit again).
    pub fn remove_hydrogens(&mut self) {
        with_ob(|| ffi::mol_delete_hydrogens(self.inner.pin_mut()));
    }

    /// Evaluate a numeric descriptor plugin by id (e.g. `"logP"`, `"TPSA"`,
    /// `"MR"`, `"MW"`). Returns `None` if OpenBabel has no such descriptor.
    pub fn descriptor(&self, id: &str) -> Option<f64> {
        let mut ok = true;
        let value = with_ob(|| ffi::descriptor(self.as_inner(), id, &mut ok));
        if ok {
            Some(value)
        } else {
            None
        }
    }

    /// Predicted octanol/water partition coefficient (logP).
    pub fn logp(&self) -> Option<f64> {
        self.descriptor("logP")
    }

    /// Topological polar surface area (TPSA).
    pub fn tpsa(&self) -> Option<f64> {
        self.descriptor("TPSA")
    }

    /// Molar refractivity (MR).
    pub fn molar_refractivity(&self) -> Option<f64> {
        self.descriptor("MR")
    }

    /// The atom at 0-based `index`, or `None` if out of range.
    pub fn atom(&self, index: u32) -> Option<Atom<'_>> {
        if index < self.num_atoms() {
            Some(Atom {
                mol: self.as_inner(),
                ob_idx: index + 1, // OpenBabel atoms are 1-based.
            })
        } else {
            None
        }
    }

    /// Iterate over the atoms in index order.
    pub fn atoms(&self) -> impl Iterator<Item = Atom<'_>> + '_ {
        let mol = self.as_inner();
        (0..self.num_atoms()).map(move |i| Atom { mol, ob_idx: i + 1 })
    }

    /// The bond at 0-based `index`, or `None` if out of range.
    pub fn bond(&self, index: u32) -> Option<Bond<'_>> {
        if index < self.num_bonds() {
            Some(Bond {
                mol: self.as_inner(),
                ob_idx: index, // OpenBabel bonds are 0-based.
            })
        } else {
            None
        }
    }

    /// Iterate over the bonds in index order.
    pub fn bonds(&self) -> impl Iterator<Item = Bond<'_>> + '_ {
        let mol = self.as_inner();
        (0..self.num_bonds()).map(move |i| Bond { mol, ob_idx: i })
    }

    pub(crate) fn as_inner(&self) -> &ffi::Molecule {
        self.inner.as_ref().expect("Molecule is never null")
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Molecule")
            .field("formula", &self.formula())
            .field("num_atoms", &self.num_atoms())
            .field("num_bonds", &self.num_bonds())
            .finish()
    }
}
