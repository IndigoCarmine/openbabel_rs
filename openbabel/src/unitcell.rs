//! A borrowed view of a molecule's crystallographic unit cell.
//!
//! Present after reading a crystal structure (CIF, etc.); obtain it from
//! [`Molecule::unit_cell`](crate::Molecule::unit_cell).

use openbabel_sys::ffi;

/// The Bravais lattice type of a [`UnitCell`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LatticeType {
    /// Type could not be determined.
    Undefined,
    /// Triclinic (a≠b≠c, α≠β≠γ).
    Triclinic,
    /// Monoclinic.
    Monoclinic,
    /// Orthorhombic.
    Orthorhombic,
    /// Tetragonal.
    Tetragonal,
    /// Rhombohedral (also called trigonal).
    Rhombohedral,
    /// Hexagonal.
    Hexagonal,
    /// Cubic.
    Cubic,
}

impl LatticeType {
    fn from_code(code: u32) -> Self {
        match code {
            1 => LatticeType::Triclinic,
            2 => LatticeType::Monoclinic,
            3 => LatticeType::Orthorhombic,
            4 => LatticeType::Tetragonal,
            5 => LatticeType::Rhombohedral,
            6 => LatticeType::Hexagonal,
            7 => LatticeType::Cubic,
            _ => LatticeType::Undefined,
        }
    }
}

/// A crystallographic unit cell, borrowed from its parent molecule.
///
/// Lengths are in ångström, angles in degrees.
#[derive(Clone, Copy)]
pub struct UnitCell<'mol> {
    mol: &'mol ffi::Molecule,
}

impl<'mol> UnitCell<'mol> {
    pub(crate) fn new(mol: &'mol ffi::Molecule) -> Self {
        UnitCell { mol }
    }

    fn parameters(&self) -> Vec<f64> {
        crate::with_ob(|| ffi::mol_cell_parameters(self.mol))
    }

    /// Cell edge lengths `(a, b, c)` in ångström.
    pub fn lengths(&self) -> (f64, f64, f64) {
        let p = self.parameters();
        (p[0], p[1], p[2])
    }

    /// Cell angles `(alpha, beta, gamma)` in degrees.
    pub fn angles(&self) -> (f64, f64, f64) {
        let p = self.parameters();
        (p[3], p[4], p[5])
    }

    /// Cell volume in ångström³.
    pub fn volume(&self) -> f64 {
        crate::with_ob(|| ffi::mol_cell_volume(self.mol))
    }

    /// Hermann–Mauguin space-group name (e.g. `"P 1"`), or empty if unset.
    pub fn space_group(&self) -> String {
        crate::with_ob(|| ffi::mol_cell_spacegroup(self.mol))
    }

    /// The Bravais [`LatticeType`].
    pub fn lattice_type(&self) -> LatticeType {
        LatticeType::from_code(crate::with_ob(|| ffi::mol_cell_lattice_type(self.mol)))
    }

    /// Convert Cartesian coordinates (Å) to fractional cell coordinates.
    pub fn to_fractional(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let v = crate::with_ob(|| ffi::mol_cell_to_fractional(self.mol, x, y, z));
        (v[0], v[1], v[2])
    }

    /// Convert fractional cell coordinates to Cartesian coordinates (Å).
    pub fn to_cartesian(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let v = crate::with_ob(|| ffi::mol_cell_to_cartesian(self.mol, x, y, z));
        (v[0], v[1], v[2])
    }
}

impl std::fmt::Debug for UnitCell<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (a, b, c) = self.lengths();
        let (alpha, beta, gamma) = self.angles();
        f.debug_struct("UnitCell")
            .field("a", &a)
            .field("b", &b)
            .field("c", &c)
            .field("alpha", &alpha)
            .field("beta", &beta)
            .field("gamma", &gamma)
            .field("lattice_type", &self.lattice_type())
            .finish()
    }
}
