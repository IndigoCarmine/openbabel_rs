//! Periodic-table element data, from OpenBabel's `OBElements` tables.
//!
//! ```
//! assert_eq!(openbabel::elements::symbol(6), "C");
//! assert_eq!(openbabel::elements::atomic_number("Cl"), 17);
//! ```

use openbabel_sys::ffi;

/// Element symbol for an atomic number, e.g. `6` → `"C"` (empty if unknown).
pub fn symbol(atomic_number: u32) -> String {
    crate::with_ob(|| ffi::element_symbol(atomic_number))
}

/// Full element name, e.g. `6` → `"Carbon"`.
pub fn name(atomic_number: u32) -> String {
    crate::with_ob(|| ffi::element_name(atomic_number))
}

/// Atomic number for an element symbol, e.g. `"Cl"` → `17` (`0` if unknown).
pub fn atomic_number(symbol: &str) -> u32 {
    crate::with_ob(|| ffi::element_atomic_number(symbol))
}

/// Standard atomic weight in g/mol.
pub fn mass(atomic_number: u32) -> f64 {
    crate::with_ob(|| ffi::element_mass(atomic_number))
}

/// Mass of the most abundant isotope.
pub fn exact_mass(atomic_number: u32) -> f64 {
    crate::with_ob(|| ffi::element_exact_mass(atomic_number))
}

/// Pauling electronegativity (`0.0` if undefined).
pub fn electronegativity(atomic_number: u32) -> f64 {
    crate::with_ob(|| ffi::element_electronegativity(atomic_number))
}

/// Covalent radius in Ångström.
pub fn covalent_radius(atomic_number: u32) -> f64 {
    crate::with_ob(|| ffi::element_covalent_radius(atomic_number))
}

/// Van der Waals radius in Ångström.
pub fn vdw_radius(atomic_number: u32) -> f64 {
    crate::with_ob(|| ffi::element_vdw_radius(atomic_number))
}

/// Maximum number of bonds typically formed by this element.
pub fn max_bonds(atomic_number: u32) -> u32 {
    crate::with_ob(|| ffi::element_max_bonds(atomic_number))
}
