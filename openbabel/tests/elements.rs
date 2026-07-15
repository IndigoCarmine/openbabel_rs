//! Integration tests for the element-data tables (T8).

use openbabel::elements;

#[test]
fn symbol_and_number_roundtrip() {
    assert_eq!(elements::symbol(6), "C");
    assert_eq!(elements::symbol(17), "Cl");
    assert_eq!(elements::atomic_number("C"), 6);
    assert_eq!(elements::atomic_number("Cl"), 17);
    assert_eq!(elements::atomic_number("Fe"), 26);
}

#[test]
fn names_and_masses() {
    assert_eq!(elements::name(8), "Oxygen");
    assert!((elements::mass(6) - 12.011).abs() < 0.01, "C weight");
    assert!((elements::exact_mass(6) - 12.0).abs() < 1e-6, "12C is exactly 12");
    assert!((elements::mass(1) - 1.008).abs() < 0.01, "H weight");
}

#[test]
fn periodic_properties() {
    // Fluorine is more electronegative than carbon.
    assert!(elements::electronegativity(9) > elements::electronegativity(6));
    // Radii are positive and vdW exceeds covalent.
    assert!(elements::covalent_radius(6) > 0.0);
    assert!(elements::vdw_radius(6) > elements::covalent_radius(6));
    // Carbon typically forms up to 4 bonds.
    assert_eq!(elements::max_bonds(6), 4);
}
