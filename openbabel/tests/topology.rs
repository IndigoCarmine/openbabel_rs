//! Integration tests for geometry & topology niche operations (T18).

use openbabel::Molecule;
use std::collections::HashSet;
use std::f64::consts::PI;

#[test]
fn set_torsion_updates_the_dihedral() {
    // Butane backbone C0–C1–C2–C3; set its dihedral to 180° (anti).
    let mut mol = Molecule::parse("CCCC", "smi").expect("butane");
    mol.add_hydrogens();
    assert!(mol.generate_3d());

    mol.set_torsion(0, 1, 2, 3, PI); // radians
    let deg = mol.torsion(0, 1, 2, 3).abs(); // torsion() reports degrees
    assert!((deg - 180.0).abs() < 1.0, "expected ~180°, got {deg}");

    // And a gauche angle (60°).
    mol.set_torsion(0, 1, 2, 3, PI / 3.0);
    let deg = mol.torsion(0, 1, 2, 3).abs();
    assert!((deg - 60.0).abs() < 1.0, "expected ~60°, got {deg}");
}

#[test]
fn find_children_walks_one_side_of_a_bond() {
    // Butane C0–C1–C2–C3 (heavy atoms 0..=3).
    let mol = Molecule::parse("CCCC", "smi").expect("butane");

    // From the C1–C2 bond, C2's side contains C3 (but not C0 or the endpoints).
    let far = mol.find_children(1, 2);
    assert!(far.contains(&3), "C3 should be on C2's side: {far:?}");
    assert!(!far.contains(&0), "C0 is on the other side: {far:?}");
    assert!(!far.contains(&1) && !far.contains(&2), "endpoints excluded: {far:?}");

    // The mirror bond direction gives C0's side.
    let near = mol.find_children(2, 1);
    assert!(near.contains(&0), "C0 should be on C1's side: {near:?}");
    assert!(!near.contains(&3), "{near:?}");
}

#[test]
fn largest_fragment_excludes_counter_ions() {
    // Ethanol plus a dissociated salt: atoms C0, C1, O2, Na3, Cl4.
    let mol = Molecule::parse("CCO.[Na+].[Cl-]", "smi").expect("parse");
    let frag: HashSet<u32> = mol.largest_fragment_atoms().into_iter().collect();
    assert_eq!(frag, HashSet::from([0, 1, 2]), "largest fragment is ethanol");
}

#[test]
fn total_charge_and_spin_round_trip() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");

    mol.set_total_charge(2);
    assert_eq!(mol.total_charge(), 2);

    mol.set_total_spin_multiplicity(3);
    assert_eq!(mol.spin_multiplicity(), 3);
}
