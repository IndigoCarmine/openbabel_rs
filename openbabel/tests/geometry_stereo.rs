//! Integration tests for inter-atom distance and 2D wedge/hash stereo (T20).

use openbabel::Molecule;

#[test]
fn distance_between_bonded_atoms_is_a_bond_length() {
    // Ethanol in 3D; the C0–C1 distance should be a typical single-bond length.
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    mol.add_hydrogens();
    assert!(mol.generate_3d());

    let d = mol.distance(0, 1);
    assert!((1.3..1.7).contains(&d), "C–C distance {d} Å out of range");

    // Distance is symmetric, and self-distance is zero.
    assert!((mol.distance(0, 1) - mol.distance(1, 0)).abs() < 1e-9);
    assert_eq!(mol.distance(0, 0), 0.0);

    // Invalid index yields 0.0.
    assert_eq!(mol.distance(0, 999), 0.0);
}

#[test]
fn wedge_and_hash_flags_round_trip() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");

    // No stereo bond markings by default.
    assert!(!mol.bond(0).unwrap().is_wedge());
    assert!(!mol.bond(0).unwrap().is_hash());

    mol.bond_mut(0).unwrap().set_wedge(true);
    assert!(mol.bond(0).unwrap().is_wedge());
    assert!(!mol.bond(0).unwrap().is_hash());

    // Marking the next bond as a hash leaves bond 0 alone.
    mol.bond_mut(1).unwrap().set_hash(true);
    assert!(mol.bond(1).unwrap().is_hash());
    assert!(mol.bond(0).unwrap().is_wedge());

    // Clearing works.
    mol.bond_mut(0).unwrap().set_wedge(false);
    assert!(!mol.bond(0).unwrap().is_wedge());
}
