//! Integration tests for perception-state flags and targeted hydrogen editing (T22).

use openbabel::Molecule;

#[test]
fn hydrogens_added_flag_tracks_state() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    assert!(!mol.has_hydrogens_added());
    mol.add_hydrogens();
    assert!(mol.has_hydrogens_added());
}

#[test]
fn nonzero_coords_flag_tracks_geometry() {
    // Straight from SMILES there are no coordinates.
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    assert!(!mol.has_nonzero_coords());

    mol.add_hydrogens();
    assert!(mol.generate_3d());
    assert!(mol.has_nonzero_coords());
}

#[test]
fn perception_flags_are_queryable() {
    // Touching ring/aromatic perception marks the corresponding flags. Reading
    // rings forces SSSR perception.
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let _ = mol.num_rings(); // forces SSSR
    assert!(mol.has_sssr_perceived());
    // Aromaticity is perceived when the aromatic SMILES is read/queried.
    let _ = mol.atom(0).unwrap().is_aromatic();
    assert!(mol.has_aromatic_perceived());
}

#[test]
fn targeted_hydrogen_editing() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse"); // C0 C1 O2
    let before = mol.num_atoms();

    // Add explicit H to just the oxygen (one hydroxyl H).
    assert!(mol.add_hydrogens_to_atom(2));
    assert_eq!(mol.num_atoms(), before + 1);
    assert_eq!(mol.atom(2).unwrap().explicit_hydrogen_count(), 1);

    // Removing them again drops the count back.
    assert!(mol.remove_hydrogens_of_atom(2));
    assert_eq!(mol.num_atoms(), before);

    // Invalid index is rejected.
    assert!(!mol.add_hydrogens_to_atom(999));
}
