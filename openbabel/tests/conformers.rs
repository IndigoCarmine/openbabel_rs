//! Integration tests for conformer search (T7, OBConformerSearch).
//!
//! Conformer search is a rotor-based genetic algorithm, so it requires a 3D
//! structure and is only meaningful for flexible molecules. It relies on the
//! Eigen support enabled in T5.

use openbabel::Molecule;

#[test]
fn conformer_search_generates_conformers() {
    let mut mol = Molecule::parse("CCCCCC", "smi").expect("hexane"); // flexible
    mol.add_hydrogens();
    assert!(mol.generate_3d(), "a 3D structure is required first");

    let n = mol.generate_conformers(5);
    assert!(n >= 1, "expected at least one conformer, got {n}");
    assert_eq!(mol.num_conformers(), n);

    // Selecting each conformer must stay in range and not panic.
    for i in 0..n {
        mol.set_conformer(i);
    }
}

#[test]
fn each_conformer_has_a_finite_energy() {
    let mut mol = Molecule::parse("CCCCCC", "smi").expect("hexane");
    mol.add_hydrogens();
    assert!(mol.generate_3d());

    let n = mol.generate_conformers(5);
    for i in 0..n {
        mol.set_conformer(i);
        if let Some(e) = mol.energy("MMFF94") {
            assert!(e.is_finite(), "conformer {i} has non-finite energy");
        }
    }
}

#[test]
fn rigid_molecule_yields_a_single_conformer() {
    // Benzene has no rotatable bonds, so a search leaves one conformer.
    let mut mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    mol.add_hydrogens();
    assert!(mol.generate_3d());

    let n = mol.generate_conformers(10);
    assert_eq!(n, 1, "a rigid molecule has a single conformer, got {n}");
}
