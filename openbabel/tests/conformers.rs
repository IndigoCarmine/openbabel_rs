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
fn conformer_coordinates_read_without_switching_active() {
    let mut mol = Molecule::parse("CCCCCC", "smi").expect("hexane");
    mol.add_hydrogens();
    assert!(mol.generate_3d());
    let n = mol.generate_conformers(5);

    mol.set_conformer(0);
    let active_before: Vec<(f64, f64, f64)> = mol.atoms().map(|a| a.coords()).collect();

    for i in 0..n {
        let coords = mol.conformer_coordinates(i).expect("in range");
        assert_eq!(coords.len() as u32, mol.num_atoms());
        assert!(coords.iter().flatten().all(|v| v.is_finite()));
    }
    assert!(mol.conformer_coordinates(n).is_none(), "out-of-range must be None");

    // Reading conformers must not disturb the active one.
    let active_after: Vec<(f64, f64, f64)> = mol.atoms().map(|a| a.coords()).collect();
    assert_eq!(active_before, active_after);
}

#[test]
fn conformer_energies_are_scored_per_conformer() {
    let mut mol = Molecule::parse("CCCCCC", "smi").expect("hexane");
    mol.add_hydrogens();
    assert!(mol.generate_3d());
    let n = mol.generate_conformers(5);
    mol.set_conformer(0);

    // One finite energy per conformer.
    let batch = mol.conformer_energies("MMFF94");
    assert_eq!(batch.len() as u32, n);
    assert!(batch.iter().all(|e| e.is_finite()));

    // Deterministic: the same molecule scores identically every time.
    assert_eq!(batch, mol.conformer_energies("MMFF94"));

    // Distinct conformers really are scored separately (not one value repeated),
    // which is only true if each conformer's coordinates are read.
    if n > 1 {
        assert!(
            batch.iter().any(|&e| (e - batch[0]).abs() > 1e-6),
            "all {n} conformer energies were identical: {batch:?}",
        );
    }

    // Scoring is non-mutating: the active conformer's coordinates are unchanged.
    let before: Vec<(f64, f64, f64)> = mol.atoms().map(|a| a.coords()).collect();
    let _ = mol.conformer_energies("MMFF94");
    let after: Vec<(f64, f64, f64)> = mol.atoms().map(|a| a.coords()).collect();
    assert_eq!(before, after);

    // An unknown force field yields no energies.
    assert!(mol.conformer_energies("NOPE").is_empty());
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
