//! Integration tests for the T3 "3D & force field" surface: coordinate
//! generation, force-field single-point energy, and geometry optimization.

use openbabel::Molecule;

#[test]
fn generate_3d_gives_coordinates() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse ethanol");
    // From SMILES the molecule has no coordinates.
    assert_eq!(mol.dimension(), 0);
    assert!(!mol.has_3d());

    assert!(mol.generate_3d(), "3D generation should succeed");
    assert_eq!(mol.dimension(), 3);
    assert!(mol.has_3d());

    // gen3d makes hydrogens explicit; ethanol becomes C2H6O = 9 atoms.
    assert_eq!(mol.num_atoms(), 9);

    // At least one atom should now have non-zero coordinates.
    let any_placed = mol
        .atoms()
        .any(|a| a.coords() != (0.0, 0.0, 0.0));
    assert!(any_placed, "expected real 3D coordinates");
}

#[test]
fn forcefield_energy_and_optimization() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse ethanol");
    assert!(mol.generate_3d());

    // MMFF94 reports energies in kcal/mol.
    assert_eq!(
        openbabel::forcefield_energy_unit("MMFF94").as_deref(),
        Some("kcal/mol")
    );

    let e0 = mol.energy("MMFF94").expect("MMFF94 energy");
    assert!(e0.is_finite(), "energy should be finite, got {e0}");

    // Optimizing should not raise the energy (it converges downhill).
    let e1 = mol
        .optimize_geometry("MMFF94", 500)
        .expect("MMFF94 optimization");
    assert!(e1.is_finite());
    assert!(
        e1 <= e0 + 1e-6,
        "optimized energy {e1} should not exceed initial {e0}"
    );
}

#[test]
fn unknown_forcefield_is_none() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    mol.generate_3d();
    assert!(mol.energy("NoSuchField").is_none());
    assert!(openbabel::forcefield_energy_unit("NoSuchField").is_none());
}

#[test]
fn uff_also_works() {
    let mut mol = Molecule::parse("c1ccccc1", "smi").expect("parse benzene");
    assert!(mol.generate_3d());
    let e = mol.energy("UFF").expect("UFF energy");
    assert!(e.is_finite());
}
