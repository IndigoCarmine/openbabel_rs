//! Integration tests for the extended atom/bond/molecule API (T8).

use openbabel::Molecule;

#[test]
fn ring_and_count_queries() {
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    assert_eq!(mol.num_heavy_atoms(), 6);
    assert_eq!(mol.num_rings(), 1);
    assert_eq!(mol.num_rotatable_bonds(), 0);

    let c = mol.atom(0).expect("atom 0");
    assert_eq!(c.ring_count(), 1);
    assert_eq!(c.smallest_ring_size(), 6);
    assert!(c.is_in_ring_size(6));
    assert!(!c.is_in_ring_size(5));
    assert_eq!(c.heavy_degree(), 2);
    assert!(!c.is_heteroatom());
    assert!(!c.is_metal());
}

#[test]
fn atom_type_mass_and_heteroatom() {
    let mol = Molecule::parse("CCO", "smi").expect("ethanol");
    let o = mol.atom(2).expect("atom 2 (O)");
    assert!(o.is_heteroatom());
    assert!(!o.type_name().is_empty(), "OB assigns an atom type");
    assert!((o.atomic_mass() - 15.999).abs() < 0.1);
    assert_eq!(o.isotope(), 0, "default isotopic mix");
}

#[test]
fn rotatable_bond_detection() {
    // Butane C-C-C-C: only the central C-C bond is a rotor.
    let mol = Molecule::parse("CCCC", "smi").expect("butane");
    assert_eq!(mol.num_rotatable_bonds(), 1);
    assert_eq!(mol.bonds().filter(|b| b.is_rotor()).count(), 1);
}

#[test]
fn functional_group_bonds() {
    // Acetamide CC(=O)N has both an amide bond and a carbonyl.
    let mol = Molecule::parse("CC(=O)N", "smi").expect("acetamide");
    assert!(mol.bonds().any(|b| b.is_amide()), "amide bond");
    assert!(mol.bonds().any(|b| b.is_carbonyl()), "carbonyl bond");
}

#[test]
fn spaced_formula_and_heavy_atoms() {
    let mut mol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert_eq!(mol.num_heavy_atoms(), 3);
    let spaced = mol.spaced_formula();
    assert!(spaced.contains("C") && spaced.contains("O") && spaced.contains(' '));
    mol.add_hydrogens();
    assert_eq!(mol.num_heavy_atoms(), 3, "heavy count ignores added H");
}

#[test]
fn bond_length_after_3d() {
    let mut mol = Molecule::parse("CC", "smi").expect("ethane");
    mol.add_hydrogens();
    assert!(mol.generate_3d());
    // The first bond should have a sensible positive length.
    let len = mol.bond(0).expect("bond 0").length();
    assert!(len > 0.5 && len < 2.0, "bond length {len} out of range");
}

#[test]
fn angle_and_torsion_geometry() {
    let mut mol = Molecule::parse("CCCC", "smi").expect("butane");
    mol.add_hydrogens();
    assert!(mol.generate_3d());
    let angle = mol.angle(0, 1, 2);
    assert!(angle > 100.0 && angle < 125.0, "C-C-C angle {angle}");
    let torsion = mol.torsion(0, 1, 2, 3);
    assert!((-180.0..=180.0).contains(&torsion), "torsion {torsion}");
}

#[test]
fn clone_is_independent() {
    let mut a = Molecule::parse("CCO", "smi").expect("ethanol");
    let b = a.clone();
    a.add_hydrogens();
    assert_eq!(a.num_atoms(), 9);
    assert_eq!(b.num_atoms(), 3, "the clone must be unaffected");
}

#[test]
fn separate_and_strip_salts() {
    let mol = Molecule::parse("CC(=O)[O-].[Na+]", "smi").expect("sodium acetate");
    let frags = mol.separate();
    assert_eq!(frags.len(), 2, "two disconnected fragments");
    assert!(frags.iter().any(|f| f.num_atoms() == 1), "one is a lone Na+");

    let mut salt = mol.clone();
    assert!(salt.strip_salts(0), "removing the counterion changes the molecule");
    assert!(
        salt.atoms().all(|a| a.atomic_number() != 11),
        "sodium should be gone"
    );
}

#[test]
fn string_properties_roundtrip() {
    let mut mol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert_eq!(mol.property("source"), None);
    mol.set_property("source", "test-suite");
    assert_eq!(mol.property("source").as_deref(), Some("test-suite"));
}
