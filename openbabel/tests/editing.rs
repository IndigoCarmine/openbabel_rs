//! Integration tests for molecule construction & editing, multi-molecule I/O,
//! and ring access (T11).

use openbabel::Molecule;

#[test]
fn build_water_from_scratch() {
    let mut mol = Molecule::new();
    mol.begin_modify();
    let o = mol.add_atom(8);
    let h1 = mol.add_atom(1);
    let h2 = mol.add_atom(1);
    assert!(mol.add_bond(o, h1, 1));
    assert!(mol.add_bond(o, h2, 1));
    mol.end_modify();

    assert_eq!(mol.num_atoms(), 3);
    assert_eq!(mol.num_bonds(), 2);
    assert_eq!(mol.formula(), "H2O");
    assert!((mol.molar_mass() - 18.02).abs() < 0.05, "mass {}", mol.molar_mass());
}

#[test]
fn edit_atom_properties() {
    // Set formal charge: ammonia N (neutral) -> ammonium-like +1 total charge.
    let mut mol = Molecule::parse("N", "smi").expect("parse");
    assert_eq!(mol.total_charge(), 0);
    mol.atom_mut(0).expect("atom 0").set_formal_charge(1);
    assert_eq!(mol.total_charge(), 1);

    // Set coordinates and read them back.
    mol.atom_mut(0).unwrap().set_position(1.5, -2.0, 3.25);
    let (x, y, z) = mol.atom(0).unwrap().coords();
    assert!((x - 1.5).abs() < 1e-9 && (y + 2.0).abs() < 1e-9 && (z - 3.25).abs() < 1e-9);

    // Change the element.
    mol.atom_mut(0).unwrap().set_atomic_number(6);
    assert_eq!(mol.atom(0).unwrap().atomic_number(), 6);
}

#[test]
fn delete_atom_and_bond() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    assert_eq!(mol.num_atoms(), 3);
    assert!(mol.delete_atom(2)); // remove the O
    assert_eq!(mol.num_atoms(), 2);
    assert!(!mol.delete_atom(9)); // out of range

    let mut eth = Molecule::parse("CC", "smi").expect("parse");
    assert_eq!(eth.num_bonds(), 1);
    assert!(eth.delete_bond(0));
    assert_eq!(eth.num_bonds(), 0);
}

#[test]
fn set_bond_order() {
    let mut mol = Molecule::parse("CC", "smi").expect("parse"); // ethane, C–C single
    assert_eq!(mol.bond(0).unwrap().order(), 1);
    mol.bond_mut(0).unwrap().set_order(2);
    assert_eq!(mol.bond(0).unwrap().order(), 2);
}

#[test]
fn connect_from_coordinates() {
    // Build a bare C and O at a bonding distance, then perceive connectivity.
    let mut mol = Molecule::new();
    mol.begin_modify();
    let c = mol.add_atom(6);
    let o = mol.add_atom(8);
    mol.end_modify();
    mol.atom_mut(c).unwrap().set_position(0.0, 0.0, 0.0);
    mol.atom_mut(o).unwrap().set_position(1.2, 0.0, 0.0);
    mol.set_dimension(3);
    assert_eq!(mol.num_bonds(), 0);

    mol.connect_the_dots();
    assert_eq!(mol.num_bonds(), 1, "covalent distance should bond C and O");
    mol.perceive_bond_orders();
    assert!(mol.bond(0).unwrap().order() >= 1);
}

#[test]
fn parse_and_write_many() {
    let input = "CCO\nc1ccccc1\nCC(=O)O\n";
    let mols = Molecule::parse_many(input, "smi");
    assert_eq!(mols.len(), 3);
    assert_eq!(mols[0].formula(), "C2H6O"); // ethanol
    assert_eq!(mols[1].formula(), "C6H6"); // benzene

    let out = openbabel::write_many(&mols, "smi").expect("write_many");
    let round = Molecule::parse_many(&out, "smi");
    assert_eq!(round.len(), 3);

    // Unknown format -> empty.
    assert!(Molecule::parse_many(input, "nosuchfmt").is_empty());
}

#[test]
fn ring_access() {
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("parse");
    assert_eq!(benzene.num_rings(), 1);
    let rings: Vec<_> = benzene.rings().collect();
    assert_eq!(rings.len(), 1);
    assert_eq!(rings[0].size(), 6);
    assert_eq!(rings[0].atom_indices().len(), 6);
    assert!(rings[0].is_aromatic());

    // Naphthalene: two fused rings.
    let naph = Molecule::parse("c1ccc2ccccc2c1", "smi").expect("parse");
    assert_eq!(naph.num_rings(), 2);
    assert!(naph.rings().all(|r| r.size() == 6 && r.is_aromatic()));

    // Acyclic molecule: no rings.
    let hexane = Molecule::parse("CCCCCC", "smi").expect("parse");
    assert_eq!(hexane.num_rings(), 0);
    assert!(hexane.ring(0).is_none());
}

#[test]
fn ph_correct_hydrogens() {
    // Acetic acid at physiological pH loses its carboxyl proton -> net -1.
    let mut mol = Molecule::parse("CC(=O)O", "smi").expect("parse");
    assert_eq!(mol.total_charge(), 0);
    assert!(mol.add_hydrogens_for_ph(7.4));
    assert_eq!(mol.total_charge(), -1, "carboxylate at pH 7.4");
}
