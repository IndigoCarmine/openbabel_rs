//! Integration tests for atom relations, persistent ids, and LSSR (T21).

use openbabel::Molecule;

#[test]
fn connectivity_relations() {
    // Butane C0–C1–C2–C3 (heavy atoms).
    let mol = Molecule::parse("CCCC", "smi").expect("butane");
    let a = |i| mol.atom(i).unwrap();

    // 1-2 (directly bonded).
    assert!(a(0).is_connected(&a(1)));
    assert!(!a(0).is_connected(&a(2)));

    // 1-3: share a common neighbour (C0 and C2 share C1).
    assert!(a(0).is_one_three(&a(2)));
    assert!(!a(0).is_one_three(&a(1)));
    assert!(!a(0).is_one_three(&a(3))); // no shared neighbour

    // 1-4: a neighbour of one is bonded to a neighbour of the other
    // (C0's C1 is bonded to C3's C2).
    assert!(a(0).is_one_four(&a(3)));

    // Atoms too far apart to be 1-4 (pentane C0…C4: C1 not bonded to C3).
    let pentane = Molecule::parse("CCCCC", "smi").expect("pentane");
    assert!(!pentane.atom(0).unwrap().is_one_four(&pentane.atom(4).unwrap()));
}

#[test]
fn persistent_id_round_trips() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    mol.atom_mut(0).unwrap().set_id(42);
    assert_eq!(mol.atom(0).unwrap().id(), 42);
}

#[test]
fn ids_survive_atom_deletion() {
    // Propane C0–C1–C2; tag each atom, then delete the first.
    let mut mol = Molecule::parse("CCC", "smi").expect("propane");
    mol.atom_mut(0).unwrap().set_id(100);
    mol.atom_mut(1).unwrap().set_id(101);
    mol.atom_mut(2).unwrap().set_id(102);

    mol.delete_atom(0); // indices shift down, ids do not

    assert_eq!(mol.num_atoms(), 2);
    // The atoms formerly at index 1 and 2 keep their ids at their new indices.
    assert_eq!(mol.atom(0).unwrap().id(), 101);
    assert_eq!(mol.atom(1).unwrap().id(), 102);
}

#[test]
fn lssr_ring_sizes_are_reported() {
    // Benzene: a single six-membered ring.
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    assert_eq!(benzene.lssr_ring_sizes(), vec![6]);

    // Naphthalene: two fused six-membered rings.
    let naph = Molecule::parse("c1ccc2ccccc2c1", "smi").expect("naphthalene");
    let mut sizes = naph.lssr_ring_sizes();
    sizes.sort_unstable();
    assert_eq!(sizes, vec![6, 6]);

    // A chain has no rings.
    let hexane = Molecule::parse("CCCCCC", "smi").expect("hexane");
    assert!(hexane.lssr_ring_sizes().is_empty());
}
