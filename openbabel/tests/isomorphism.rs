//! Integration tests for subgraph isomorphism and automorphisms (T15).

use openbabel::Molecule;

#[test]
fn identical_molecule_matches_itself_once() {
    // Ethanol has no symmetry, so it maps onto itself exactly one way.
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    let query = Molecule::parse("CCO", "smi").expect("parse");
    let maps = mol.substructure_search(&query);
    assert_eq!(maps.len(), 1);
    assert_eq!(maps[0].len(), 3);
    // The one mapping is the identity.
    assert_eq!(maps[0], vec![0, 1, 2]);
}

#[test]
fn substructure_maps_to_correct_elements() {
    // The C–O fragment occurs in ethanol (its C1–O2 bond).
    let mol = Molecule::parse("CCO", "smi").expect("parse"); // C0–C1–O2
    let query = Molecule::parse("CO", "smi").expect("parse"); // C–O
    let maps = mol.substructure_search(&query);
    assert!(!maps.is_empty());
    for m in &maps {
        assert_eq!(m.len(), 2);
        // query atom 0 is carbon, atom 1 is oxygen — the mapped target atoms
        // must have those elements.
        assert_eq!(mol.atom(m[0]).unwrap().atomic_number(), 6);
        assert_eq!(mol.atom(m[1]).unwrap().atomic_number(), 8);
    }
}

#[test]
fn has_substructure_positive_and_negative() {
    let ethanol = Molecule::parse("CCO", "smi").expect("parse");
    let methyl = Molecule::parse("C", "smi").expect("parse");
    let butane = Molecule::parse("CCCC", "smi").expect("parse");

    assert!(ethanol.has_substructure(&methyl)); // contains a carbon
    assert!(!ethanol.has_substructure(&butane)); // no 4-carbon chain
}

#[test]
fn benzene_ring_occurs_in_toluene() {
    let toluene = Molecule::parse("Cc1ccccc1", "smi").expect("parse");
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("parse");
    let maps = toluene.substructure_search(&benzene);
    assert!(!maps.is_empty(), "benzene should be found in toluene");
    assert!(maps.iter().all(|m| m.len() == 6));
}

#[test]
fn benzene_has_twelve_automorphisms() {
    // The six-membered carbon ring (implicit H) has the dihedral symmetry of a
    // hexagon: |Aut| = 12.
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("parse");
    let auts = benzene.automorphisms();
    assert_eq!(auts.len(), 12, "got {} automorphisms", auts.len());
    assert!(auts.iter().all(|a| a.len() == 6));
}

#[test]
fn asymmetric_molecule_has_only_the_identity_automorphism() {
    let ethanol = Molecule::parse("CCO", "smi").expect("parse");
    let auts = ethanol.automorphisms();
    assert_eq!(auts.len(), 1);
    assert_eq!(auts[0], vec![0, 1, 2]);
}
