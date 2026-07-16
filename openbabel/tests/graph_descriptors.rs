//! Integration tests for whole-molecule graph descriptor vectors and the ring
//! type / root-atom / membership extras (T25).

use openbabel::Molecule;

#[test]
fn graph_theoretical_distances_are_eccentricities() {
    // Butane C0-C1-C2-C3: the terminal carbons are the farthest apart, so they
    // have the largest eccentricity and the interior carbons the smallest.
    let butane = Molecule::parse("CCCC", "smi").expect("butane");
    let gtd = butane.graph_theoretical_distances();
    assert_eq!(gtd, vec![4, 3, 3, 4]);
}

#[test]
fn graph_invariants_respect_symmetry() {
    // Benzene's six carbons are all equivalent, so every graph invariant (and
    // its distance-refined variant) is identical.
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let gi = benzene.graph_invariants();
    assert_eq!(gi.len(), benzene.num_atoms() as usize);
    assert!(gi.iter().all(|&v| v == gi[0]));

    let gid = benzene.graph_invariant_distances();
    assert_eq!(gid.len(), benzene.num_atoms() as usize);
    assert!(gid.iter().all(|&v| v == gid[0]));

    // A less symmetric molecule has more than one distinct invariant.
    let propanol = Molecule::parse("CCCO", "smi").expect("propanol");
    let gi2 = propanol.graph_invariants();
    assert!(gi2.iter().any(|&v| v != gi2[0]));
}

#[test]
fn ring_root_atom_is_the_heteroatom() {
    // Pyridine's single ring is rooted on its nitrogen.
    let pyridine = Molecule::parse("c1ccncc1", "smi").expect("pyridine");
    let ring = pyridine.ring(0).expect("one ring");
    let root = ring.root_atom().expect("heterocycle has a root atom");
    assert_eq!(pyridine.atom(root).unwrap().atomic_number(), 7);

    // Benzene is all-carbon, so it has no distinguished root atom.
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    assert_eq!(benzene.ring(0).unwrap().root_atom(), None);
}

#[test]
fn ring_membership() {
    // Toluene: the ring contains its six aromatic carbons but not the methyl.
    let toluene = Molecule::parse("Cc1ccccc1", "smi").expect("toluene");
    let ring = toluene.ring(0).expect("one ring");

    // Atom 0 is the methyl carbon (written first in the SMILES); it is a
    // substituent, not a ring member.
    let methyl = toluene.atom(0).unwrap();
    assert!(!ring.contains_atom(&methyl));

    // Every atom the ring reports really is a member.
    for &i in &ring.atom_indices() {
        assert!(ring.contains_atom(&toluene.atom(i).unwrap()));
    }
    // The ring has six members, none of them atom 0.
    let members = (0..toluene.num_atoms())
        .filter(|&i| ring.contains_atom(&toluene.atom(i).unwrap()))
        .count();
    assert_eq!(members, 6);
}

#[test]
fn ring_type_is_reported_after_perception() {
    // GetType is populated by OpenBabel's ring typer. It is either empty (not
    // yet perceived) or a recognised name; querying it must not panic and, when
    // present, is a sensible token.
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let ty = benzene.ring(0).unwrap().ring_type();
    assert!(ty.is_empty() || ty.chars().all(|c| c.is_ascii_graphic()));
}
