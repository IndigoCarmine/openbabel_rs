//! Port of OpenBabel's `test/isomorphismtest.cpp` — subgraph isomorphism.
//!
//! `MapUnique` maps to [`Molecule::substructure_search`] (unique mappings). Only
//! the parts that use `MapUnique` (or `FindAutomorphisms` without a mask) are
//! ported; the C++'s `MapAll`/`MapFirst` counts and the `OBBitVec`-masked
//! subtests (parts 1, 5, 6, 7, 9) need APIs the binding does not expose.

mod common;

use openbabel::Molecule;

/// testIsomorphism2 (part 2): benzene occurs once (uniquely) in p-xylene.
#[test]
fn benzene_in_p_xylene() {
    let mol = Molecule::parse("Cc1ccc(C)cc1", "smi").expect("mol");
    let query = Molecule::parse("c1ccccc1", "smi").expect("query");
    assert_eq!(mol.substructure_search(&query).len(), 1);
}

/// testIsomorphism3 (part 3): two cyclopropane rings.
#[test]
fn cyclopropane_twice() {
    let mol = Molecule::parse("C1CC1CC1CC1", "smi").expect("mol");
    let query = Molecule::parse("C1CC1", "smi").expect("query");
    assert_eq!(mol.substructure_search(&query).len(), 2);
}

/// testIsomorphism4 (part 4): bicyclobutane contains cyclopropane twice.
#[test]
fn cyclopropane_in_bicyclobutane() {
    let mol = Molecule::parse("C12C(C2)C1", "smi").expect("mol");
    let query = Molecule::parse("C1CC1", "smi").expect("query");
    assert_eq!(mol.substructure_search(&query).len(), 2);
}

/// testAutomorphismPreMapping (part 8).
#[test]
fn pentamethylbenzene_automorphisms() {
    let mol = Molecule::parse("c1(C)c(C)c(C)c(C)c(C)c1", "smi").expect("mol");
    assert_eq!(mol.automorphisms().len(), 2);
}
