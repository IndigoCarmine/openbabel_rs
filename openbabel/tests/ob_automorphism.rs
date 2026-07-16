//! Port of OpenBabel's `test/automorphismtest.cpp` — automorphism counting.
//!
//! `FindAutomorphisms(mol).size()` maps directly to `Molecule::automorphisms().len()`.
//! Part 1 is an inline SMILES; parts 2–10 are the Hao & Xu (Computers &
//! Chemistry 26 (2002) 119-123, fig. 2) test structures `hao_xu_1..9.mol`.

mod common;

use common::ob_test_file;
use openbabel::Molecule;

/// Part 1.
#[test]
fn inline_ring_has_eight_automorphisms() {
    let mol = Molecule::parse("C1C(CC2CC2)C1", "smi").expect("parse");
    assert_eq!(mol.automorphisms().len(), 8);
}

/// Parts 2–10: Hao & Xu fig. 2 structures 1–9.
#[test]
fn hao_xu_structures() {
    let cases = [
        ("hao_xu_1.mol", 8),
        ("hao_xu_2.mol", 2),
        ("hao_xu_3.mol", 48),
        ("hao_xu_4.mol", 2),
        ("hao_xu_5.mol", 2),
        ("hao_xu_6.mol", 6),
        ("hao_xu_7.mol", 1),
        ("hao_xu_8.mol", 1),
        ("hao_xu_9.mol", 20),
    ];
    for (file, expected) in cases {
        let path = ob_test_file(file);
        let mol = Molecule::read_file(path.to_str().unwrap(), Some("mol"))
            .unwrap_or_else(|_| panic!("read {file}"));
        assert_eq!(
            mol.automorphisms().len(),
            expected,
            "automorphism count for {file}"
        );
    }
}
