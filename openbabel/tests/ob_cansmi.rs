//! Port of OpenBabel's `test/cansmi.cpp` — canonical SMILES generation.
//!
//! For every molecule in `nci.smi`: read it, write it as canonical SMILES, read
//! that back, and assert the heavy-atom counts match (a write→read→write
//! invariant). `nci.smi` holds ~1000 molecules, so this is one of the heavier
//! ports.

mod common;

use common::read_test_file;
use openbabel::Molecule;

#[test]
fn canonical_smiles_roundtrip_preserves_heavy_atoms() {
    let content = read_test_file("nci.smi");

    let mut molecules = 0usize;
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        // A failed read is a `not ok` in the C++ (i.e. a test failure).
        let mol = Molecule::parse(line, "smi")
            .unwrap_or_else(|_| panic!("SMILES read failed for {line:?}"));
        if mol.num_atoms() == 0 {
            continue;
        }
        molecules += 1;

        let canonical = mol.write("can").expect("write canonical SMILES");
        let reparsed = Molecule::parse(canonical.trim(), "smi")
            .unwrap_or_else(|_| panic!("re-reading canonical SMILES failed: {canonical:?}"));

        assert_eq!(
            mol.num_heavy_atoms(),
            reparsed.num_heavy_atoms(),
            "heavy-atom count changed on canonical roundtrip of {line:?} (canonical {canonical:?})"
        );
    }

    assert!(molecules > 0, "no molecules were read from nci.smi");
}
