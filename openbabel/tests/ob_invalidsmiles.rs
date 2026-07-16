//! Port of OpenBabel's `test/invalidsmiles.cpp` — invalid SMILES must be
//! rejected.
//!
//! Every line of `invalid-smiles.txt` (including its blank line) must fail to
//! produce a molecule. The C++'s commented-out "random garbage" section
//! (PR#1730132) is not ported, matching the original.

mod common;

use common::read_test_file;
use openbabel::Molecule;

#[test]
fn invalid_smiles_are_rejected() {
    let content = read_test_file("invalid-smiles.txt");

    let mut lines = 0usize;
    for line in content.lines() {
        lines += 1;
        // OpenBabel's `conv.Read` returns false (no molecule) for these; the
        // safe API surfaces that as either a parse error or an empty molecule.
        let rejected = match Molecule::parse(line, "smi") {
            Err(_) => true,
            Ok(mol) => mol.num_atoms() == 0,
        };
        assert!(
            rejected,
            "line {lines} parsed but should be invalid SMILES: {line:?}"
        );
    }

    assert!(lines > 0, "invalid-smiles.txt was empty");
}
