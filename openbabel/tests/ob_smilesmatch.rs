//! Port of OpenBabel's `test/smilesmatch.cpp` — a SMILES, used as a SMARTS
//! query, matches the molecule it came from.
//!
//! For each molecule in `nci.smi`: strip the title, skip disconnected SMILES
//! (those containing `.`, as the C++ does), then compile the SMILES as a SMARTS
//! pattern and assert it matches its own molecule.

mod common;

use common::read_test_file;
use openbabel::{Molecule, SmartsPattern};

#[test]
fn smiles_match_themselves_as_smarts() {
    let content = read_test_file("nci.smi");

    let mut tested = 0usize;
    for line in content.lines() {
        let mol = match Molecule::parse(line, "smi") {
            Ok(m) if m.num_atoms() > 0 => m,
            _ => continue,
        };

        // Trim off any title (everything from the first whitespace).
        let smiles = line
            .split([' ', '\t', '\r', '\n'])
            .next()
            .unwrap_or("");
        // Skip disconnected structures, exactly as the C++ does.
        if smiles.contains('.') {
            continue;
        }

        let pattern = SmartsPattern::new(smiles)
            .unwrap_or_else(|_| panic!("could not compile SMILES as SMARTS: {smiles}"));
        assert!(
            pattern.matches(&mol),
            "SMARTS did not match its own SMILES molecule: {smiles}"
        );
        tested += 1;
    }

    assert!(tested > 0, "no molecules were tested from nci.smi");
}
