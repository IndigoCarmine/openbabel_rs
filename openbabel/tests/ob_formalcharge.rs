//! Port of OpenBabel's `test/formalcharge.cpp` — molecular formal charges.
//!
//! **This is a faithful port of a no-op.** OpenBabel's `formalcharge.cpp` reads
//! `attype.00.smi` and `formalchargeresults.txt` in lockstep, tokenizes each
//! reference line, and then — where a charge comparison was clearly intended
//! (`// check charges`) — does nothing, always reporting a single passing test.
//!
//! Per the "replicate literally" fidelity rule, this port reproduces that exact
//! behaviour: it advances both files together (so a truncated reference would
//! still surface) but makes no assertion about the charges themselves.

mod common;

use common::{read_test_file, tokenize};
use openbabel::Molecule;

#[test]
fn formalcharge_is_a_literal_noop() {
    let smiles = read_test_file("attype.00.smi");
    let results = read_test_file("formalchargeresults.txt");
    let mut ref_lines = results.lines();

    let mut molecules = 0usize;
    for line in smiles.lines() {
        // Skip anything that does not parse to a non-empty molecule, exactly as
        // OpenBabel's `if (mol.Empty()) continue;` does (the molecule itself is
        // otherwise unused — the charge check was never implemented).
        match Molecule::parse(line, "smi") {
            Ok(m) if m.num_atoms() > 0 => {}
            _ => continue,
        }
        molecules += 1;
        let ref_line = ref_lines
            .next()
            .expect("ran out of reference data (formal charge)");
        // Tokenized exactly as the C++ does, then — exactly as the C++ does —
        // left unchecked.
        let _ = tokenize(ref_line);
    }

    assert!(molecules > 0, "no molecules were read from attype.00.smi");
}
