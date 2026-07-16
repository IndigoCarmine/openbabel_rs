//! Port of OpenBabel's `test/pdbreadfile.cpp` — reading a ligand from PDB.
//!
//! All four files are the same 00T ligand (ideal / nonstandard, with and
//! without the `HETATM`-only variant), so the checks are identical: 22 atoms,
//! and specific elements at fixed positions. OpenBabel's `GetAtom(n)` is 1-based;
//! the Rust `atom(i)` accessor is 0-based, so each index is `n - 1`.

mod common;

use common::ob_test_file;
use openbabel::Molecule;

fn check(file: &str) {
    let path = ob_test_file(file);
    let mol = Molecule::read_file(path.to_str().unwrap(), Some("pdb"))
        .unwrap_or_else(|_| panic!("read {file}"));

    assert_eq!(mol.num_atoms(), 22, "{file}: atom count");
    assert_eq!(mol.atom(9).unwrap().atomic_number(), 17, "{file}: atom 10 is Cl");
    assert_eq!(mol.atom(5).unwrap().atomic_number(), 7, "{file}: atom 6 is N");
    assert_eq!(mol.atom(11).unwrap().atomic_number(), 1, "{file}: atom 12 is H");
    assert_eq!(mol.atom(12).unwrap().atomic_number(), 1, "{file}: atom 13 is H");
    assert_eq!(mol.atom(13).unwrap().atomic_number(), 1, "{file}: atom 14 is H");
}

#[test]
fn pdb_00t_ideal() {
    check("00T_ideal.pdb");
}

#[test]
fn pdb_00t_nonstandard() {
    check("00T_nonstandard.pdb");
}

#[test]
fn pdb_00t_ideal_het() {
    check("00T_ideal_het.pdb");
}

#[test]
fn pdb_00t_nonstandard_het() {
    check("00T_nonstandard_het.pdb");
}
