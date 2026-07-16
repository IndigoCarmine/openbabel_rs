//! Port of OpenBabel's `test/canonconsistenttest.cpp` — canonical SMILES
//! round-trip stability.
//!
//! For every molecule in an SDF, writing canonical SMILES, re-reading it, and
//! writing canonical SMILES again must reproduce the identical string (with the
//! title cleared first, as the C++ does). This is the faithfully portable slice
//! of OpenBabel's canonical-stability suite — the shuffle-based tests require
//! `RenumberAtoms`, which the binding does not expose.

mod common;

use common::ob_test_file;
use openbabel::Molecule;

fn check_file(file: &str) {
    let path = ob_test_file(file);
    let molecules =
        Molecule::read_file_many(path.to_str().unwrap(), Some("sdf")).unwrap_or_else(|_| panic!("read {file}"));

    let mut checked = 0usize;
    for mut mol in molecules {
        mol.set_title("");
        let output = mol.write("can").expect("write canonical SMILES");
        let output = output.trim();

        let mut round2 = Molecule::parse(output, "smi")
            .unwrap_or_else(|_| panic!("re-read canonical SMILES {output:?}"));
        round2.set_title("");
        let roundtrip = round2.write("can").expect("write canonical SMILES");
        let roundtrip = roundtrip.trim();

        assert_eq!(
            output, roundtrip,
            "canonical SMILES not stable across roundtrip in {file}"
        );
        checked += 1;
    }

    assert!(checked > 0, "no molecules were read from {file}");
}

#[test]
fn forcefield_sdf() {
    check_file("forcefield.sdf");
}

#[test]
fn filterset_sdf() {
    check_file("filterset.sdf");
}

#[test]
fn cantest_sdf() {
    check_file("cantest.sdf");
}
