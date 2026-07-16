//! Port of OpenBabel's `test/ffgaff.cpp` — GAFF force-field single-point
//! energies.
//!
//! Reads the molecules from `gaff.sdf` (each with 3D coordinates) and compares
//! the GAFF energy against `gaffresults.txt` (one energy per molecule), with the
//! C++ tolerance of 1e-3.
//!
//! The `ValidateGradients()` check is **skipped** (gradient validation is not
//! exposed by the binding).

mod common;

use common::{ob_test_file, read_test_file};
use openbabel::Molecule;

#[test]
fn gaff_single_point_energies_match_reference() {
    let path = ob_test_file("gaff.sdf");
    let molecules =
        Molecule::read_file_many(path.to_str().unwrap(), Some("sdf")).expect("read gaff.sdf");
    let results = read_test_file("gaffresults.txt");
    let mut ref_energies = results.lines();

    let mut checked = 0usize;
    for mol in molecules {
        let expected: f64 = ref_energies
            .next()
            .expect("ran out of reference energies")
            .trim()
            .parse()
            .expect("energy value");
        let got = mol.energy("GAFF").expect("GAFF energy");
        assert!(
            (got - expected).abs() <= 1.0e-3,
            "GAFF energy for {:?}: expected {expected}, got {got}",
            mol.title()
        );
        checked += 1;
    }

    assert!(checked > 0, "no molecules were checked");
}
