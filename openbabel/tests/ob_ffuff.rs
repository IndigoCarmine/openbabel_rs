//! Port of OpenBabel's `test/ffuff.cpp` — UFF force-field single-point energies.
//!
//! Reads the molecules from `forcefield.sdf` (each with 3D coordinates) and
//! compares the UFF energy against `uffresults.txt` (one energy per molecule),
//! with the C++ tolerance of 1e-3.
//!
//! The C++ also calls `ValidateGradients()`; the binding does not expose
//! gradient validation, so that check is **skipped**.

mod common;

use common::{ob_test_file, read_test_file};
use openbabel::Molecule;

#[test]
fn uff_single_point_energies_match_reference() {
    let path = ob_test_file("forcefield.sdf");
    let molecules = Molecule::read_file_many(path.to_str().unwrap(), Some("sdf"))
        .expect("read forcefield.sdf");
    let results = read_test_file("uffresults.txt");
    let mut ref_energies = results.lines();

    let mut checked = 0usize;
    for mol in molecules {
        let expected: f64 = ref_energies
            .next()
            .expect("ran out of reference energies")
            .trim()
            .parse()
            .expect("energy value");
        let got = mol.energy("UFF").expect("UFF energy");
        assert!(
            (got - expected).abs() <= 1.0e-3,
            "UFF energy for {:?}: expected {expected}, got {got}",
            mol.title()
        );
        checked += 1;
    }

    assert!(checked > 0, "no molecules were checked");
}
