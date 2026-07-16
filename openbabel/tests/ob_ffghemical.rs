//! Port of OpenBabel's `test/ffghemical.cpp` — Ghemical force-field single-point
//! energies.
//!
//! Reads the molecules from `forcefield.sdf` (each with 3D coordinates) and
//! compares the Ghemical energy against `ghemicalresults.txt` (one energy per
//! molecule), with the C++ tolerance of 1e-3.
//!
//! The `ValidateGradients()` check is **skipped** (gradient validation is not
//! exposed by the binding).

mod common;

use common::{ob_test_file, read_test_file};
use openbabel::Molecule;

#[test]
fn ghemical_single_point_energies_match_reference() {
    let path = ob_test_file("forcefield.sdf");
    let molecules = Molecule::read_file_many(path.to_str().unwrap(), Some("sdf"))
        .expect("read forcefield.sdf");
    let results = read_test_file("ghemicalresults.txt");
    let mut ref_energies = results.lines();

    let mut checked = 0usize;
    for mol in molecules {
        let expected: f64 = ref_energies
            .next()
            .expect("ran out of reference energies")
            .trim()
            .parse()
            .expect("energy value");
        let got = mol.energy("Ghemical").expect("Ghemical energy");
        assert!(
            (got - expected).abs() <= 1.0e-3,
            "Ghemical energy for {:?}: expected {expected}, got {got}",
            mol.title()
        );
        checked += 1;
    }

    assert!(checked > 0, "no molecules were checked");
}
