//! Port of OpenBabel's `test/ffmmff94.cpp` — MMFF94 / MMFF94s single-point
//! energies.
//!
//! Compares force-field energies against the reference files, with the C++
//! tolerance of 1e-3. Ports the four dielectric-constant-1.0 cases (the C++'s
//! choices 1–4):
//!   * MMFF94  on `forcefield.sdf`  vs `mmff94results.txt`
//!   * MMFF94  on `more-mmff94.sdf` vs `more-mmff94results.txt`
//!   * MMFF94s on `forcefield.sdf`  vs `mmff94sresults.txt`
//!   * MMFF94s on `more-mmff94.sdf` vs `more-mmff94sresults.txt`
//!
//! **Skipped** (per the "skip & document" policy): the `ValidateGradients()`
//! checks, and the C++'s choices 5–6, which call `SetDielectricConstant(4.0)`
//! (not exposed) — and choice 6 additionally references a results file that
//! does not exist in the OpenBabel tree.

mod common;

use common::{ob_test_file, read_test_file};
use openbabel::Molecule;

fn check_energies(sdf: &str, results: &str, method: &str) {
    let path = ob_test_file(sdf);
    let molecules =
        Molecule::read_file_many(path.to_str().unwrap(), Some("sdf")).expect("read sdf molecules");
    let content = read_test_file(results);
    let mut ref_energies = content.lines();

    let mut checked = 0usize;
    for mol in molecules {
        let expected: f64 = ref_energies
            .next()
            .expect("ran out of reference energies")
            .trim()
            .parse()
            .expect("energy value");
        let got = mol
            .energy(method)
            .unwrap_or_else(|| panic!("{method} energy failed for {:?}", mol.title()));
        assert!(
            (got - expected).abs() <= 1.0e-3,
            "{method} energy for {:?}: expected {expected}, got {got}",
            mol.title()
        );
        checked += 1;
    }

    assert!(checked > 0, "no molecules were checked for {method}/{sdf}");
}

#[test]
fn mmff94_energies_forcefield_sdf() {
    check_energies("forcefield.sdf", "mmff94results.txt", "MMFF94");
}

#[test]
fn mmff94_energies_more_sdf() {
    check_energies("more-mmff94.sdf", "more-mmff94results.txt", "MMFF94");
}

#[test]
fn mmff94s_energies_forcefield_sdf() {
    check_energies("forcefield.sdf", "mmff94sresults.txt", "MMFF94s");
}

#[test]
fn mmff94s_energies_more_sdf() {
    check_energies("more-mmff94.sdf", "more-mmff94sresults.txt", "MMFF94s");
}
