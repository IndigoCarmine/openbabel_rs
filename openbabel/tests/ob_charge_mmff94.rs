//! Port of OpenBabel's `test/charge_mmff94.cpp` — the MMFF94 charge model.
//!
//! Reads the molecules from `forcefield.sdf`, computes MMFF94 charges, and
//! compares each atom's partial charge against `charge-mmff94.txt` (one charge
//! per line, all atoms of all molecules in order), with the C++ tolerance of
//! 1e-3.
//!
//! The dipole-moment comparison against `dipole-mmff94.txt` is **skipped**:
//! `GetDipoleMoment` is not exposed by the binding.

mod common;

use common::{ob_test_file, read_test_file};
use openbabel::Molecule;

#[test]
fn mmff94_partial_charges_match_reference() {
    let path = ob_test_file("forcefield.sdf");
    let molecules = Molecule::read_file_many(path.to_str().unwrap(), Some("sdf"))
        .expect("read forcefield.sdf");
    let charges = read_test_file("charge-mmff94.txt");
    let mut ref_charges = charges.lines();

    let mut checked = 0usize;
    for mut mol in molecules {
        let title = mol.title();
        assert!(
            mol.compute_charges("mmff94"),
            "could not compute MMFF94 charges on {title:?}"
        );
        for atom in mol.atoms() {
            let expected: f64 = ref_charges
                .next()
                .expect("ran out of reference charges")
                .trim()
                .parse()
                .expect("charge value");
            let got = atom.partial_charge();
            assert!(
                (got - expected).abs() <= 1.0e-3,
                "MMFF94 charge for atom {} of {title:?}: expected {expected}, got {got}",
                atom.index()
            );
            checked += 1;
        }
    }

    assert!(checked > 0, "no atoms were checked");
}
