//! Port of OpenBabel's `test/aromatest.cpp` — aromaticity perception.
//!
//! For every molecule in `aromatics.smi` the C++ test makes hydrogens explicit,
//! then removes them, and after each step asserts that *every heavy atom* is
//! aromatic. It then checks a handful of hand-picked negative cases where not
//! every atom should be perceived as aromatic.

mod common;

use common::read_test_file;
use openbabel::Molecule;

/// Every heavy atom of every molecule in `aromatics.smi` is aromatic, both with
/// explicit hydrogens added and with them removed again.
#[test]
fn aromatics_file_all_heavy_atoms_aromatic() {
    let content = read_test_file("aromatics.smi");

    let mut molecules = 0usize;
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let mut mol = match Molecule::parse(line, "smi") {
            Ok(m) if m.num_atoms() > 0 => m,
            _ => continue,
        };
        molecules += 1;
        let title = mol.title();

        for pass in 0..2 {
            if pass == 0 {
                mol.add_hydrogens();
            } else {
                mol.remove_hydrogens();
            }
            for atom in mol.atoms() {
                if atom.atomic_number() == 1 {
                    continue;
                }
                assert!(
                    atom.is_aromatic(),
                    "atom {} is not aromatic in molecule {title:?} (pass {pass})",
                    atom.index()
                );
            }
        }
    }

    assert!(molecules > 0, "no molecules were read from aromatics.smi");
}

/// Molecules that were aromatized in error during earlier development: not every
/// heavy atom should be aromatic.
#[test]
fn negative_cases_not_all_aromatic() {
    let cases = [
        "c1ccc2[N+]=c3ccccc3=Nc2c1", // N radical found in eMolecules
        "N1S[SH+]C=C1",
        "S1C=[NH+]=[NH+]=C1",
        "C1(N23)=CC=CC2=CC=CC3=CC=C1", // pyrido[2,1,6-de]quinolizine
    ];
    for smiles in cases {
        let mol = Molecule::parse(smiles, "smi").expect("parse negative case");
        let all_aromatic = mol
            .atoms()
            .filter(|a| a.atomic_number() != 1)
            .all(|a| a.is_aromatic());
        assert!(
            !all_aromatic,
            "every atom was aromatic (unexpected) in {smiles}"
        );
    }
}
