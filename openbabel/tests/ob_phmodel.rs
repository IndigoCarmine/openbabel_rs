//! Port of OpenBabel's `test/phmodel.cpp` — pH-dependent (de)protonation.
//!
//! Each case parses an amino acid / peptide from SMILES, adds hydrogens for a
//! given pH (`AddHydrogens(false, true, pH)` == [`Molecule::add_hydrogens_for_ph`]),
//! and asserts the resulting atom count. The two cases the C++ marks as known
//! bugs (Glu at pH 4.15, Lys at pH 9.0) are commented out upstream and are
//! likewise omitted here.

mod common;

use openbabel::Molecule;

fn atoms_at_ph(smiles: &str, ph: f64) -> u32 {
    let mut mol = Molecule::parse(smiles, "smi").expect("parse SMILES");
    mol.add_hydrogens_for_ph(ph);
    mol.num_atoms()
}

#[test]
fn aspartic_acid() {
    let s = "NC(CC(O)=O)C(O)=O";
    assert_eq!(atoms_at_ph(s, 1.0), 17);
    assert_eq!(atoms_at_ph(s, 3.9), 16);
    assert_eq!(atoms_at_ph(s, 7.4), 15);
    assert_eq!(atoms_at_ph(s, 13.0), 14);
}

#[test]
fn glutamic_acid() {
    let s = "NC(CCC(O)=O)C(O)=O";
    assert_eq!(atoms_at_ph(s, 1.0), 20);
    assert_eq!(atoms_at_ph(s, 7.4), 18);
    assert_eq!(atoms_at_ph(s, 13.0), 17);
}

#[test]
fn histidine() {
    let s = "NC(Cc1nc[nH]c1)C(O)=O";
    assert_eq!(atoms_at_ph(s, 1.0), 22);
    assert_eq!(atoms_at_ph(s, 5.0), 21);
    assert_eq!(atoms_at_ph(s, 7.4), 20);
    assert_eq!(atoms_at_ph(s, 13.0), 19);
}

#[test]
fn lysine() {
    let s = "NC(CCCCN)C(O)=O";
    assert_eq!(atoms_at_ph(s, 1.0), 26);
    assert_eq!(atoms_at_ph(s, 7.4), 25);
    assert_eq!(atoms_at_ph(s, 13.0), 23);
}

#[test]
fn tyrosine() {
    let s = "NC(Cc1ccc(O)cc1)C(O)=O";
    assert_eq!(atoms_at_ph(s, 1.0), 25);
    assert_eq!(atoms_at_ph(s, 7.4), 24);
    assert_eq!(atoms_at_ph(s, 10.05), 23);
    assert_eq!(atoms_at_ph(s, 13.0), 22);
}

#[test]
fn arginine() {
    let s = "NC(CCCNC(N)=N)C(O)=O";
    assert_eq!(atoms_at_ph(s, 1.0), 28);
    assert_eq!(atoms_at_ph(s, 7.4), 27);
    assert_eq!(atoms_at_ph(s, 11.0), 26);
    assert_eq!(atoms_at_ph(s, 13.0), 25);
}

#[test]
fn gly_gly() {
    let s = "NCC(=O)NCC(=O)O";
    assert_eq!(atoms_at_ph(s, 1.0), 18);
    assert_eq!(atoms_at_ph(s, 7.4), 17);
    assert_eq!(atoms_at_ph(s, 13.0), 16);
}
