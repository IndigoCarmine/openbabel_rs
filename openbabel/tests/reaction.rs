//! Integration tests for SMARTS-based transformations (T7, OBChemTsfm).

use openbabel::{Molecule, Transform};

#[test]
fn deprotonate_carboxylic_acid() {
    // OpenBabel's own pH-model transform: a neutral carboxyl -OH -> -O(-).
    let t = Transform::new("O=C[OD1-0:1]", "O=C[O-:1]").expect("transform");
    let mut mol = Molecule::parse("CC(=O)O", "smi").expect("acetic acid");
    assert_eq!(mol.total_charge(), 0);

    assert!(t.apply(&mut mol), "transform should match acetic acid");
    assert_eq!(mol.total_charge(), -1, "carboxylate carries -1");
}

#[test]
fn protonate_amine_adds_positive_charge() {
    // The opposite direction: a neutral amine nitrogen gains a +1 charge.
    let t = Transform::new("[NX3:1]", "[NX3+:1]").expect("transform");
    let mut mol = Molecule::parse("CN", "smi").expect("methylamine");
    assert_eq!(mol.total_charge(), 0);

    assert!(t.apply(&mut mol), "should protonate the amine");
    assert_eq!(mol.total_charge(), 1, "ammonium carries +1");
}

#[test]
fn transform_without_match_is_a_noop() {
    let t = Transform::new("[NX3:1]", "[NX3+:1]").expect("transform");
    let mut mol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert!(!t.apply(&mut mol), "no amine nitrogen to transform");
    assert_eq!(mol.total_charge(), 0);
}

#[test]
fn invalid_transform_is_rejected() {
    assert!(Transform::new("[N:1", "[N+:1]").is_err());
}
