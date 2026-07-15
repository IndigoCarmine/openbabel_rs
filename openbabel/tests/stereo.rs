//! Integration tests for stereochemistry perception (T7).

use openbabel::Molecule;

#[test]
fn alanine_has_one_tetrahedral_center() {
    // L-alanine: the second atom (index 1) is the chiral carbon.
    let mol = Molecule::parse("C[C@@H](N)C(=O)O", "smi").expect("alanine");
    assert_eq!(mol.tetrahedral_stereo_count(), 1);
    assert_eq!(mol.cistrans_stereo_count(), 0);

    let chiral = mol.atom(1).expect("atom 1");
    assert_eq!(chiral.atomic_number(), 6, "chiral center is carbon");
    assert!(chiral.is_tetrahedral_stereo());
    assert!(
        chiral.stereo_winding().is_some(),
        "a specified stereocenter has a winding"
    );
}

#[test]
fn achiral_molecule_has_no_stereocenters() {
    let mol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert_eq!(mol.tetrahedral_stereo_count(), 0);
    assert_eq!(mol.cistrans_stereo_count(), 0);
    assert!(!mol.atom(1).unwrap().is_tetrahedral_stereo());
}

#[test]
fn trans_difluoroethene_has_a_cistrans_bond() {
    let mol = Molecule::parse("F/C=C/F", "smi").expect("trans-1,2-difluoroethene");
    assert_eq!(mol.cistrans_stereo_count(), 1);
    assert_eq!(mol.tetrahedral_stereo_count(), 0);

    let double_bond_has_stereo = mol
        .bonds()
        .any(|b| b.order() == 2 && b.is_cistrans_stereo());
    assert!(double_bond_has_stereo, "the C=C bond should carry stereo");
}
