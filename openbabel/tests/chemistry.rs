//! Integration tests for the T4 surface: partial charges and richer atom
//! chemistry. (InChI identifiers are covered in `identity.rs`.)

use openbabel::Molecule;

#[test]
fn gasteiger_partial_charges_sum_to_total() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse ethanol");
    mol.add_hydrogens();
    assert!(mol.compute_charges("gasteiger"), "gasteiger charges");

    // Partial charges of a neutral molecule sum to ~0.
    let sum: f64 = mol.atoms().map(|a| a.partial_charge()).sum();
    assert!(sum.abs() < 0.05, "partial charges should sum to ~0, got {sum}");

    // The hydroxyl oxygen should carry a negative partial charge.
    let min = mol
        .atoms()
        .map(|a| a.partial_charge())
        .fold(f64::INFINITY, f64::min);
    assert!(min < -0.2, "expected a clearly negative atom, min was {min}");
}

#[test]
fn unknown_charge_model_fails_cleanly() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    assert!(!mol.compute_charges("no-such-model"));
}

#[test]
fn atom_chemistry_of_ethanol() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse ethanol"); // C-C-O
    // Before adding H, the terminal methyl carbon has 3 implicit H and degree 1.
    let c0 = mol.atom(0).expect("atom 0");
    assert_eq!(c0.atomic_number(), 6);
    assert_eq!(c0.degree(), 1, "methyl C connects to one heavy atom");
    assert_eq!(c0.implicit_hydrogens(), 3, "methyl C has 3 implicit H");
    assert_eq!(c0.total_valence(), 4, "carbon is tetravalent");
    assert_eq!(c0.hybridization(), 3, "sp3 carbon");

    // The hydroxyl oxygen accepts hydrogen bonds.
    let o = mol.atom(2).expect("atom 2");
    assert_eq!(o.atomic_number(), 8);
    assert!(o.is_hbond_acceptor(), "hydroxyl O accepts H-bonds");

    // After making H explicit, degree reflects the added hydrogens, and the
    // hydroxyl oxygen (now with an explicit H neighbour) is a donor.
    mol.add_hydrogens();
    let c0 = mol.atom(0).expect("atom 0");
    assert_eq!(c0.degree(), 4, "methyl C now has 4 explicit neighbours");
    assert_eq!(c0.implicit_hydrogens(), 0);
    let o = mol.atom(2).expect("atom 2");
    assert!(o.is_hbond_donor(), "hydroxyl O donates H-bonds");
}

#[test]
fn aromatic_carbon_is_sp2() {
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let c = mol.atom(0).expect("atom 0");
    assert!(c.is_aromatic());
    assert_eq!(c.hybridization(), 2, "aromatic carbon is sp2");
}
