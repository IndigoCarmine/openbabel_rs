//! Integration tests for per-atom and per-bond string data (T19).

use openbabel::Molecule;

#[test]
fn atom_data_round_trips() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");

    // Absent by default.
    assert_eq!(mol.atom(0).unwrap().data("label"), None);

    mol.atom_mut(0).unwrap().set_data("label", "alpha carbon");
    assert_eq!(
        mol.atom(0).unwrap().data("label").as_deref(),
        Some("alpha carbon"),
    );
    // A different atom is unaffected.
    assert_eq!(mol.atom(1).unwrap().data("label"), None);

    // Setting again replaces the value.
    mol.atom_mut(0).unwrap().set_data("label", "changed");
    assert_eq!(mol.atom(0).unwrap().data("label").as_deref(), Some("changed"));
}

#[test]
fn bond_data_round_trips() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");

    assert_eq!(mol.bond(0).unwrap().data("role"), None);
    mol.bond_mut(0).unwrap().set_data("role", "scissile");
    assert_eq!(mol.bond(0).unwrap().data("role").as_deref(), Some("scissile"));
    assert_eq!(mol.bond(1).unwrap().data("role"), None);
}

#[test]
fn chained_atom_data_setters() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse");
    mol.atom_mut(2)
        .unwrap()
        .set_data("kind", "hydroxyl O")
        .set_data("charge_note", "acceptor");

    let o = mol.atom(2).unwrap();
    assert_eq!(o.data("kind").as_deref(), Some("hydroxyl O"));
    assert_eq!(o.data("charge_note").as_deref(), Some("acceptor"));
}
