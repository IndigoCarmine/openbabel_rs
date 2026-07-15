//! Integration tests for the safe `openbabel` API, mirroring the kind of
//! checks OpenBabel's own suite makes on core molecule handling.

use openbabel::Molecule;

/// Ethanol from SMILES: formula, weight, atom/bond counts, canonical roundtrip.
#[test]
fn ethanol_from_smiles() {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse CCO");

    // Before adding H: 3 heavy atoms (C, C, O), 2 bonds.
    assert_eq!(mol.num_atoms(), 3, "heavy-atom count");
    assert_eq!(mol.num_bonds(), 2, "heavy-atom bond count");

    // Formula reflects implicit hydrogens even before they are made explicit.
    assert_eq!(mol.formula(), "C2H6O");

    // Molecular weight of ethanol ~ 46.07 g/mol.
    let mw = mol.molar_mass();
    assert!((mw - 46.07).abs() < 0.05, "unexpected molar mass: {mw}");

    // Neutral molecule.
    assert_eq!(mol.total_charge(), 0);

    // After making hydrogens explicit: 9 atoms (3 heavy + 6 H).
    mol.add_hydrogens();
    assert_eq!(mol.num_atoms(), 9, "atom count after AddHydrogens");

    // Canonical SMILES should roundtrip to an equivalent molecule.
    let can = mol.write("can").expect("write canonical smiles");
    assert!(!can.trim().is_empty(), "canonical SMILES was empty");
    let mol2 = Molecule::parse(can.trim(), "smi").expect("reparse canonical");
    assert_eq!(mol2.formula(), "C2H6O", "roundtrip formula mismatch");
}

/// Atom-level access: elements and aromaticity of benzene.
#[test]
fn benzene_atoms_and_aromaticity() {
    let mol = Molecule::parse("c1ccccc1", "smi").expect("parse benzene");
    assert_eq!(mol.num_atoms(), 6);
    assert_eq!(mol.num_bonds(), 6);
    assert_eq!(mol.formula(), "C6H6");

    // Every ring atom is aromatic carbon.
    for atom in mol.atoms() {
        assert_eq!(atom.atomic_number(), 6, "expected carbon");
        assert!(atom.is_aromatic(), "benzene carbons are aromatic");
        assert!(atom.is_in_ring(), "benzene carbons are in a ring");
    }
    // Every bond is aromatic and in a ring.
    for bond in mol.bonds() {
        assert!(bond.is_aromatic(), "benzene bonds are aromatic");
        assert!(bond.is_in_ring());
    }
}

/// An unknown format id is a clean error, not a panic.
#[test]
fn unknown_format_errors() {
    let err = Molecule::parse("CCO", "not-a-format");
    assert!(err.is_err(), "unknown input format should error");
}

/// Garbage SMILES fails to parse rather than producing a bogus molecule.
#[test]
fn invalid_smiles_errors() {
    let err = Molecule::parse("this is not smiles !!!", "smi");
    assert!(err.is_err(), "invalid SMILES should error");
}
