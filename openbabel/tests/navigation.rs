//! Integration tests for graph navigation and crystallographic unit cells (T12).

use openbabel::{LatticeType, Molecule};

#[test]
fn atom_neighbors_and_bonds() {
    // Benzene: each aromatic carbon has 3 neighbours (2 C + 1 H) and 3 bonds
    // once the implicit hydrogens are made explicit.
    let mut mol = Molecule::parse("c1ccccc1", "smi").expect("parse");
    mol.add_hydrogens();
    let c0 = mol.atom(0).expect("atom 0");
    assert_eq!(c0.atomic_number(), 6);
    assert_eq!(c0.neighbors().len(), 3);
    assert_eq!(c0.bonds().len(), 3);
    assert_eq!(c0.degree(), 3);

    // Two of the neighbours are carbons, one is hydrogen.
    let mut carbons = 0;
    let mut hydrogens = 0;
    for n in c0.neighbors() {
        match n.atomic_number() {
            6 => carbons += 1,
            1 => hydrogens += 1,
            _ => {}
        }
    }
    assert_eq!((carbons, hydrogens), (2, 1));
}

#[test]
fn bond_between_and_other_atom() {
    let mol = Molecule::parse("CCO", "smi").expect("parse"); // C0–C1–O2
    // C0 and C1 are bonded; C0 and O2 are not.
    let b = mol.bond_between(0, 1).expect("C0–C1 bond");
    assert!(mol.bond_between(0, 2).is_none());

    // other_atom walks across the bond.
    let c0 = mol.atom(0).unwrap();
    let across = b.other_atom(&c0).expect("other end");
    assert_eq!(across.index(), 1);

    // bond_to mirrors bond_between.
    let c1 = mol.atom(1).unwrap();
    assert!(c0.bond_to(&c1).is_some());
    let o2 = mol.atom(2).unwrap();
    assert!(c0.bond_to(&o2).is_none());
}

#[test]
fn bond_orders_and_explicit_hydrogens() {
    // Ethene C=C: each carbon has one double bond and (after adding H) two H.
    let mut mol = Molecule::parse("C=C", "smi").expect("parse");
    let c0 = mol.atom(0).unwrap();
    assert_eq!(c0.count_bonds_of_order(2), 1);
    assert_eq!(c0.count_bonds_of_order(1), 0);

    mol.add_hydrogens();
    assert_eq!(mol.atom(0).unwrap().explicit_hydrogen_count(), 2);
}

/// A minimal orthorhombic CIF (a=5, b=6, c=7, all angles 90°), space group P 1.
const CELL_CIF: &str = r#"data_test
_cell_length_a 5.0000
_cell_length_b 6.0000
_cell_length_c 7.0000
_cell_angle_alpha 90.0000
_cell_angle_beta 90.0000
_cell_angle_gamma 90.0000
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0000 0.0000 0.0000
O1 O 0.5000 0.5000 0.5000
"#;

#[test]
fn unit_cell_parameters() {
    let mol = Molecule::parse(CELL_CIF, "cif").expect("parse CIF");
    assert!(mol.has_unit_cell());
    let cell = mol.unit_cell().expect("unit cell");

    let (a, b, c) = cell.lengths();
    assert!((a - 5.0).abs() < 1e-3 && (b - 6.0).abs() < 1e-3 && (c - 7.0).abs() < 1e-3);
    let (alpha, beta, gamma) = cell.angles();
    assert!((alpha - 90.0).abs() < 1e-3 && (beta - 90.0).abs() < 1e-3 && (gamma - 90.0).abs() < 1e-3);
    assert!((cell.volume() - 210.0).abs() < 0.1, "volume {}", cell.volume());
    assert_ne!(cell.lattice_type(), LatticeType::Undefined);
    assert!(!cell.space_group().is_empty());
}

#[test]
fn unit_cell_coordinate_transforms() {
    let mol = Molecule::parse(CELL_CIF, "cif").expect("parse CIF");
    let cell = mol.unit_cell().expect("unit cell");

    // Fractional (0.5, 0.5, 0.5) is the box centre: Cartesian (2.5, 3.0, 3.5).
    let (cx, cy, cz) = cell.to_cartesian(0.5, 0.5, 0.5);
    assert!((cx - 2.5).abs() < 1e-3 && (cy - 3.0).abs() < 1e-3 && (cz - 3.5).abs() < 1e-3);

    // Round-trip back to fractional.
    let (fx, fy, fz) = cell.to_fractional(cx, cy, cz);
    assert!((fx - 0.5).abs() < 1e-6 && (fy - 0.5).abs() < 1e-6 && (fz - 0.5).abs() < 1e-6);
}

#[test]
fn no_unit_cell_for_smiles() {
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    assert!(!mol.has_unit_cell());
    assert!(mol.unit_cell().is_none());
}
