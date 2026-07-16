//! Port of OpenBabel's `test/mol.cpp` — unit tests for the core molecule type.
//!
//! The C++ file runs 15 numbered checks. This port reproduces each one that the
//! safe Rust API can express, keeping the same construction steps, inputs and
//! expected values.
//!
//! Two of the original checks exercise C++ API that the Rust binding does not
//! expose and are therefore **skipped** (per the "skip & document" policy):
//!   * ok 3  — `OBMol::ReserveAtoms(-1/0/2)`: no `reserve_atoms` binding; it is
//!     purely a capacity hint with no observable effect.
//!   * ok 13 — `OBMol::SetFormula`/`GetFormula` round-trip: there is no
//!     `set_formula` binding (`Molecule::formula` is derived from the atoms).

mod common;

use common::ob_test_file;
use openbabel::Molecule;

/// ok 1, 2, 4 — an empty molecule constructs and has no atoms.
/// (ok 3, `ReserveAtoms`, is skipped: no binding.)
#[test]
fn new_molecule_is_empty() {
    let mol = Molecule::new();
    assert_eq!(mol.num_atoms(), 0);
}

/// ok 5, 6 — adding an atom then a bond updates the counts.
#[test]
fn new_atom_and_bond() {
    let mut mol = Molecule::new();
    mol.add_atom(0); // OBMol::NewAtom() — atomic number defaults to 0
    assert_eq!(mol.num_atoms(), 1);

    mol.add_atom(0);
    // C++ `AddBond(1, 2, 1)` uses 1-based atom ids; the safe API is 0-based.
    mol.add_bond(0, 1, 1);
    assert_eq!(mol.num_bonds(), 1);
}

/// ok 7 — `Clear` empties the molecule.
#[test]
fn clear_empties_molecule() {
    let mut mol = Molecule::new();
    mol.add_atom(6);
    mol.add_atom(6);
    mol.add_bond(0, 1, 1);
    mol.clear();
    assert_eq!(mol.num_atoms(), 0);
}

/// ok 8 — read a 3D structure from `files/test3d.xyz`, then re-centre it.
#[test]
fn read_3d_xyz() {
    let path = ob_test_file("test3d.xyz");
    let mut mol = Molecule::read_file(path.to_str().unwrap(), Some("xyz"))
        .expect("read files/test3d.xyz");
    assert!(mol.num_atoms() > 0, "3D molecule should have atoms");
    mol.center(); // OBMol::Center()
}

/// ok 9 — low-level bond insertion inside a modify block must not crash
/// (regression PR#1665649). The C++ builds two carbons and joins them.
#[test]
fn bond_insertion_in_modify_block() {
    let mut mol = Molecule::new();
    mol.begin_modify();
    let a1 = mol.add_atom(6);
    mol.atom_mut(a1).unwrap().set_position(0.0, 0.0, 0.0);
    let a2 = mol.add_atom(6);
    mol.atom_mut(a2).unwrap().set_position(1.6, 0.0, 0.0);
    mol.add_bond(a1, a2, 1);
    mol.end_modify();
    assert_eq!(mol.num_atoms(), 2);
    assert_eq!(mol.num_bonds(), 1);
}

/// ok 10 — `AddHydrogens` fills in the implicit H count set inside a modify
/// block: one carbon with 4 implicit H becomes CH4 (5 atoms).
#[test]
fn add_hydrogens_from_implicit_count() {
    let mut mol = Molecule::new();
    mol.begin_modify();
    let c = mol.add_atom(6);
    {
        let mut a = mol.atom_mut(c).unwrap();
        a.set_position(0.5, 0.5, 0.5);
        a.set_implicit_hydrogens(4);
    }
    mol.end_modify();
    mol.add_hydrogens();
    assert_eq!(mol.num_atoms(), 5);
}

/// ok 11 — same as above but without a modify block (regression PR#1665519).
#[test]
fn add_hydrogens_without_modify_block() {
    let mut mol = Molecule::new();
    let c = mol.add_atom(6);
    {
        let mut a = mol.atom_mut(c).unwrap();
        a.set_position(0.5, 0.5, 0.5);
        a.set_implicit_hydrogens(4);
    }
    mol.add_hydrogens();
    assert_eq!(mol.num_atoms(), 5);
}

/// ok 12 — writing an empty molecule to InChI succeeds (regression PR#2864334).
#[test]
fn write_empty_inchi() {
    let mol = Molecule::new();
    assert!(
        mol.write("inchi").is_ok(),
        "writing an empty InChI must not fail"
    );
}

/// ok 14 — molecular formula with large / undefined atomic numbers: element 118
/// (Og), an undefined element 200 (ignored), and a deuterium (H, isotope 2)
/// give the formula `DOg`.
#[test]
fn formula_with_large_atomic_numbers() {
    let mut mol = Molecule::new();
    mol.begin_modify();
    mol.add_atom(118);
    mol.add_atom(200); // undefined — ignored, not a crash
    let h = mol.add_atom(1);
    mol.atom_mut(h).unwrap().set_isotope(2);
    mol.end_modify();
    assert_eq!(mol.formula(), "DOg");
}

/// ok 15 — the dihedral of four coplanar points spanning a straight backbone is
/// ±180°. The C++ calls the free `CalcTorsionAngle`; here we place the same four
/// points as atoms and read `Molecule::torsion`, which computes the identical
/// geometric dihedral.
#[test]
fn torsion_angle_of_planar_points() {
    let pts = [
        (-1.0, -1.0, 0.0),
        (-1.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
    ];
    let mut mol = Molecule::new();
    for (x, y, z) in pts {
        let i = mol.add_atom(6);
        mol.atom_mut(i).unwrap().set_position(x, y, z);
    }
    let dihedral = mol.torsion(0, 1, 2, 3);
    assert!(
        (dihedral.abs() - 180.0).abs() < 0.001,
        "expected ±180°, got {dihedral}"
    );
}
