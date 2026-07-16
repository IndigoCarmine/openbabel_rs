//! Port of OpenBabel's `test/graphsymtest.cpp` — topological symmetry classes.
//!
//! `OBGraphSym::GetSymmetry` maps to [`Molecule::symmetry_classes`].
//!   * Parts 1–2 (`genericGraphSymTest`): the number of distinct symmetry
//!     classes is invariant across a canonical-SMILES roundtrip, and atoms that
//!     share a class value across the two molecules agree on element / degrees /
//!     implicit-H count.
//!   * Parts 3–5 (`countGraphSymClassesTest`): the number of distinct symmetry
//!     classes for each `stereo/razinger_*.mol` equals a known value. (Only the
//!     files the C++ actually references — the ones it comments out as "missing"
//!     are omitted here too.)

mod common;

use common::ob_test_file;
use openbabel::Molecule;

fn distinct_classes(classes: &[u32]) -> usize {
    let mut v = classes.to_vec();
    v.sort_unstable();
    v.dedup();
    v.len()
}

fn generic_graphsym_test(smiles: &str) {
    let mol1 = Molecule::parse(smiles, "smi").expect("parse smiles");
    let sym1 = mol1.symmetry_classes();

    let canonical = mol1.write("can").expect("write canonical smiles");
    let mol2 = Molecule::parse(canonical.trim(), "smi").expect("reparse canonical");
    let sym2 = mol2.symmetry_classes();

    assert_eq!(
        distinct_classes(&sym1),
        distinct_classes(&sym2),
        "distinct symmetry-class count changed across roundtrip for {smiles}"
    );

    // For every atom of mol1, the first atom of mol2 sharing its symmetry-class
    // value must agree on the basic graph invariants.
    for a1 in mol1.atoms() {
        let class1 = sym1[a1.index() as usize];
        if let Some(i2) = sym2.iter().position(|&c| c == class1) {
            let a2 = mol2.atom(i2 as u32).unwrap();
            assert_eq!(a1.atomic_number(), a2.atomic_number(), "atomic number ({smiles})");
            assert_eq!(a1.degree(), a2.degree(), "degree ({smiles})");
            assert_eq!(a1.heavy_degree(), a2.heavy_degree(), "heavy degree ({smiles})");
            assert_eq!(a1.hetero_degree(), a2.hetero_degree(), "hetero degree ({smiles})");
            assert_eq!(
                a1.implicit_hydrogens(),
                a2.implicit_hydrogens(),
                "implicit H count ({smiles})"
            );
        }
    }
}

fn count_classes_test(file: &str, expected: usize) {
    let path = ob_test_file(file);
    let mol = Molecule::read_file(path.to_str().unwrap(), Some("mol"))
        .unwrap_or_else(|_| panic!("read {file}"));
    assert_eq!(
        distinct_classes(&mol.symmetry_classes()),
        expected,
        "distinct symmetry classes for {file}"
    );
}

/// Part 1.
#[test]
fn generic_smiles_roundtrip_symmetry() {
    for smiles in [
        "C[C@H](O)N",
        "Cl[C@@](CCl)(I)Br",
        "Cl/C=C/F",
        r"CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO",
        "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5",
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1",
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2",
        r"CC(=O)OCCC(/C)=C\C[C@H](C(C)=C)CCC=C",
        "CC[C@H](O1)CC[C@@]12CCCO2",
        "CN1CCC[C@H]1c2cccnc2",
        "C(CS[14CH2][14C@@H]1[14C@H]([14C@H]([14CH](O1)O)O)O)[C@@H](C(=O)O)N",
        "CCC[C@@H]1C[C@H](N(C1)C)C(=O)NC([C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)SC)OP(=O)(O)O)O)O)C(C)Cl",
        "CC(C)[C@H]1CC[C@]([C@@H]2[C@@H]1C=C(COC2=O)C(=O)O)(CCl)O",
        "CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2",
    ] {
        generic_graphsym_test(smiles);
    }
}

/// Part 2.
#[test]
fn generic_smiles_roundtrip_symmetry_aromatic() {
    generic_graphsym_test("Cc1cn(c(=O)[nH]c1=O)[C@H]1C[C@@H]([C@H](O1)CNCc1ccccc1)O");
}

/// Part 3.
#[test]
fn razinger_class_counts_part3() {
    let cases = [
        ("stereo/razinger_fig3.mol", 4),
        ("stereo/razinger_fig7_1.mol", 5),
        ("stereo/razinger_fig7_2.mol", 7),
        ("stereo/razinger_fig7_3.mol", 7),
        ("stereo/razinger_fig7_4.mol", 7),
        ("stereo/razinger_fig7_5.mol", 5),
        ("stereo/razinger_fig7_6.mol", 4),
        ("stereo/razinger_fig7_7.mol", 7),
        ("stereo/razinger_fig7_8.mol", 7),
        ("stereo/razinger_fig7_9.mol", 10),
        ("stereo/razinger_fig7_10.mol", 13),
        ("stereo/razinger_fig7_11.mol", 9),
        ("stereo/razinger_fig7_12.mol", 17),
        ("stereo/razinger_fig7_13.mol", 4),
        ("stereo/razinger_fig7_14.mol", 4),
        ("stereo/razinger_fig7_15.mol", 5),
        ("stereo/razinger_fig7_16.mol", 4),
        ("stereo/razinger_fig7_17.mol", 5),
        ("stereo/razinger_fig7_18.mol", 3),
        ("stereo/razinger_fig7_19.mol", 6),
    ];
    for (file, n) in cases {
        count_classes_test(file, n);
    }
}

/// Part 4.
#[test]
fn razinger_class_counts_part4() {
    let cases = [
        ("stereo/razinger_fig7_20.mol", 5),
        ("stereo/razinger_fig7_21.mol", 7),
        ("stereo/razinger_fig7_22.mol", 2),
        ("stereo/razinger_fig7_23.mol", 4),
        ("stereo/razinger_fig7_24.mol", 4),
        ("stereo/razinger_fig7_25.mol", 4),
        ("stereo/razinger_fig7_26.mol", 4),
        ("stereo/razinger_fig7_27.mol", 8),
        ("stereo/razinger_fig7_28.mol", 3),
        ("stereo/razinger_fig7_29.mol", 6),
    ];
    for (file, n) in cases {
        count_classes_test(file, n);
    }
}

/// Part 5.
#[test]
fn razinger_class_counts_part5() {
    let cases = [
        ("stereo/razinger_fig7_48.mol", 3),
        ("stereo/razinger_fig7_49.mol", 11),
        ("stereo/razinger_fig7_50.mol", 5),
        ("stereo/razinger_fig7_51.mol", 7),
        ("stereo/razinger_fig7_52.mol", 4),
        ("stereo/razinger_fig7_58.mol", 5),
        ("stereo/razinger_fig7_59.mol", 4),
        ("stereo/razinger_fig7_60.mol", 7),
        ("stereo/razinger_fig7_64.mol", 8),
        ("stereo/razinger_fig7_69.mol", 2),
    ];
    for (file, n) in cases {
        count_classes_test(file, n);
    }
}
