//! Integration tests for topological symmetry and canonical atom ranking (T13).

use openbabel::Molecule;
use std::collections::HashSet;

fn distinct(values: &[u32]) -> usize {
    values.iter().copied().collect::<HashSet<_>>().len()
}

#[test]
fn symmetry_classes_have_one_entry_per_atom() {
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    assert_eq!(mol.symmetry_classes().len() as u32, mol.num_atoms());
    assert_eq!(mol.canonical_ranks().len() as u32, mol.num_atoms());
}

#[test]
fn benzene_carbons_are_all_equivalent() {
    // All six aromatic carbons are related by symmetry → one symmetry class.
    let mol = Molecule::parse("c1ccccc1", "smi").expect("parse");
    let classes = mol.symmetry_classes();
    assert_eq!(classes.len(), 6);
    assert_eq!(distinct(&classes), 1, "classes: {classes:?}");
}

#[test]
fn propane_methyls_are_equivalent_middle_is_not() {
    // CCC: atoms 0 and 2 are the equivalent methyls; atom 1 is the middle CH2.
    let mol = Molecule::parse("CCC", "smi").expect("parse");
    let c = mol.symmetry_classes();
    assert_eq!(c.len(), 3);
    assert_eq!(c[0], c[2], "methyl carbons should match: {c:?}");
    assert_ne!(c[0], c[1], "middle carbon should differ: {c:?}");
    assert_eq!(distinct(&c), 2);
}

#[test]
fn ethanol_atoms_are_all_distinct() {
    // C–C–O: no two heavy atoms are equivalent.
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    assert_eq!(distinct(&mol.symmetry_classes()), 3);
}

#[test]
fn canonical_ranks_are_a_permutation() {
    // Even though benzene's carbons are all symmetry-equivalent, canonical
    // labelling assigns each a distinct rank: a permutation of 1..=6.
    let mol = Molecule::parse("c1ccccc1", "smi").expect("parse");
    let mut ranks = mol.canonical_ranks();
    assert_eq!(ranks.len(), 6);
    ranks.sort_unstable();
    assert_eq!(ranks, vec![1, 2, 3, 4, 5, 6]);
}

#[test]
fn canonical_ranks_are_input_order_independent() {
    // The same molecule written two ways must produce the same canonical
    // atom sequence: order atoms by canonical rank, then compare elements.
    fn elements_by_rank(smiles: &str) -> Vec<u32> {
        let mol = Molecule::parse(smiles, "smi").expect("parse");
        let ranks = mol.canonical_ranks();
        let mut pairs: Vec<(u32, u32)> = mol
            .atoms()
            .enumerate()
            .map(|(i, a)| (ranks[i], a.atomic_number()))
            .collect();
        pairs.sort_unstable_by_key(|&(rank, _)| rank);
        pairs.into_iter().map(|(_, z)| z).collect()
    }

    // Ethanol, atoms listed C,C,O vs O,C,C.
    assert_eq!(elements_by_rank("CCO"), elements_by_rank("OCC"));
    // 2-chloropropane written from either end.
    assert_eq!(elements_by_rank("CC(Cl)C"), elements_by_rank("C(Cl)(C)C"));
}
