//! Port of OpenBabel's `test/iterators.cpp` — molecule iterators.
//!
//! For every molecule in `attype.00.smi`, the atom-, bond- and ring-iterator
//! counts must equal `num_atoms()`, `num_bonds()` and the SSSR size
//! respectively (the ring-iterator check covers regression PR#2815025).
//!
//! The C++ also exercises `FOR_DFS_OF_MOL` / `FOR_BFS_OF_MOL`; the safe API
//! exposes no depth-/breadth-first atom iterator, so those two checks are
//! **skipped** (per the "skip & document" policy).

mod common;

use common::read_test_file;
use openbabel::Molecule;

#[test]
fn atom_bond_and_ring_iterator_counts() {
    let content = read_test_file("attype.00.smi");

    let mut molecules = 0usize;
    for line in content.lines() {
        let mol = match Molecule::parse(line, "smi") {
            Ok(m) if m.num_atoms() > 0 => m,
            _ => continue,
        };
        molecules += 1;
        let title = mol.title();

        assert_eq!(
            mol.atoms().count() as u32,
            mol.num_atoms(),
            "atom iterator count for {title:?}"
        );
        assert_eq!(
            mol.bonds().count() as u32,
            mol.num_bonds(),
            "bond iterator count for {title:?}"
        );
        assert_eq!(
            mol.rings().count() as u32,
            mol.num_rings(),
            "ring iterator count for {title:?}"
        );
    }

    assert!(molecules > 0, "no molecules were read from attype.00.smi");
}
