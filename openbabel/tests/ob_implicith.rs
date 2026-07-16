//! Port of OpenBabel's `test/implicitHtest.cpp` — hydrogen round-trip.
//!
//! For every molecule in `implicitH.sdf`, deleting then re-adding hydrogens must
//! leave the atom count unchanged (no hydrogens lost or gained).

mod common;

use common::ob_test_file;
use openbabel::Molecule;

#[test]
fn delete_then_add_hydrogens_preserves_atom_count() {
    let path = ob_test_file("implicitH.sdf");
    let molecules = Molecule::read_file_many(path.to_str().unwrap(), Some("sdf"))
        .expect("read implicitH.sdf");

    let mut checked = 0usize;
    for mut mol in molecules {
        let before = mol.num_atoms();
        mol.remove_hydrogens();
        mol.add_hydrogens();
        assert_eq!(
            before,
            mol.num_atoms(),
            "atom count changed for {:?} after delete+add hydrogens",
            mol.title()
        );
        checked += 1;
    }

    assert!(checked > 0, "no molecules were read from implicitH.sdf");
}
