//! Integration tests for InChI / InChIKey identifiers (T4).
//!
//! These require InChI support to be compiled into the linked OpenBabel
//! (`WITH_INCHI=ON` in `openbabel-sys/build.rs`).

use openbabel::Molecule;

#[test]
fn ethanol_inchi_and_key() {
    let mol = Molecule::parse("CCO", "smi").expect("parse ethanol");

    let inchi = mol.inchi().expect("InChI support should be compiled in");
    assert_eq!(inchi, "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");

    let key = mol.inchikey().expect("InChIKey");
    assert_eq!(key, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N");
}

#[test]
fn inchi_is_canonical_across_input_order() {
    // The same molecule written two ways yields the same InChI.
    let a = Molecule::parse("OCC", "smi").expect("parse");
    let b = Molecule::parse("CCO", "smi").expect("parse");
    assert_eq!(a.inchi(), b.inchi());
    assert_eq!(a.inchikey(), b.inchikey());
}

#[test]
fn inchikey_has_standard_shape() {
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let key = mol.inchikey().expect("InChIKey");
    // Standard InChIKey: 14 + '-' + 10 + '-' + 1 = 27 chars.
    assert_eq!(key.len(), 27, "unexpected InChIKey {key:?}");
    let parts: Vec<&str> = key.split('-').collect();
    assert_eq!(parts.len(), 3);
    assert_eq!(parts[0].len(), 14);
    assert_eq!(parts[1].len(), 10);
    assert_eq!(parts[2].len(), 1);
}
