//! Port of OpenBabel's `test/logp_psa.cpp` — the logP and TPSA descriptors.
//!
//! Three inline molecules with hard-coded reference values (from JOELib2) and
//! the C++ test's tight tolerance of 2e-6. `OBDescriptor::Predict` for `"logP"`
//! and `"TPSA"` maps to [`Molecule::logp`] / [`Molecule::tpsa`].

use openbabel::Molecule;

const TOL: f64 = 2e-6;

fn prepared(smiles: &str) -> Molecule {
    let mut mol = Molecule::parse(smiles, "smi").expect("parse SMILES");
    mol.add_hydrogens();
    mol
}

#[test]
fn logp_and_tpsa_reference_values() {
    // Oc1ccccc1OC — guaiacol
    let m = prepared("Oc1ccccc1OC");
    assert!((m.logp().expect("logP") - 1.4008).abs() < TOL);
    assert!((m.tpsa().expect("TPSA") - 29.46).abs() < TOL);

    // c1ccccc1CBr — benzyl bromide
    let m = prepared("c1ccccc1CBr");
    assert!((m.logp().expect("logP") - 2.5815).abs() < TOL);
    assert!((m.tpsa().expect("TPSA") - 0.0).abs() < TOL);

    // Cc1ccccc1NC(=O)C — an acetanilide
    let m = prepared("Cc1ccccc1NC(=O)C");
    assert!((m.logp().expect("logP") - 2.0264).abs() < TOL);
    assert!((m.tpsa().expect("TPSA") - 29.1).abs() < TOL);
}
