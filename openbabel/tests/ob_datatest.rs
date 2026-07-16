//! Port of OpenBabel's `test/datatest.cpp` — element data tables.
//!
//! The C++ test's one substantive check is that helium's standard atomic mass
//! is 4.0026 (within 2e-3). `OBElements::GetMass(2)` maps to
//! [`openbabel::elements::mass`].

use openbabel::elements;

#[test]
fn helium_mass() {
    let mass = elements::mass(2);
    assert!(
        (mass - 4.0026).abs() < 2e-3,
        "helium mass expected ~4.0026, got {mass}"
    );
}
