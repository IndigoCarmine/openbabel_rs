//! Port of OpenBabel's `test/spectrophoretest.cpp` — the Spectrophore™ descriptor.
//!
//! Only `test01` is portable: it uses `OBSpectrophore`'s default settings
//! (NoNormalization, resolution 3.0, AngStepSize20, NoStereoSpecificProbes),
//! which is exactly what the binding's [`Molecule::spectrophore`] computes.
//! test02–test09 configure resolution / accuracy / normalization / stereo via
//! setters the binding does not expose, so they are not ported. The molecule
//! (C, H, F, Br, Cl around a central carbon) is built atom-by-atom with explicit
//! coordinates, exactly as the C++ does.
//!
//! # Marked `#[ignore]` — the vendored expected bounds are stale
//!
//! The first property block (`r[0..12]`, the partial-charge probes) matches the
//! C++ bounds exactly, but later blocks do not: e.g. `r[12]` is `2.2932`, while
//! `test01` expects `(2.707, 2.709)`. This is an upstream data staleness, not a
//! binding bug — running OpenBabel 3.2.1's own `obspectrophore` tool on this
//! exact molecule yields `2.2932` at index 12, matching the binding. OpenBabel's
//! own `spectrophore` test would therefore fail against these bounds with this
//! build. The faithful port is kept (runs under `cargo test -- --ignored`) to
//! document the discrepancy.

mod common;

use openbabel::Molecule;

#[test]
#[ignore = "vendored spectrophoretest bounds do not match OpenBabel 3.2.1's own \
            obspectrophore output (r[12]=2.2932, not 2.707..2.709); upstream reference is stale"]
fn spectrophore_default_settings() {
    let mut mol = Molecule::new();
    let a0 = mol.add_atom(6);
    mol.atom_mut(a0).unwrap().set_position(-0.013, 1.086, 0.008);
    let a1 = mol.add_atom(1);
    mol.atom_mut(a1).unwrap().set_position(0.002, -0.004, 0.002);
    let a2 = mol.add_atom(9);
    mol.atom_mut(a2).unwrap().set_position(1.300, 1.570, -0.002);
    let a3 = mol.add_atom(35);
    mol.atom_mut(a3).unwrap().set_position(-0.964, 1.737, -1.585);
    let a4 = mol.add_atom(17);
    mol.atom_mut(a4).unwrap().set_position(-0.857, 1.667, 1.491);
    for i in [a1, a2, a3, a4] {
        mol.add_bond(a0, i, 1);
    }

    let r = mol.spectrophore();
    assert_eq!(r.len(), 48, "default Spectrophore has 48 values");

    // (low, high) bounds for r[0..48], verbatim from spectrophoretest.cpp test01.
    let bounds: [(f64, f64); 48] = [
        (1.598, 1.600), (1.576, 1.578), (1.169, 1.171), (3.760, 3.762),
        (5.174, 5.176), (5.780, 5.782), (3.796, 3.798), (3.712, 3.714),
        (4.650, 4.652), (7.736, 7.738), (7.949, 7.951), (4.868, 4.870),
        (2.707, 2.709), (3.470, 3.472), (6.698, 6.700), (9.485, 9.487),
        (7.667, 7.669), (8.881, 8.883), (4.899, 4.901), (7.478, 7.480),
        (9.323, 9.325), (10.292, 10.294), (12.955, 12.957), (10.334, 10.336),
        (4.020, 4.022), (3.813, 3.815), (2.946, 2.948), (6.380, 6.382),
        (11.003, 11.005), (8.278, 8.280), (6.548, 6.550), (7.135, 7.137),
        (8.612, 8.614), (13.181, 13.183), (13.743, 13.745), (9.083, 9.085),
        (0.458, 0.460), (0.641, 0.643), (2.171, 2.173), (2.752, 2.754),
        (2.347, 2.349), (2.604, 2.606), (1.613, 1.615), (3.165, 3.167),
        (3.390, 3.392), (3.131, 3.133), (4.104, 4.106), (2.874, 2.876),
    ];
    for (i, (low, high)) in bounds.iter().enumerate() {
        assert!(
            r[i] > *low && r[i] < *high,
            "spectrophore[{i}] = {} not in ({low}, {high})",
            r[i]
        );
    }
}
