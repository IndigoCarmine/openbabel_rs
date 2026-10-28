//! Force-field results must not depend on what was evaluated before.
//!
//! OpenBabel keeps one shared plugin instance per force field and caches its
//! last setup. Before the fix, `energy` / `optimize_geometry*` / `minimize` on
//! the same geometry returned different numbers depending on which molecule had
//! been set up, optimized, or built in 3D just before (`generate_3d` leaves the
//! shared MMFF94 instance with distance cut-offs switched on and a stale
//! non-bonded pair list; same-shaped molecules reused each other's terms;
//! constraints leaked between calls). These tests pin every evaluation to the
//! value a fresh, history-free evaluation gives.

use openbabel::{Algorithm, Constraints, Minimizer, Molecule};

const FORCEFIELDS: [&str; 5] = ["MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical"];

/// Bit-for-bit is what we expect; 1e-9 leaves room only for nothing at all.
const TOL: f64 = 1e-9;

struct Fixture {
    name: &'static str,
    mol: Molecule,
}

fn build(name: &'static str, smiles: &str) -> Fixture {
    let mut mol = Molecule::parse(smiles, "smi").expect("parse");
    mol.add_hydrogens();
    assert!(mol.generate_3d(), "3D generation failed for {smiles}");
    Fixture { name, mol }
}

fn fixtures() -> Vec<Fixture> {
    vec![
        build("butane", "CCCC"),
        build("ethanol", "CCO"),
        build("benzene", "c1ccccc1"),
        build("methylammonium", "C[NH3+]"),
        build("acetate", "CC(=O)[O-]"),
    ]
}

fn coords(m: &Molecule) -> Vec<f64> {
    m.atoms()
        .flat_map(|a| {
            let (x, y, z) = a.coords();
            [x, y, z]
        })
        .collect()
}

fn assert_close(what: &str, got: Option<f64>, want: Option<f64>) {
    match (got, want) {
        (Some(g), Some(w)) => assert!(
            (g - w).abs() <= TOL,
            "{what}: {g} != {w} (diff {:e})",
            (g - w).abs()
        ),
        (None, None) => {}
        _ => panic!("{what}: {got:?} vs {want:?}"),
    }
}

/// A clone of `mol` evaluated with nothing else in between — but after a
/// *different* molecule, so no cached setup can be reused.
fn fresh_energy(mol: &Molecule, ff: &str) -> Option<f64> {
    let mut other = Molecule::parse("O", "smi").unwrap();
    other.add_hydrogens();
    let _ = other.energy(ff);
    mol.clone().energy(ff)
}

#[test]
fn energies_do_not_depend_on_history() {
    let fx = fixtures();
    // Reference values, each taken right after an unrelated molecule.
    let reference: Vec<Vec<Option<f64>>> = fx
        .iter()
        .map(|f| {
            FORCEFIELDS
                .iter()
                .map(|ff| fresh_energy(&f.mol, ff))
                .collect()
        })
        .collect();
    // Every force field handles every fixture, the ions included.
    for (i, f) in fx.iter().enumerate() {
        for (j, ff) in FORCEFIELDS.iter().enumerate() {
            assert!(reference[i][j].is_some(), "{}/{ff}: no energy", f.name);
        }
    }

    // Interleave in several orders, with 3D builds and optimizations of other
    // molecules (which reconfigure OpenBabel's shared instances) in between.
    let n = fx.len();
    for round in 0..3 {
        for step in 0..n * FORCEFIELDS.len() {
            let i = (step * (round + 2) + round) % n;
            let j = (step + round * 3) % FORCEFIELDS.len();
            let ff = FORCEFIELDS[j];
            if step % 4 == 0 {
                let mut junk = build("junk", ["CCCCO", "c1ccncc1", "CC(C)C"][step % 3]).mol;
                let _ = junk.optimize_geometry(FORCEFIELDS[(j + 1) % FORCEFIELDS.len()], 25);
            }
            let got = fx[i].mol.energy(ff);
            assert_close(
                &format!("round {round}: {}/{ff}", fx[i].name),
                got,
                reference[i][j],
            );
        }
    }
}

#[test]
fn energy_after_generate_3d_matches_fresh() {
    // `generate_3d` runs OpenBabel's gen3d op, which enables cut-offs on the
    // shared MMFF94 instance; it must not leak into later energies.
    let butane = build("butane", "CCCC").mol;
    for ff in FORCEFIELDS {
        let want = fresh_energy(&butane, ff);
        let _ = build("other", "CCCCCCCCCCCC").mol; // a long chain: many far pairs
        assert_close(&format!("butane/{ff} after gen3d"), butane.energy(ff), want);
    }
}

#[test]
fn optimization_does_not_depend_on_history() {
    let butane = build("butane", "CCCC").mol;
    for ff in FORCEFIELDS {
        // Fresh run.
        let mut a = butane.clone();
        let _ = fresh_energy(&butane, ff);
        let ea = a.optimize_geometry(ff, 50);
        // Same start, after unrelated work on the shared instance.
        let mut junk = build("junk", "CCCCCCO").mol;
        let _ = junk.optimize_geometry(ff, 30);
        let mut b = butane.clone();
        let eb = b.optimize_geometry(ff, 50);
        assert_close(&format!("optimize_geometry/{ff}"), eb, ea);
        for (x, y) in coords(&a).iter().zip(coords(&b)) {
            assert!(
                (x - y).abs() <= TOL,
                "optimize_geometry/{ff}: coordinates differ"
            );
        }

        let mut cfg = Minimizer::new(ff);
        cfg.algorithm(Algorithm::ConjugateGradients)
            .max_steps(40)
            .steps_per_frame(10);
        let mut c = butane.clone();
        let _ = fresh_energy(&butane, ff);
        let ec = c.optimize_geometry_with(&cfg);
        let _ = build("junk", "CC(=O)O").mol.energy(ff);
        let mut d = butane.clone();
        let ed = d.optimize_geometry_with(&cfg);
        assert_close(&format!("optimize_geometry_with/{ff}"), ed, ec);

        let mut e = butane.clone();
        let fe: Vec<f64> = e.minimize(&cfg).map(|s| s.energy).collect();
        let _ = build("junk", "c1ccccc1O").mol.energy(ff);
        let mut f = butane.clone();
        let ff_frames: Vec<f64> = f.minimize(&cfg).map(|s| s.energy).collect();
        assert_eq!(fe.len(), ff_frames.len(), "minimize/{ff}: frame count");
        for (x, y) in fe.iter().zip(&ff_frames) {
            assert!((x - y).abs() <= TOL, "minimize/{ff}: {x} vs {y}");
        }
    }
}

#[test]
fn constraints_do_not_leak_into_later_calls() {
    let butane = build("butane", "CCCC").mol;
    for ff in FORCEFIELDS {
        let want_energy = fresh_energy(&butane, ff);
        let mut want_opt = butane.clone();
        let want_e = want_opt.optimize_geometry(ff, 40);

        // A constrained run on a same-topology molecule: fixed atom, distance
        // and torsion restraints, and an ignored atom (which changes the term
        // lists themselves).
        let mut cons = Constraints::new();
        cons.fix_atom(0)
            .distance(0, 3, 3.2)
            .torsion(0, 1, 2, 3, 60.0)
            .ignore(5);
        let mut cfg = Minimizer::new(ff);
        cfg.algorithm(Algorithm::SteepestDescent)
            .max_steps(30)
            .constraints(cons);
        let mut constrained = butane.clone();
        let _ = constrained.optimize_geometry_with(&cfg);
        assert_close(
            &format!("{ff}: energy after constrained run"),
            butane.energy(ff),
            want_energy,
        );

        let _ = constrained.clone().minimize(&cfg).count();
        let mut again = butane.clone();
        let got_e = again.optimize_geometry(ff, 40);
        assert_close(
            &format!("{ff}: optimize after constrained run"),
            got_e,
            want_e,
        );
        for (x, y) in coords(&want_opt).iter().zip(coords(&again)) {
            assert!(
                (x - y).abs() <= TOL,
                "{ff}: optimized coordinates differ after constrained run"
            );
        }
    }
}

#[test]
fn scan_loop_matches_fresh_evaluations() {
    // set_coordinates + energy on one molecule, repeatedly — the fast path —
    // must agree with evaluating each geometry from scratch.
    let mut butane = build("butane", "CCCC").mol;
    for ff in FORCEFIELDS {
        let mut scan = Vec::new();
        let mut geoms = Vec::new();
        for k in 0..12 {
            butane.set_torsion(0, 1, 2, 3, (k as f64) * 30.0_f64.to_radians());
            let coords = coords(&butane);
            assert!(butane.set_coordinates(&coords));
            scan.push(butane.energy(ff));
            geoms.push(coords);
        }
        for (k, coords) in geoms.iter().enumerate() {
            let mut m = butane.clone();
            assert!(m.set_coordinates(coords));
            assert_close(
                &format!("scan {ff} step {k}"),
                scan[k],
                fresh_energy(&m, ff),
            );
        }
        // And the torsion actually changes the energy (the scan is not trivial).
        let e: Vec<f64> = scan.iter().map(|e| e.unwrap()).collect();
        assert!(
            e.iter().any(|x| (x - e[0]).abs() > 1e-3),
            "{ff}: flat scan {e:?}"
        );
    }
}

#[test]
fn same_shape_different_connectivity_is_not_reused() {
    // Butane and isobutane look identical to OpenBabel's own IsSetupNeeded():
    // C4H10 with the carbons first, every carbon of degree 4, and bond i of the
    // same order between the same elements in both. Only the wiring differs, so
    // OpenBabel would happily reuse one's terms for the other.
    let a = build("butane", "CCCC").mol;
    let b = build("isobutane", "CC(C)C").mol;
    let shape = |m: &Molecule| {
        let atoms: Vec<_> = m.atoms().map(|x| (x.atomic_number(), x.degree())).collect();
        let bonds: Vec<_> = m
            .bonds()
            .map(|x| {
                let z = |i| m.atom(i).unwrap().atomic_number();
                (x.order(), z(x.begin_atom_index()), z(x.end_atom_index()))
            })
            .collect();
        (atoms, bonds)
    };
    assert_eq!(
        shape(&a),
        shape(&b),
        "precondition: same shape to IsSetupNeeded"
    );
    for ff in FORCEFIELDS {
        let wa = fresh_energy(&a, ff);
        let wb = fresh_energy(&b, ff);
        let _ = a.energy(ff);
        assert_close(&format!("{ff}: isobutane after butane"), b.energy(ff), wb);
        assert_close(&format!("{ff}: butane after isobutane"), a.energy(ff), wa);
    }
}

/// `coords` (flat x, y, z) turned rigidly by `angle` radians about the unit
/// `axis` through the origin.
fn turned(coords: &[f64], axis: [f64; 3], angle: f64) -> Vec<f64> {
    let (s, c) = angle.sin_cos();
    let [kx, ky, kz] = axis;
    coords
        .chunks(3)
        .flat_map(|p| {
            let (x, y, z) = (p[0], p[1], p[2]);
            let kdp = kx * x + ky * y + kz * z;
            let kxp = [ky * z - kz * y, kz * x - kx * z, kx * y - ky * x];
            [
                x * c + kxp[0] * s + kx * kdp * (1.0 - c),
                y * c + kxp[1] * s + ky * kdp * (1.0 - c),
                z * c + kxp[2] * s + kz * kdp * (1.0 - c),
            ]
        })
        .collect()
}

#[test]
fn exactly_anti_and_perpendicular_torsions_are_not_charged_as_syn() {
    // UFF's gradient-free torsion term takes the angle from acos and falls back
    // to 1e-3 rad (syn) when the cosine rounds past -1 -- at an exactly anti
    // torsion, depending on the coordinates' last bits -- or when the dot
    // product is tiny, at an exactly perpendicular one. The energy then jumped
    // by 2V (8.87 kJ/mol for butane's C-C-C-C) between geometries that differ
    // only by a rigid rotation, e.g. a scan's -180 and +180 points. Every force
    // field must give a smooth, rotation-invariant energy there.
    let butane = build("butane", "CCCC").mol;
    let axis = {
        let a: [f64; 3] = [0.3, -0.7, 0.648];
        let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        [a[0] / n, a[1] / n, a[2] / n]
    };
    const DELTA: f64 = 1e-4;
    for ff in FORCEFIELDS {
        for deg in [180.0_f64, -180.0, 90.0, -90.0] {
            let at = |offset: f64| {
                let mut m = butane.clone();
                m.set_torsion(0, 1, 2, 3, deg.to_radians() + offset);
                coords(&m)
            };
            let energy = |c: &[f64]| {
                let mut m = butane.clone();
                assert!(m.set_coordinates(c));
                m.energy(ff).unwrap()
            };
            let base = at(0.0);
            // The same geometry seen from 32 orientations: OpenBabel's cosine
            // rounds to either side of -1 for some of them.
            let turns: Vec<f64> = (0..32)
                .map(|k| energy(&turned(&base, axis, 0.2 * k as f64)))
                .collect();
            // Off the exact torsion OpenBabel's own angle is used; their mean
            // is the energy at it to O(DELTA^2).
            let smooth = 0.5 * (energy(&at(DELTA)) + energy(&at(-DELTA)));
            for (k, e) in turns.iter().enumerate() {
                assert!(
                    (e - smooth).abs() < 1e-5,
                    "{ff} at {deg} deg, turn {k}: {e} but {smooth} either side ({turns:?})"
                );
            }
        }
    }
}
