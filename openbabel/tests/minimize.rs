//! Integration tests for configurable, constrained, and streamed force-field
//! minimization (T10).

use openbabel::{Algorithm, Axis, Constraints, Minimizer, Molecule};

fn ethanol_3d() -> Molecule {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse ethanol");
    assert!(mol.generate_3d(), "generate 3D");
    mol
}

fn distance(mol: &Molecule, a: u32, b: u32) -> f64 {
    let (ax, ay, az) = mol.atom(a).unwrap().coords();
    let (bx, by, bz) = mol.atom(b).unwrap().coords();
    ((ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2)).sqrt()
}

#[test]
fn each_algorithm_minimizes_without_raising_energy() {
    for algo in [
        Algorithm::SteepestDescent,
        Algorithm::ConjugateGradients,
        Algorithm::Lbfgs,
    ] {
        // Compare energies within a single minimization run (one atomic
        // operation): the last frame must not sit above the first.
        let mut mol = ethanol_3d();
        let mut cfg = Minimizer::new("MMFF94");
        cfg.algorithm(algo).max_steps(200).steps_per_frame(20);
        let frames: Vec<_> = mol.minimize(&cfg).collect();
        assert!(!frames.is_empty(), "{algo:?}: no frames");
        let first = frames.first().unwrap().energy;
        let last = frames.last().unwrap().energy;
        assert!(last.is_finite(), "{algo:?}: non-finite energy");
        assert!(last <= first + 1e-3, "{algo:?}: energy rose {first} -> {last}");
    }
}

#[test]
fn unknown_forcefield_returns_none() {
    let mut mol = ethanol_3d();
    let cfg = Minimizer::new("NoSuchFF");
    assert!(mol.optimize_geometry_with(&cfg).is_none());
    // The trajectory iterator yields nothing when setup fails.
    assert_eq!(mol.minimize(&cfg).count(), 0);
}

#[test]
fn trajectory_streams_frames() {
    let mut mol = ethanol_3d();
    let n_atoms = mol.num_atoms() as usize;

    let mut cfg = Minimizer::new("MMFF94");
    cfg.algorithm(Algorithm::SteepestDescent)
        .max_steps(60)
        .steps_per_frame(5);

    let frames: Vec<_> = mol.minimize(&cfg).collect();
    assert!(!frames.is_empty(), "should yield at least one frame");

    for f in &frames {
        assert_eq!(f.coordinates.len(), n_atoms, "frame carries full geometry");
        assert!(f.energy.is_finite());
    }
    // The cumulative step counter advances by steps_per_frame.
    assert_eq!(frames[0].step, 5);
    // Overall the energy does not increase across the trajectory.
    let first = frames.first().unwrap().energy;
    let last = frames.last().unwrap().energy;
    assert!(last <= first + 1e-3, "energy rose across trajectory {first} -> {last}");
}

#[test]
fn fixed_atom_stays_put() {
    let mut mol = ethanol_3d();
    let before = mol.atom(0).unwrap().coords();

    let mut c = Constraints::new();
    c.fix_atom(0);
    let mut cfg = Minimizer::new("MMFF94");
    cfg.algorithm(Algorithm::ConjugateGradients)
        .max_steps(300)
        .constraints(c);
    mol.optimize_geometry_with(&cfg).expect("optimize");

    let after = mol.atom(0).unwrap().coords();
    let moved = ((before.0 - after.0).powi(2)
        + (before.1 - after.1).powi(2)
        + (before.2 - after.2).powi(2))
    .sqrt();
    assert!(moved < 1e-3, "fixed atom moved by {moved} Å");
}

#[test]
fn distance_constraint_stretches_bond() {
    let mut base = Molecule::parse("CC", "smi").expect("parse ethane");
    assert!(base.generate_3d());
    assert_eq!(base.atom(0).unwrap().atomic_number(), 6);
    assert_eq!(base.atom(1).unwrap().atomic_number(), 6);

    // Unconstrained minimum C–C distance.
    let mut free = base.clone();
    let mut cfg_free = Minimizer::new("MMFF94");
    cfg_free.max_steps(500);
    free.optimize_geometry_with(&cfg_free).expect("free optimize");
    let d_free = distance(&free, 0, 1);

    // Restrain the C–C distance to a clearly stretched 1.9 Å.
    let mut con = base.clone();
    let mut c = Constraints::new();
    c.distance(0, 1, 1.9).force_factor(100_000.0);
    let mut cfg_con = Minimizer::new("MMFF94");
    cfg_con.max_steps(500).constraints(c);
    con.optimize_geometry_with(&cfg_con).expect("constrained optimize");
    let d_con = distance(&con, 0, 1);

    assert!(
        d_con > d_free + 0.1,
        "distance constraint did not stretch the bond: free={d_free:.3}, con={d_con:.3}"
    );
}

#[test]
fn lbfgs_uff_is_refused_not_crash() {
    // L-BFGS + UFF triggers a heap-corruption bug in the vendored OpenBabel
    // 3.2.1, so the API refuses that one pairing rather than invoking it.
    let mut mol = ethanol_3d();
    let mut bad = Minimizer::new("UFF");
    bad.algorithm(Algorithm::Lbfgs).max_steps(100);
    assert!(mol.optimize_geometry_with(&bad).is_none());
    assert_eq!(mol.minimize(&bad).count(), 0);

    // UFF with the other algorithms is fine.
    let mut ok = Minimizer::new("UFF");
    ok.algorithm(Algorithm::ConjugateGradients).max_steps(100);
    assert!(mol.optimize_geometry_with(&ok).is_some());
}

#[test]
fn axis_fix_and_ignore_run() {
    // Exercise single-axis fix + ignore together, with L-BFGS on a force field
    // where L-BFGS is safe (MMFF94).
    let mut mol = ethanol_3d();
    let last = mol.num_atoms() - 1;
    let mut c = Constraints::new();
    c.fix_atom_axis(0, Axis::Z).ignore(last);

    let mut cfg = Minimizer::new("MMFF94");
    cfg.algorithm(Algorithm::Lbfgs).max_steps(100).constraints(c);
    assert!(mol.optimize_geometry_with(&cfg).is_some());
}
