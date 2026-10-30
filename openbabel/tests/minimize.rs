//! Integration tests for configurable, constrained, and streamed force-field
//! minimization (T10).

use openbabel::{Algorithm, Axis, Constraints, Minimizer, Molecule, StopReason};

fn ethanol_3d() -> Molecule {
    let mut mol = Molecule::parse("CCO", "smi").expect("parse ethanol");
    assert!(mol.generate_3d(), "generate 3D");
    mol
}

/// `[NH4+].CC(=O)[O-]` after add_hydrogens + generate_3d: an ion pair whose
/// UFF minimum lies near -0.64, which the minimizers used to call converged at
/// about +1.9 after a few hundred conjugate-gradient steps.
fn ammonium_acetate() -> Molecule {
    Molecule::parse(include_str!("data/ammonium_acetate.sdf"), "sdf").expect("fixture")
}

/// Taxol after generate_3d_with("med"), where GAFF L-BFGS used to stop after one
/// step (its gradient was NaN).
fn taxol() -> Molecule {
    Molecule::parse(include_str!("data/taxol_med.sdf"), "sdf").expect("fixture")
}

/// The largest per-atom gradient, |dE/dx_atom|, at the molecule's geometry.
fn max_atom_gradient(mol: &Molecule, ff: &str) -> f64 {
    let (_, grad) = mol.energy_gradient(ff, &Constraints::new()).expect("gradient");
    grad.iter()
        .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt())
        .fold(0.0, f64::max)
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
fn lbfgs_runs_with_uff() {
    // L-BFGS + UFF used to be refused: the pairing corrupted the heap. With
    // that fixed it is an ordinary pairing, and minimizes like the others.
    for (name, mut mol) in [("ethanol", ethanol_3d()), ("ammonium acetate", ammonium_acetate())] {
        let start = mol.energy("UFF").unwrap();
        let mut cfg = Minimizer::new("UFF");
        cfg.algorithm(Algorithm::Lbfgs).max_steps(2500);
        let run = mol.minimize(&cfg);
        assert_eq!(run.stop_reason(), StopReason::Converged, "{name}");
        let end = mol.energy("UFF").unwrap();
        assert!(end < start, "{name}: {start} -> {end}");
        assert!(mol.optimize_geometry_with(&cfg).is_some(), "{name}");
    }
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

#[test]
fn stop_reason_distinguishes_budget_from_convergence() {
    // A budget of one step cannot minimize anything: the run must report that it
    // was cut short rather than passing an unfinished geometry off as minimized.
    let mut mol = ethanol_3d();
    let mut tight = Minimizer::new("MMFF94");
    tight.max_steps(1).steps_per_frame(1).energy_convergence(1e-30);
    let run = mol.minimize(&tight);
    assert_eq!(run.stop_reason(), StopReason::MaxSteps);
    assert_eq!(run.steps_taken(), 1);

    // Given room and a reachable threshold, the same molecule converges — and
    // says so, having stopped before burning the whole budget.
    let mut loose = Minimizer::new("MMFF94");
    loose.max_steps(10_000).steps_per_frame(10).energy_convergence(1e-4);
    let run = mol.minimize(&loose);
    assert_eq!(run.stop_reason(), StopReason::Converged);
    assert!(run.steps_taken() < 10_000, "burned the budget: {}", run.steps_taken());
}

#[test]
fn steps_taken_never_exceeds_the_budget() {
    // steps_per_frame does not divide max_steps, so the final chunk is clamped
    // to what is left (2 + 2 + 1). The reported step counts must respect that
    // rather than overshoot, the way the old synthesized `step` counter did —
    // it just added steps_per_frame per frame and sailed past the budget.
    //
    // Both knobs matter: at the default threshold ethanol converges within a
    // handful of steps, so the budget would never be reached. A tight econv
    // pushes convergence out to ~16 steps, and a budget of 5 lands well short.
    let mut mol = ethanol_3d();
    let mut cfg = Minimizer::new("MMFF94");
    cfg.max_steps(5).steps_per_frame(2).energy_convergence(1e-30);
    let run = mol.minimize(&cfg);
    assert_eq!(run.stop_reason(), StopReason::MaxSteps);
    assert_eq!(run.steps_taken(), 5);
    let frames: Vec<_> = run.collect();
    let steps: Vec<u32> = frames.iter().map(|f| f.step).collect();
    assert_eq!(steps, vec![2, 4, 5], "final chunk was not clamped to the budget");
}

#[test]
fn unknown_forcefield_reports_failed() {
    let mut mol = ethanol_3d();
    let mut cfg = Minimizer::new("NoSuchFF");
    cfg.algorithm(Algorithm::Lbfgs);
    assert_eq!(mol.minimize(&cfg).stop_reason(), StopReason::Failed);
}

/// OpenBabel's convergence test is "energy change below econv and per-atom
/// gradient below 0.1". It used to take the *smallest* atom's gradient, so one
/// atom at rest -- a hydrogen, typically -- was enough: ammonium acetate
/// "converged" under UFF at E ~ +1.9 with ions still pulling on each other,
/// against a minimum near -0.64. A converged run must leave every atom at rest.
#[test]
fn converged_means_every_atom_is_at_rest() {
    for algo in [Algorithm::SteepestDescent, Algorithm::ConjugateGradients, Algorithm::Lbfgs] {
        let mut mol = ammonium_acetate();
        let mut cfg = Minimizer::new("UFF");
        cfg.algorithm(algo).max_steps(10_000);
        let run = mol.minimize(&cfg);
        let reason = run.stop_reason();
        let energy = mol.energy("UFF").unwrap();
        if reason == StopReason::Converged {
            let worst = max_atom_gradient(&mol, "UFF");
            assert!(
                worst < GRADIENT_AT_REST,
                "{algo:?}: converged with an atom still under a gradient of {worst:.3} (E = {energy:.3})"
            );
        }
        // Conjugate gradients from this geometry runs on into the minimum
        // (L-BFGS settles in a shallower one, at rest too).
        if algo == Algorithm::ConjugateGradients {
            assert_eq!(reason, StopReason::Converged, "{algo:?}: E = {energy:.3}");
            assert!(energy < -0.6, "{algo:?}: stopped at E = {energy:.3}, above the minimum near -0.64");
        }
    }
}

/// OpenBabel's per-atom gradient threshold (0.1), with room for the one line
/// search step taken after it was measured.
const GRADIENT_AT_REST: f64 = 0.12;

/// ConjugateGradientsInitialize took a line-search step of its own, which no
/// step count saw: a budget of one ran two. The first conjugate-gradient step
/// is a steepest-descent one (there is no previous direction), so a one-step
/// budget must land exactly where one steepest-descent step does.
#[test]
fn conjugate_gradients_counts_its_first_step() {
    for ff in ["MMFF94", "UFF"] {
        let mut sd = ethanol_3d();
        let mut cg = sd.clone();
        let mut cfg = Minimizer::new(ff);
        cfg.max_steps(1).energy_convergence(1e-30);
        cfg.algorithm(Algorithm::SteepestDescent);
        let sd_run = sd.minimize(&cfg);
        cfg.algorithm(Algorithm::ConjugateGradients);
        let cg_run = cg.minimize(&cfg);
        assert_eq!(cg_run.steps_taken(), 1);
        assert_eq!(cg_run.stop_reason(), sd_run.stop_reason());
        for i in 0..sd.num_atoms() {
            let (a, b) = (sd.atom(i).unwrap().coords(), cg.atom(i).unwrap().coords());
            let moved = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) + (a.2 - b.2).powi(2)).sqrt();
            assert!(moved < 1e-9, "{ff}: atom {i} is {moved:.2e} Å from the steepest-descent step");
        }
    }
}

/// GAFF's gradient on generated taxol was NaN on some atoms. L-BFGS (which,
/// unlike the other minimizers' line searches, does not zero a NaN direction)
/// then found no lower energy along it and stopped after one step, reporting
/// convergence at E ~ 3716 where conjugate gradients went on to ~1400.
#[test]
fn gaff_lbfgs_minimizes_taxol() {
    let mut mol = taxol();
    let start = mol.energy("GAFF").unwrap();
    let mut cfg = Minimizer::new("GAFF");
    cfg.algorithm(Algorithm::Lbfgs).max_steps(2500);
    let run = mol.minimize(&cfg);
    let end = mol.energy("GAFF").unwrap();
    assert!(run.steps_taken() > 10, "stopped after {} steps at E = {end:.1}", run.steps_taken());
    assert!(end < start - 50.0, "E {start:.1} -> {end:.1}");
}

/// Conjugate gradients mixed beta atom by atom (each atom's |g|^2 over its own
/// previous direction's), which builds neither conjugate nor necessarily
/// downhill directions: on taxol it was still descending after 5000 steps
/// under every force field. One Polak-Ribiere+ beta over all coordinates, with
/// a steepest-descent restart whenever the direction is not downhill,
/// converges in under 2500.
#[test]
fn conjugate_gradients_converges_on_taxol() {
    for ff in ["UFF", "MMFF94", "GAFF"] {
        let mut mol = taxol();
        let mut cfg = Minimizer::new(ff);
        cfg.algorithm(Algorithm::ConjugateGradients).max_steps(4000);
        let run = mol.minimize(&cfg);
        assert_eq!(
            run.stop_reason(),
            StopReason::Converged,
            "{ff}: still descending after {} steps, E = {:.1}",
            run.steps_taken(),
            mol.energy(ff).unwrap()
        );
    }
}

/// An angle restraint's energy is 0.0002 f (theta - target)^2 in degrees, but
/// its gradient was taken per radian -- a pull 57x weaker than the energy, so
/// minimizers stalled well short of the target. Under every force field the
/// default stiffness must now hold the angle at its target.
#[test]
fn angle_restraint_reaches_its_target() {
    let mut base = Molecule::parse("CCCC", "smi").unwrap();
    base.add_hydrogens();
    assert!(base.generate_3d());
    for ff in ["UFF", "GAFF", "MMFF94"] {
        let mut mol = base.clone();
        let mut c = Constraints::new();
        c.angle(0, 1, 2, 95.0);
        let mut cfg = Minimizer::new(ff);
        cfg.max_steps(2500).constraints(c);
        let run = mol.minimize(&cfg);
        let angle = mol.angle(0, 1, 2);
        assert_eq!(run.stop_reason(), StopReason::Converged, "{ff}: angle {angle:.3}");
        // The restraint is stiff but finite, so it settles a fraction of a
        // degree towards the force field's own angle (UFF: ~0.4 deg). Before,
        // the runs stalled ~1.2 deg out and never converged.
        assert!((angle - 95.0).abs() < 0.5, "{ff}: angle {angle:.3}, wanted 95");
    }
}
