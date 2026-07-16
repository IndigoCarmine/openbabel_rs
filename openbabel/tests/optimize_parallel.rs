//! The Rust force-field core runs its minimization loop **without** holding the
//! OpenBabel lock, so many molecules can be optimized at the same time. These
//! tests exercise that: the same optimizations spread across threads must match
//! the sequential results exactly and must not deadlock or corrupt memory.
//!
//! `Molecule: Send` (an `unsafe impl`, justified by the global OpenBabel lock)
//! is what lets a molecule move onto a worker thread in the first place.

use openbabel::Molecule;
use std::thread;

/// A spread of small molecules, each exercising different MMFF94 terms.
const SMILES: &[&str] = &["CCO", "CCC", "CCN", "CC(=O)O", "c1ccccc1", "CCCCO", "CC(C)C", "OCCO"];

fn built(smiles: &str) -> Molecule {
    let mut m = Molecule::parse(smiles, "smi").expect("parse");
    assert!(m.generate_3d(), "3D generation failed for {smiles}");
    m
}

/// Optimizing on many threads yields exactly the sequential result. Each thread
/// works on its own clone of a shared starting geometry (built once, so the
/// comparison is deterministic despite `generate_3d`'s randomness).
#[test]
fn parallel_optimize_matches_sequential() {
    let starts: Vec<Molecule> = SMILES.iter().map(|s| built(s)).collect();

    // Sequential reference, each from a fresh clone of the starting geometry.
    let sequential: Vec<f64> = starts
        .iter()
        .map(|m| m.clone().optimize_geometry_rs("MMFF94", 500).expect("sequential optimize"))
        .collect();

    // The same work, one thread per molecule. The clone (an OpenBabel call)
    // happens on this thread; only the owned clone crosses to the worker, where
    // the lock-free minimizer runs concurrently with the others.
    let handles: Vec<_> = starts
        .iter()
        .map(|m| {
            let mut clone = m.clone();
            thread::spawn(move || clone.optimize_geometry_rs("MMFF94", 500).expect("parallel optimize"))
        })
        .collect();
    let parallel: Vec<f64> = handles.into_iter().map(|h| h.join().expect("thread panicked")).collect();

    for (i, (s, p)) in sequential.iter().zip(&parallel).enumerate() {
        assert!(
            (s - p).abs() < 1e-6,
            "{}: sequential {s} vs parallel {p} (lock-free optimize must be deterministic)",
            SMILES[i]
        );
    }
}

/// Many rounds hammering the shared OpenBabel lock from several threads must not
/// deadlock or corrupt state: every optimization still returns a finite energy.
#[test]
fn parallel_optimize_stress_is_stable() {
    let handles: Vec<_> = (0..16)
        .map(|i| {
            let smiles = SMILES[i % SMILES.len()].to_string();
            thread::spawn(move || {
                let mut m = built(&smiles);
                m.optimize_geometry_rs("MMFF94", 300).expect("optimize")
            })
        })
        .collect();

    for h in handles {
        let e = h.join().expect("thread panicked");
        assert!(e.is_finite(), "non-finite energy from a parallel optimization");
    }
}

/// The async wrapper drives the same lock-free optimization on Tokio's blocking
/// pool; awaiting many at once optimizes them concurrently.
#[cfg(feature = "async")]
#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn async_optimize_runs_concurrently() {
    let starts: Vec<Molecule> = SMILES.iter().map(|s| built(s)).collect();

    let tasks: Vec<_> = starts
        .into_iter()
        .map(|m| tokio::spawn(m.optimize_geometry_rs_async("MMFF94", 500)))
        .collect();

    for task in tasks {
        let (mol, energy) = task.await.expect("join");
        assert!(energy.expect("optimize").is_finite());
        assert!(mol.num_atoms() > 0);
    }
}
