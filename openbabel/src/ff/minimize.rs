//! Force-field-agnostic energy minimizer over an [`EnergyModel`].
//!
//! Runs entirely in Rust on a coordinate buffer and a term list — no OpenBabel
//! calls — so it holds no global state and is safe to run on many molecules in
//! parallel. Implements steepest descent and Polak–Ribière conjugate gradients
//! with a backtracking (Armijo) line search; the L-BFGS selector currently maps
//! to conjugate gradients.
//!
//! The stopping criterion matches OpenBabel's public knob: stop once the
//! step-to-step energy change drops below `econv`. A trajectory frame is
//! recorded every `steps_per_frame` steps (plus the final one).

use super::EnergyModel;
use crate::Algorithm;

/// One recorded step of a minimization trajectory.
pub(crate) struct Frame {
    pub step: u32,
    pub energy: f64,
    pub coords: Vec<f64>,
}

/// Result of a minimization: the final energy/geometry and the recorded frames.
pub(crate) struct MinOutcome {
    pub energy: f64,
    pub coords: Vec<f64>,
    pub frames: Vec<Frame>,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn neg(a: &[f64]) -> Vec<f64> {
    a.iter().map(|x| -x).collect()
}

fn inf_norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0, |m, x| m.max(x.abs()))
}

/// Backtracking Armijo line search along `dir` from `x`. Returns the accepted
/// point and its energy, or `None` if no sufficient decrease is found.
fn line_search<M: EnergyModel + ?Sized>(
    model: &M,
    x: &[f64],
    dir: &[f64],
    e0: f64,
    g: &[f64],
) -> Option<(Vec<f64>, f64)> {
    let slope = dot(g, dir); // directional derivative; < 0 for a descent dir
    if slope >= 0.0 {
        return None;
    }
    const C1: f64 = 1.0e-4;
    // Cap the first trial so the largest atom moves at most ~0.1 Å; this keeps
    // stiff terms (large force constants) from overshooting on step one.
    let scale = inf_norm(dir).max(1.0e-12);
    let mut alpha = (0.1 / scale).min(1.0);
    for _ in 0..40 {
        let trial: Vec<f64> = x.iter().zip(dir).map(|(xi, di)| xi + alpha * di).collect();
        let e = model.energy(&trial);
        if e <= e0 + C1 * alpha * slope {
            return Some((trial, e));
        }
        alpha *= 0.5;
    }
    None
}

/// Polak–Ribière⁺ CG update coefficient.
fn polak_ribiere(g_prev: &[f64], g: &[f64]) -> f64 {
    let denom = dot(g_prev, g_prev);
    if denom <= 0.0 {
        return 0.0;
    }
    let num: f64 = g.iter().zip(g_prev).map(|(gi, gp)| gi * (gi - gp)).sum();
    (num / denom).max(0.0)
}

/// Minimize `initial` under `model`. See the module docs for semantics.
///
/// `?Sized` so a `&dyn EnergyModel` (force-field-dispatched at runtime) works as
/// well as a concrete model.
pub(crate) fn minimize<M: EnergyModel + ?Sized>(
    model: &M,
    initial: &[f64],
    algorithm: Algorithm,
    max_steps: u32,
    econv: f64,
    steps_per_frame: u32,
) -> MinOutcome {
    let n = initial.len();
    let mut x = initial.to_vec();
    let mut g = vec![0.0; n];
    let mut energy = model.energy_gradient(&x, &mut g);
    let mut dir = neg(&g);
    let spf = steps_per_frame.max(1);
    let use_cg = matches!(algorithm, Algorithm::ConjugateGradients | Algorithm::Lbfgs);
    let mut frames = Vec::new();

    for step in 1..=max_steps {
        if dot(&dir, &g) >= 0.0 {
            dir = neg(&g); // fall back to steepest descent if not a descent dir
        }
        let (new_x, new_e) = match line_search(model, &x, &dir, energy, &g) {
            Some(v) => v,
            None => {
                // No further decrease: record where we are and stop.
                frames.push(Frame { step, energy, coords: x.clone() });
                break;
            }
        };
        let mut new_g = vec![0.0; n];
        model.energy_gradient(&new_x, &mut new_g);

        if use_cg {
            let beta = polak_ribiere(&g, &new_g);
            for i in 0..n {
                dir[i] = -new_g[i] + beta * dir[i];
            }
        } else {
            dir = neg(&new_g);
        }

        let converged = (energy - new_e).abs() < econv;
        x = new_x;
        g = new_g;
        energy = new_e;

        if step % spf == 0 || converged || step == max_steps {
            frames.push(Frame { step, energy, coords: x.clone() });
        }
        if converged {
            break;
        }
    }

    MinOutcome { energy, coords: x, frames }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ff::{BondTerm, Terms};
    use crate::ff::uff::UffModel;

    #[test]
    fn harmonic_bond_relaxes_to_equilibrium() {
        let terms = Terms {
            n_atoms: 2,
            bonds: vec![BondTerm { a: 0, b: 1, kb: 100.0, r0: 1.5 }],
            ..Default::default()
        };
        let m = UffModel::new(terms);
        // Start stretched at 2.2 Å.
        let start = [0.0, 0.0, 0.0, 2.2, 0.0, 0.0];
        let out = minimize(&m, &start, Algorithm::ConjugateGradients, 500, 1e-8, 1);
        // Energy should fall essentially to zero and the bond reach r0.
        assert!(out.energy < 1e-4, "final energy {}", out.energy);
        let dx = out.coords[3] - out.coords[0];
        let dy = out.coords[4] - out.coords[1];
        let dz = out.coords[5] - out.coords[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((r - 1.5).abs() < 1e-3, "final bond length {r}");
        assert!(!out.frames.is_empty());
    }

    #[test]
    fn steepest_descent_also_converges() {
        let terms = Terms {
            n_atoms: 2,
            bonds: vec![BondTerm { a: 0, b: 1, kb: 50.0, r0: 1.2 }],
            ..Default::default()
        };
        let m = UffModel::new(terms);
        let start = [0.0, 0.0, 0.0, 0.7, 0.0, 0.0]; // compressed
        let out = minimize(&m, &start, Algorithm::SteepestDescent, 1000, 1e-9, 10);
        assert!(out.energy < 1e-3, "final energy {}", out.energy);
    }
}
