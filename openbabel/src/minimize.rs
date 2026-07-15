//! Configurable, constraint-aware force-field geometry minimization with an
//! optional step-by-step trajectory.
//!
//! [`Minimizer`] bundles the choices OpenBabel exposes — force field,
//! [`Algorithm`], step budget, energy-convergence threshold, and a
//! [`Constraints`] set — into one config. Feed it to
//! [`Molecule::optimize_geometry_with`](crate::Molecule::optimize_geometry_with)
//! to minimize in one shot, or to [`Molecule::minimize`](crate::Molecule::minimize)
//! to get the whole trajectory back as an [`Optimization`] iterator of
//! [`OptStep`] frames (step count, energy, coordinates).
//!
//! ```no_run
//! use openbabel::{Algorithm, Minimizer, Molecule};
//!
//! let mut mol = Molecule::parse("CCO", "smi").unwrap();
//! mol.generate_3d();
//!
//! let mut cfg = Minimizer::new("MMFF94");
//! cfg.algorithm(Algorithm::ConjugateGradients)
//!    .max_steps(500)
//!    .energy_convergence(1e-6)
//!    .steps_per_frame(10);
//!
//! for step in mol.minimize(&cfg) {
//!     println!("step {:>4}  E = {:.4}", step.step, step.energy);
//! }
//! ```
//!
//! Because OpenBabel is not thread-safe and keeps force-field state in *static*
//! (shared) members, an optimization must run as one atomic unit under the
//! crate's global lock. [`minimize`](crate::Molecule::minimize) therefore runs
//! the whole minimization eagerly and hands back the *recorded* trajectory: the
//! molecule is left at the final geometry, and each [`OptStep`] carries the
//! geometry captured at that frame.
//!
//! Convergence note: OpenBabel 3.2.1 lets you set only the **energy** criterion
//! (`energy_convergence`). Its steepest-descent / conjugate-gradient / L-BFGS
//! minimizers combine it with a *fixed* internal gradient criterion; there is no
//! public setter for a gradient tolerance, so this API does not expose one.
//!
//! L-BFGS caveat: in the vendored OpenBabel 3.2.1 the L-BFGS minimizer corrupts
//! the heap when paired with the **UFF** force field (a bug in that
//! locally-added minimizer). That one pairing is refused — it returns `None` /
//! an empty trajectory rather than running. Use conjugate gradients or steepest
//! descent with UFF; L-BFGS is fine with MMFF94, MMFF94s, GAFF, and Ghemical.

use openbabel_sys::ffi;

use crate::{with_ob, Constraints, Molecule};

/// Minimization algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Algorithm {
    /// Steepest descent: robust, slow to converge.
    SteepestDescent,
    /// Conjugate gradients: a good general-purpose default.
    ConjugateGradients,
    /// Limited-memory BFGS: often the fastest to converge.
    Lbfgs,
}

impl Algorithm {
    fn code(self) -> u32 {
        match self {
            Algorithm::SteepestDescent => 0,
            Algorithm::ConjugateGradients => 1,
            Algorithm::Lbfgs => 2,
        }
    }
}

/// Configuration for a force-field minimization.
///
/// Builder methods take `&mut self` and return `&mut Self` for chaining.
pub struct Minimizer {
    forcefield: String,
    algorithm: Algorithm,
    max_steps: u32,
    energy_convergence: f64,
    steps_per_frame: u32,
    constraints: Constraints,
}

impl Minimizer {
    /// A new config for the named force field (`"MMFF94"`, `"MMFF94s"`,
    /// `"UFF"`, `"GAFF"`, `"Ghemical"`). Defaults: conjugate gradients, 500
    /// steps, `1e-6` energy convergence, one step per trajectory frame, no
    /// constraints.
    pub fn new(forcefield: &str) -> Self {
        Minimizer {
            forcefield: forcefield.to_string(),
            algorithm: Algorithm::ConjugateGradients,
            max_steps: 500,
            energy_convergence: 1e-6,
            steps_per_frame: 1,
            constraints: Constraints::new(),
        }
    }

    /// Choose the minimization [`Algorithm`].
    pub fn algorithm(&mut self, algorithm: Algorithm) -> &mut Self {
        self.algorithm = algorithm;
        self
    }

    /// Maximum number of optimization steps.
    pub fn max_steps(&mut self, steps: u32) -> &mut Self {
        self.max_steps = steps;
        self
    }

    /// Energy-convergence threshold: minimization stops once the step-to-step
    /// energy change drops below this (in the force field's energy unit).
    pub fn energy_convergence(&mut self, econv: f64) -> &mut Self {
        self.energy_convergence = econv;
        self
    }

    /// How many optimization steps each [`OptStep`] frame advances when
    /// iterating with [`Molecule::minimize`](crate::Molecule::minimize). Larger
    /// values mean coarser, cheaper trajectories. Clamped to at least 1. Has no
    /// effect on [`optimize_geometry_with`](crate::Molecule::optimize_geometry_with).
    pub fn steps_per_frame(&mut self, steps: u32) -> &mut Self {
        self.steps_per_frame = steps.max(1);
        self
    }

    /// Apply a [`Constraints`] set (fixed atoms, distance/angle/torsion
    /// restraints, ignored atoms).
    pub fn constraints(&mut self, constraints: Constraints) -> &mut Self {
        self.constraints = constraints;
        self
    }

    /// Whether this force-field / algorithm pairing is safe to run.
    ///
    /// The L-BFGS minimizer in the vendored OpenBabel 3.2.1 (a local addition,
    /// not part of upstream 3.2.1) corrupts the heap when driven with the UFF
    /// force field specifically. Since that is undefined behaviour we cannot
    /// contain, this combination is refused; every other pairing is fine (UFF
    /// works with steepest descent / conjugate gradients, and L-BFGS works with
    /// the MMFF94, MMFF94s, GAFF, and Ghemical force fields).
    fn is_supported(&self) -> bool {
        !(matches!(self.algorithm, Algorithm::Lbfgs) && self.forcefield.eq_ignore_ascii_case("uff"))
    }

    /// Run the minimization to completion (no trajectory), returning the final
    /// energy or `None` on an unsupported pairing / unknown force field / setup
    /// failure. Backs `Molecule::optimize_geometry_with`.
    pub(crate) fn run(&self, mol: &mut Molecule) -> Option<f64> {
        if !self.is_supported() {
            return None;
        }
        let mut ok = true;
        let e = with_ob(|| {
            ffi::optimizer_run_to_end(
                mol.as_inner_pin_mut(),
                &self.forcefield,
                self.algorithm.code(),
                self.max_steps,
                self.energy_convergence,
                self.constraints.as_inner(),
                &mut ok,
            )
        });
        if ok {
            Some(e)
        } else {
            None
        }
    }
}

/// One frame of an optimization [trajectory](Optimization): the geometry and
/// energy after another chunk of steps.
#[derive(Clone, Debug, PartialEq)]
pub struct OptStep {
    /// Cumulative optimization step count at this frame.
    pub step: u32,
    /// Potential energy at this frame, in the force field's energy unit.
    pub energy: f64,
    /// Atom coordinates at this frame, one `[x, y, z]` per atom in index order.
    pub coordinates: Vec<[f64; 3]>,
}

/// A recorded optimization trajectory: an [`Iterator`] (and
/// [`ExactSizeIterator`]) over the [`OptStep`] frames captured while minimizing.
///
/// Produced by [`Molecule::minimize`](crate::Molecule::minimize), which runs the
/// whole minimization eagerly (see the [module docs](self)) and leaves the
/// molecule at the final geometry. Iterate, `collect()`, or index the frames to
/// inspect the path.
pub struct Optimization {
    frames: std::vec::IntoIter<OptStep>,
}

impl Optimization {
    pub(crate) fn new(mol: &mut Molecule, config: &Minimizer) -> Self {
        if !config.is_supported() {
            return Optimization {
                frames: Vec::new().into_iter(),
            };
        }
        let flat = with_ob(|| {
            ffi::optimizer_run_trajectory(
                mol.as_inner_pin_mut(),
                &config.forcefield,
                config.algorithm.code(),
                config.max_steps,
                config.energy_convergence,
                config.constraints.as_inner(),
                config.steps_per_frame,
            )
        });
        // Each frame is [energy, x0, y0, z0, ...] — 1 + 3·num_atoms doubles.
        let frame_len = 1 + 3 * mol.num_atoms() as usize;
        let mut frames = Vec::new();
        let mut cumulative = 0u32;
        for chunk in flat.chunks_exact(frame_len) {
            cumulative += config.steps_per_frame;
            let coordinates = chunk[1..].chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            frames.push(OptStep {
                step: cumulative,
                energy: chunk[0],
                coordinates,
            });
        }
        Optimization {
            frames: frames.into_iter(),
        }
    }
}

impl Iterator for Optimization {
    type Item = OptStep;

    fn next(&mut self) -> Option<OptStep> {
        self.frames.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.frames.size_hint()
    }
}

impl ExactSizeIterator for Optimization {
    fn len(&self) -> usize {
        self.frames.len()
    }
}

impl std::fmt::Debug for Optimization {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Optimization")
            .field("remaining_frames", &self.frames.len())
            .finish()
    }
}
