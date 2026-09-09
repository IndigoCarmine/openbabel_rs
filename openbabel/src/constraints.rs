//! Force-field constraints for constrained geometry optimization.
//!
//! Build a [`Constraints`] set and hand it to a
//! [`Minimizer`](crate::Minimizer) to restrain a minimization: hold atoms in
//! place, pin bond lengths / angles / torsions to target values, or exclude
//! atoms from the calculation entirely.
//!
//! All atom indices are 0-based (matching [`Atom::index`](crate::Atom::index)).
//!
//! ```no_run
//! use openbabel::Constraints;
//! let mut c = Constraints::new();
//! c.fix_atom(0)                 // hold atom 0 in place
//!  .distance(1, 2, 1.54)        // pin the 1–2 bond to 1.54 Å
//!  .force_factor(50000.0);      // stiffen the restraints
//! ```

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::with_ob;

/// A Cartesian axis, for pinning a single coordinate of an atom in place with
/// [`Constraints::fix_atom_axis`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    /// The x axis.
    X,
    /// The y axis.
    Y,
    /// The z axis.
    Z,
}

/// A set of force-field restraints (`OBFFConstraints`).
///
/// The builder methods take `&mut self` and return `&mut Self`, so they can be
/// chained or applied conditionally.
pub struct Constraints {
    inner: UniquePtr<ffi::Constraints>,
    /// Every restraint added so far, so the set can be rebuilt when the force
    /// factor changes. See [`Constraints::force_factor`].
    added: Vec<Restraint>,
    factor: Option<f64>,
}

/// A restraint as it was requested, kept so it can be replayed onto a fresh
/// `OBFFConstraints` after a force-factor change.
#[derive(Clone, Copy)]
enum Restraint {
    Ignore(u32),
    Atom(u32),
    AtomAxis(u32, Axis),
    Distance(u32, u32, f64),
    Angle(u32, u32, u32, f64),
    Torsion(u32, u32, u32, u32, f64),
}

impl Constraints {
    /// Create an empty constraint set (no restraints).
    pub fn new() -> Self {
        Constraints {
            inner: with_ob(ffi::constraints_new),
            added: Vec::new(),
            factor: None,
        }
    }

    /// Replay every restraint onto a fresh set, with the current factor applied
    /// first so it reaches all of them.
    fn rebuild(&mut self) {
        let mut inner = with_ob(ffi::constraints_new);
        if let Some(factor) = self.factor {
            with_ob(|| ffi::constraints_set_factor(inner.pin_mut(), factor));
        }
        for r in &self.added {
            with_ob(|| match *r {
                Restraint::Ignore(a) => ffi::constraints_add_ignore(inner.pin_mut(), a),
                Restraint::Atom(a) => ffi::constraints_add_atom(inner.pin_mut(), a),
                Restraint::AtomAxis(a, Axis::X) => {
                    ffi::constraints_add_atom_x(inner.pin_mut(), a)
                }
                Restraint::AtomAxis(a, Axis::Y) => {
                    ffi::constraints_add_atom_y(inner.pin_mut(), a)
                }
                Restraint::AtomAxis(a, Axis::Z) => {
                    ffi::constraints_add_atom_z(inner.pin_mut(), a)
                }
                Restraint::Distance(a, b, v) => {
                    ffi::constraints_add_distance(inner.pin_mut(), a, b, v)
                }
                Restraint::Angle(a, b, c, v) => {
                    ffi::constraints_add_angle(inner.pin_mut(), a, b, c, v)
                }
                Restraint::Torsion(a, b, c, d, v) => {
                    ffi::constraints_add_torsion(inner.pin_mut(), a, b, c, d, v)
                }
            });
        }
        self.inner = inner;
    }

    /// Exclude `atom` from the force-field calculation entirely (its
    /// interactions are ignored, as if it weren't there).
    pub fn ignore(&mut self, atom: u32) -> &mut Self {
        with_ob(|| ffi::constraints_add_ignore(self.inner.pin_mut(), atom));
        self.added.push(Restraint::Ignore(atom));
        self
    }

    /// Fix `atom` at its current position (all three coordinates).
    pub fn fix_atom(&mut self, atom: u32) -> &mut Self {
        with_ob(|| ffi::constraints_add_atom(self.inner.pin_mut(), atom));
        self.added.push(Restraint::Atom(atom));
        self
    }

    /// Fix a single coordinate ([`Axis`]) of `atom`, leaving the others free.
    pub fn fix_atom_axis(&mut self, atom: u32, axis: Axis) -> &mut Self {
        with_ob(|| match axis {
            Axis::X => ffi::constraints_add_atom_x(self.inner.pin_mut(), atom),
            Axis::Y => ffi::constraints_add_atom_y(self.inner.pin_mut(), atom),
            Axis::Z => ffi::constraints_add_atom_z(self.inner.pin_mut(), atom),
        });
        self.added.push(Restraint::AtomAxis(atom, axis));
        self
    }

    /// Restrain the distance between atoms `a` and `b` to `length` (Å).
    pub fn distance(&mut self, a: u32, b: u32, length: f64) -> &mut Self {
        with_ob(|| ffi::constraints_add_distance(self.inner.pin_mut(), a, b, length));
        self.added.push(Restraint::Distance(a, b, length));
        self
    }

    /// Restrain the `a`–`b`–`c` valence angle to `degrees`.
    pub fn angle(&mut self, a: u32, b: u32, c: u32, degrees: f64) -> &mut Self {
        with_ob(|| ffi::constraints_add_angle(self.inner.pin_mut(), a, b, c, degrees));
        self.added.push(Restraint::Angle(a, b, c, degrees));
        self
    }

    /// Restrain the `a`–`b`–`c`–`d` torsion angle to `degrees`.
    pub fn torsion(&mut self, a: u32, b: u32, c: u32, d: u32, degrees: f64) -> &mut Self {
        with_ob(|| ffi::constraints_add_torsion(self.inner.pin_mut(), a, b, c, d, degrees));
        self.added.push(Restraint::Torsion(a, b, c, d, degrees));
        self
    }

    /// Set the force constant applied to the restraints (larger = stiffer;
    /// OpenBabel's default of 50 000 is used until set).
    ///
    /// Order does not matter: OpenBabel copies its factor into each restraint
    /// as that restraint is added, so `SetFactor` on the underlying
    /// `OBFFConstraints` reaches only the restraints added *after* it — a set
    /// built the obvious way, `c.torsion(..).force_factor(1e8)`, silently keeps
    /// the default and the caller sees a restraint that does not hold. This
    /// rebuilds the set so the factor always applies to all of them.
    ///
    /// How stiff is stiff enough depends on the restraint. Distance and angle
    /// restraints are harmonic and unbounded, so the default holds them. A
    /// torsion restraint is a cosine well whose depth is only `0.002 * factor`,
    /// which the force field's own torsional term can overwhelm: pinning
    /// butane's C–C–C–C dihedral to 90° needs about 1e8 (at 1e6 it lands within
    /// 5°, at the default it does not move at all).
    pub fn force_factor(&mut self, factor: f64) -> &mut Self {
        self.factor = Some(factor);
        self.rebuild();
        self
    }

    pub(crate) fn as_inner(&self) -> &ffi::Constraints {
        self.inner.as_ref().expect("Constraints is never null")
    }
}

impl Default for Constraints {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for Constraints {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Constraints").finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Algorithm, Minimizer, Molecule};

    fn butane() -> Molecule {
        let mut mol = Molecule::parse("CCCC", "smi").expect("parse");
        mol.add_hydrogens();
        assert!(mol.generate_3d(), "gen3d");
        mol
    }

    fn pin_torsion(factor_first: bool) -> f64 {
        const TARGET: f64 = 90.0;
        const FACTOR: f64 = 1.0e8;

        let mut mol = butane();
        let mut c = Constraints::new();
        if factor_first {
            c.force_factor(FACTOR);
            c.torsion(0, 1, 2, 3, TARGET);
        } else {
            c.torsion(0, 1, 2, 3, TARGET);
            c.force_factor(FACTOR);
        }

        let mut cfg = Minimizer::new("UFF");
        cfg.algorithm(Algorithm::ConjugateGradients)
            .max_steps(5000)
            .constraints(c);
        let _: Vec<_> = mol.minimize(&cfg).collect();
        mol.torsion(0, 1, 2, 3).abs()
    }

    /// The bug this rebuild exists for.
    ///
    /// `OBFFConstraints` copies its factor into each restraint as the restraint
    /// is added, so raising the factor afterwards reaches nothing already there.
    /// Built the obvious way — add the restraint, then set the stiffness — a
    /// torsion restraint kept the default 50 000, which its bounded cosine well
    /// cannot hold a dihedral with: butane stayed at anti (~170°) no matter what
    /// factor was asked for, all the way to 1e10.
    ///
    /// Both orders must now reach the target.
    #[test]
    fn the_force_factor_applies_whatever_the_order() {
        crate::init();
        let factor_first = pin_torsion(true);
        let factor_last = pin_torsion(false);

        for (label, got) in [("factor first", factor_first), ("factor last", factor_last)] {
            assert!(
                (got - 90.0).abs() < 5.0,
                "{label}: the torsion should be held at 90°, reads {got:.1}°                  (an unrestrained butane relaxes to ~180°)"
            );
        }
    }

    /// Distance and angle restraints are harmonic and unbounded, so they hold at
    /// the default factor — the counterpart to the torsion case above, and the
    /// reason the ordering bug went unnoticed.
    #[test]
    fn distance_and_angle_hold_at_the_default_factor() {
        crate::init();
        let mut mol = butane();

        let mut c = Constraints::new();
        c.distance(0, 1, 2.2);
        c.angle(0, 1, 2, 95.0);

        let mut cfg = Minimizer::new("UFF");
        cfg.algorithm(Algorithm::ConjugateGradients)
            .max_steps(5000)
            .constraints(c);
        let _: Vec<_> = mol.minimize(&cfg).collect();

        let d = mol.distance(0, 1);
        let a = mol.angle(0, 1, 2);
        assert!((d - 2.2).abs() < 0.2, "distance restraint: {d:.3} Å, wanted 2.2");
        assert!((a - 95.0).abs() < 5.0, "angle restraint: {a:.1}°, wanted 95");
    }

    /// Replaying the restraints must not lose or duplicate any of them.
    #[test]
    fn a_rebuild_keeps_every_restraint() {
        crate::init();
        let mut mol = butane();

        let mut c = Constraints::new();
        c.distance(0, 1, 2.2);
        c.angle(0, 1, 2, 95.0);
        // Triggers the rebuild, after both restraints are in.
        c.force_factor(1.0e8);

        let mut cfg = Minimizer::new("UFF");
        cfg.algorithm(Algorithm::ConjugateGradients)
            .max_steps(5000)
            .constraints(c);
        let _: Vec<_> = mol.minimize(&cfg).collect();

        let d = mol.distance(0, 1);
        let a = mol.angle(0, 1, 2);
        assert!((d - 2.2).abs() < 0.2, "distance lost in the rebuild: {d:.3} Å");
        assert!((a - 95.0).abs() < 5.0, "angle lost in the rebuild: {a:.1}°");
    }
}
