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
}

impl Constraints {
    /// Create an empty constraint set (no restraints).
    pub fn new() -> Self {
        Constraints {
            inner: with_ob(ffi::constraints_new),
        }
    }

    /// Exclude `atom` from the force-field calculation entirely (its
    /// interactions are ignored, as if it weren't there).
    pub fn ignore(&mut self, atom: u32) -> &mut Self {
        with_ob(|| ffi::constraints_add_ignore(self.inner.pin_mut(), atom));
        self
    }

    /// Fix `atom` at its current position (all three coordinates).
    pub fn fix_atom(&mut self, atom: u32) -> &mut Self {
        with_ob(|| ffi::constraints_add_atom(self.inner.pin_mut(), atom));
        self
    }

    /// Fix a single coordinate ([`Axis`]) of `atom`, leaving the others free.
    pub fn fix_atom_axis(&mut self, atom: u32, axis: Axis) -> &mut Self {
        with_ob(|| match axis {
            Axis::X => ffi::constraints_add_atom_x(self.inner.pin_mut(), atom),
            Axis::Y => ffi::constraints_add_atom_y(self.inner.pin_mut(), atom),
            Axis::Z => ffi::constraints_add_atom_z(self.inner.pin_mut(), atom),
        });
        self
    }

    /// Restrain the distance between atoms `a` and `b` to `length` (Å).
    pub fn distance(&mut self, a: u32, b: u32, length: f64) -> &mut Self {
        with_ob(|| ffi::constraints_add_distance(self.inner.pin_mut(), a, b, length));
        self
    }

    /// Restrain the `a`–`b`–`c` valence angle to `degrees`.
    pub fn angle(&mut self, a: u32, b: u32, c: u32, degrees: f64) -> &mut Self {
        with_ob(|| ffi::constraints_add_angle(self.inner.pin_mut(), a, b, c, degrees));
        self
    }

    /// Restrain the `a`–`b`–`c`–`d` torsion angle to `degrees`.
    pub fn torsion(&mut self, a: u32, b: u32, c: u32, d: u32, degrees: f64) -> &mut Self {
        with_ob(|| ffi::constraints_add_torsion(self.inner.pin_mut(), a, b, c, d, degrees));
        self
    }

    /// Set the harmonic force constant applied to the restraints (larger =
    /// stiffer; OpenBabel's default is used until set).
    pub fn force_factor(&mut self, factor: f64) -> &mut Self {
        with_ob(|| ffi::constraints_set_factor(self.inner.pin_mut(), factor));
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
