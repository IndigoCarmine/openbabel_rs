//! Ghemical force field energy evaluation.
//!
//! Reproduces `OBForceFieldGhemical`'s `Compute<>` formulas
//! (`forcefieldghemical.cpp`) over the coefficients OpenBabel precomputes:
//! harmonic bonds and angles, a three-term Fourier torsion, a Lennard-Jones
//! 12-6 van der Waals term, and Coulomb electrostatics (no out-of-plane term).
//! Angles are handled in **degrees**, matching OpenBabel (its `theta0` and the
//! `(θ − θ0)` are in degrees). Energies are in kJ/mol.

use super::geom::{angle, dihedral, distance};
use super::{BondTerm, ElecTerm, EnergyModel};

const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
struct AngleTerm {
    a: usize,
    b: usize,
    c: usize,
    ka: f64,
    theta0_deg: f64,
}

#[derive(Clone, Copy, Debug)]
struct TorsionTerm {
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    k1: f64,
    k2: f64,
    k3: f64,
}

#[derive(Clone, Copy, Debug)]
struct VdwTerm {
    a: usize,
    b: usize,
    sigma12: f64,
    sigma6: f64,
}

/// A Ghemical potential-energy surface for one molecule.
pub(crate) struct GhemicalModel {
    n_atoms: usize,
    bonds: Vec<BondTerm>,
    angles: Vec<AngleTerm>,
    torsions: Vec<TorsionTerm>,
    vdws: Vec<VdwTerm>,
    elecs: Vec<ElecTerm>,
}

impl GhemicalModel {
    /// Parse the flat buffer produced by the shim's `ff_export_terms` for the
    /// Ghemical force field. `None` if `format_ok` is 0 or the buffer is short.
    pub(crate) fn from_flat(f: &[f64]) -> Option<GhemicalModel> {
        let mut i = 0usize;
        macro_rules! g {
            () => {{
                let v = *f.get(i)?;
                i += 1;
                v
            }};
        }
        if g!() < 0.5 {
            return None;
        }
        let n_atoms = g!() as usize;
        let mut m = GhemicalModel {
            n_atoms,
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            vdws: Vec::new(),
            elecs: Vec::new(),
        };
        let n = g!() as usize;
        for _ in 0..n {
            m.bonds.push(BondTerm { a: g!() as usize, b: g!() as usize, kb: g!(), r0: g!() });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.angles.push(AngleTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                ka: g!(),
                theta0_deg: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.torsions.push(TorsionTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                k1: g!(),
                k2: g!(),
                k3: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.vdws.push(VdwTerm { a: g!() as usize, b: g!() as usize, sigma12: g!(), sigma6: g!() });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.elecs.push(ElecTerm { a: g!() as usize, b: g!() as usize, qq: g!() });
        }
        Some(m)
    }
}

impl EnergyModel for GhemicalModel {
    fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    fn energy(&self, coords: &[f64]) -> f64 {
        let mut e = 0.0;
        for t in &self.bonds {
            let d = distance(coords, t.a, t.b) - t.r0;
            e += t.kb * d * d;
        }
        for t in &self.angles {
            let theta_deg = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG;
            let d = theta_deg - t.theta0_deg;
            e += t.ka * d * d;
        }
        for t in &self.torsions {
            let phi = dihedral(coords, t.a, t.b, t.c, t.d);
            e += t.k1 * (1.0 + phi.cos())
                + t.k2 * (1.0 - (2.0 * phi).cos())
                + t.k3 * (1.0 + (3.0 * phi).cos());
        }
        for t in &self.vdws {
            let mut r = distance(coords, t.a, t.b);
            if r < 1.0e-4 {
                r = 1.0e-4;
            }
            let a = t.sigma12 / r;
            let b = t.sigma6 / r;
            e += a.powi(12) - b.powi(6);
        }
        for t in &self.elecs {
            let mut r = distance(coords, t.a, t.b);
            if r.abs() < 1.0e-3 {
                r = 1.0e-3;
            }
            e += t.qq / r;
        }
        e
    }
}
