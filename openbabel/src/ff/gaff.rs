//! GAFF (General Amber Force Field) energy evaluation.
//!
//! Reproduces `OBForceFieldGaff`'s `Compute<>` formulas (`forcefieldgaff.cpp`)
//! over the coefficients OpenBabel precomputes: harmonic bonds and angles, an
//! AMBER cosine torsion and improper (out-of-plane) torsion of the same form, a
//! Lennard-Jones 12-6 van der Waals term in AMBER `R_min` parameterization, and
//! Coulomb electrostatics. Energies are in kJ/mol.
//!
//! Angles and torsions are evaluated in radians here; OpenBabel stores the
//! reference angle `theta0` and torsion phase `gamma` in degrees, so both are
//! converted to radians at load time. OpenBabel's `n·tor − gamma` (with `tor`
//! in degrees) is thus `n·φ − γ` with everything in radians.

use super::geom::{angle, dihedral, distance};
use super::{BondTerm, ElecTerm, EnergyModel};

const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

#[derive(Clone, Copy, Debug)]
struct AngleTerm {
    a: usize,
    b: usize,
    c: usize,
    kth: f64,
    theta0_rad: f64,
}

/// AMBER cosine torsion / improper torsion: `vn_half·(1 + cos(n·φ − γ))`.
/// Both proper torsions and out-of-plane (improper) terms use this form and are
/// evaluated over the dihedral of their four atoms.
#[derive(Clone, Copy, Debug)]
struct TorsionTerm {
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    vn_half: f64,
    gamma_rad: f64,
    n: f64,
}

/// Van der Waals in AMBER `R_min` form: `Eab·((RVDWab/r)¹² − 2(RVDWab/r)⁶)`.
#[derive(Clone, Copy, Debug)]
struct VdwTerm {
    a: usize,
    b: usize,
    eab: f64,
    rvdwab: f64,
}

/// A GAFF potential-energy surface for one molecule.
pub(crate) struct GaffModel {
    n_atoms: usize,
    bonds: Vec<BondTerm>,
    angles: Vec<AngleTerm>,
    torsions: Vec<TorsionTerm>,
    oops: Vec<TorsionTerm>,
    vdws: Vec<VdwTerm>,
    elecs: Vec<ElecTerm>,
}

impl GaffModel {
    /// Parse the flat buffer produced by the shim's `ff_export_terms` for the
    /// GAFF force field. `None` if `format_ok` is 0 or the buffer is short.
    pub(crate) fn from_flat(f: &[f64]) -> Option<GaffModel> {
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
        let mut m = GaffModel {
            n_atoms,
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            oops: Vec::new(),
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
                kth: g!(),
                theta0_rad: g!() * DEG_TO_RAD,
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.torsions.push(TorsionTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                vn_half: g!(),
                gamma_rad: g!() * DEG_TO_RAD,
                n: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.oops.push(TorsionTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                vn_half: g!(),
                gamma_rad: g!() * DEG_TO_RAD,
                n: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.vdws.push(VdwTerm { a: g!() as usize, b: g!() as usize, eab: g!(), rvdwab: g!() });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.elecs.push(ElecTerm { a: g!() as usize, b: g!() as usize, qq: g!() });
        }
        Some(m)
    }
}

impl EnergyModel for GaffModel {
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
            let d = angle(coords, t.a, t.b, t.c) - t.theta0_rad;
            e += t.kth * d * d;
        }
        // Proper torsions and improper (out-of-plane) torsions share the form.
        for t in self.torsions.iter().chain(self.oops.iter()) {
            let phi = dihedral(coords, t.a, t.b, t.c, t.d);
            e += t.vn_half * (1.0 + (t.n * phi - t.gamma_rad).cos());
        }
        for t in &self.vdws {
            let mut r = distance(coords, t.a, t.b);
            if r < 1.0e-4 {
                r = 1.0e-4;
            }
            let term = t.rvdwab / r;
            let term6 = term.powi(6);
            let term12 = term6 * term6;
            e += t.eab * (term12 - 2.0 * term6);
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
