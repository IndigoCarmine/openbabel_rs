//! UFF (Universal Force Field) energy evaluation.
//!
//! Reproduces `OBForceFieldUFF`'s `Compute<>` formulas (`forcefielduff.cpp`)
//! term-for-term, reading the coefficients OpenBabel precomputed into the
//! exported [`Terms`]. Energy is analytic; the gradient currently uses the
//! finite-difference default (a correct baseline — analytic gradients are the
//! next optimization, and OpenBabel's per-term derivative formulas are already
//! captured for the port).
//!
//! UFF energies are in kJ/mol (OpenBabel converts from kcal/mol internally).

use super::geom::{angle, dihedral, distance, distance_sq, wilson};
use super::{EnergyModel, Terms};

/// A UFF potential-energy surface for one molecule: the precomputed term list.
pub(crate) struct UffModel {
    terms: Terms,
}

impl UffModel {
    pub(crate) fn new(terms: Terms) -> Self {
        UffModel { terms }
    }

    fn bond_energy(&self, coords: &[f64]) -> f64 {
        self.terms
            .bonds
            .iter()
            .map(|t| {
                let delta = distance(coords, t.a, t.b) - t.r0;
                t.kb * delta * delta
            })
            .sum()
    }

    fn angle_energy(&self, coords: &[f64]) -> f64 {
        self.terms
            .angles
            .iter()
            .map(|t| {
                let theta = angle(coords, t.a, t.b, t.c);
                let n = t.n as f64;
                match t.coord {
                    // sp: linear, minimum at 180°.
                    1 => t.ka * (1.0 + theta.cos()),
                    // sp2 / square planar / octahedral: periodic + a soft wall
                    // near θ=0 (ESFF-style penalty), exactly as OpenBabel.
                    2 | 4 | 6 => {
                        t.ka * (1.0 - (n * theta).cos())
                            + (-20.0 * (theta - t.theta0 + 0.25)).exp()
                    }
                    // pentagonal bipyramidal (verbatim, including OpenBabel's
                    // literal constants).
                    7 => {
                        let cos_t = theta.cos();
                        t.ka * t.c1
                            * (cos_t - 0.309_016_99)
                            * (cos_t - 0.309_061_99)
                            * (cos_t + 0.809_016_99)
                            * (cos_t + 0.809_169_9)
                    }
                    // general sp3: Fourier expansion, cos2θ = 2cos²θ − 1.
                    _ => {
                        let cos_t = theta.cos();
                        t.ka * (t.c0 + t.c1 * cos_t + t.c2 * (2.0 * cos_t * cos_t - 1.0))
                    }
                }
            })
            .sum()
    }

    fn torsion_energy(&self, coords: &[f64]) -> f64 {
        self.terms
            .torsions
            .iter()
            .map(|t| {
                let tor = dihedral(coords, t.a, t.b, t.c, t.d);
                t.v * (1.0 - t.cos_nphi0 * ((t.n as f64) * tor).cos())
            })
            .sum()
    }

    fn oop_energy(&self, coords: &[f64]) -> f64 {
        self.terms
            .oops
            .iter()
            .map(|t| {
                let psi = wilson(coords, t.a, t.b, t.c, t.d);
                t.koop * (t.c0 + t.c1 * psi.cos() + t.c2 * (2.0 * psi).cos())
            })
            .sum()
    }

    fn vdw_energy(&self, coords: &[f64]) -> f64 {
        self.terms
            .vdws
            .iter()
            .map(|t| {
                // Clamp like OpenBabel so a near-coincident pair can't blow up.
                let mut rab2 = distance_sq(coords, t.a, t.b);
                if rab2 < 1.0e-5 {
                    rab2 = 1.0e-5;
                }
                let term6 = {
                    let s = t.ka_squared / rab2; // (x_ij / r)²
                    s * s * s // ^6
                };
                let term12 = term6 * term6;
                t.kab * (term12 - 2.0 * term6)
            })
            .sum()
    }

    fn elec_energy(&self, coords: &[f64]) -> f64 {
        self.terms
            .elecs
            .iter()
            .map(|t| {
                let mut rab = distance(coords, t.a, t.b);
                if rab.abs() < 1.0e-3 {
                    rab = 1.0e-3;
                }
                t.qq / rab
            })
            .sum()
    }
}

impl EnergyModel for UffModel {
    fn n_atoms(&self) -> usize {
        self.terms.n_atoms
    }

    fn energy(&self, coords: &[f64]) -> f64 {
        self.bond_energy(coords)
            + self.angle_energy(coords)
            + self.torsion_energy(coords)
            + self.oop_energy(coords)
            + self.vdw_energy(coords)
            + self.elec_energy(coords)
    }
}

#[cfg(test)]
mod tests {
    use super::super::{BondTerm, Terms, VdwTerm};
    use super::*;

    fn model(terms: Terms) -> UffModel {
        UffModel::new(terms)
    }

    #[test]
    fn bond_at_equilibrium_is_zero() {
        let terms = Terms {
            n_atoms: 2,
            bonds: vec![BondTerm { a: 0, b: 1, kb: 100.0, r0: 1.5 }],
            ..Default::default()
        };
        let m = model(terms);
        let coords = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        assert!(m.energy(&coords).abs() < 1e-9);
    }

    #[test]
    fn bond_harmonic_matches_formula() {
        let terms = Terms {
            n_atoms: 2,
            bonds: vec![BondTerm { a: 0, b: 1, kb: 100.0, r0: 1.5 }],
            ..Default::default()
        };
        let m = model(terms);
        // stretched to 1.7 → kb·(0.2)² = 100·0.04 = 4.0
        let coords = [0.0, 0.0, 0.0, 1.7, 0.0, 0.0];
        assert!((m.energy(&coords) - 4.0).abs() < 1e-9);
    }

    #[test]
    fn vdw_minimum_near_equilibrium() {
        // LJ minimum at r = x_ij (where ka_squared = x_ij²): energy = -kab.
        let x_ij = 3.5;
        let terms = Terms {
            n_atoms: 2,
            vdws: vec![VdwTerm { a: 0, b: 1, kab: 0.2, ka_squared: x_ij * x_ij }],
            ..Default::default()
        };
        let m = model(terms);
        let at_min = [0.0, 0.0, 0.0, x_ij, 0.0, 0.0];
        assert!((m.energy(&at_min) - (-0.2)).abs() < 1e-9);
        // Slightly compressed/expanded is higher.
        let closer = [0.0, 0.0, 0.0, x_ij - 0.3, 0.0, 0.0];
        assert!(m.energy(&closer) > m.energy(&at_min));
    }

    #[test]
    fn gradient_matches_finite_difference() {
        // A tiny bond+vdw system in a generic geometry; the analytic-or-default
        // gradient must agree with an independent finite difference.
        let terms = Terms {
            n_atoms: 3,
            bonds: vec![
                BondTerm { a: 0, b: 1, kb: 120.0, r0: 1.4 },
                BondTerm { a: 1, b: 2, kb: 90.0, r0: 1.5 },
            ],
            vdws: vec![VdwTerm { a: 0, b: 2, kab: 0.15, ka_squared: 3.2 * 3.2 }],
            ..Default::default()
        };
        let m = model(terms);
        let coords = [0.1, 0.0, 0.0, 1.5, 0.2, 0.0, 2.9, -0.3, 0.1];
        let mut g = vec![0.0; coords.len()];
        m.energy_gradient(&coords, &mut g);
        let mut fd = vec![0.0; coords.len()];
        super::super::fd_gradient(&m, &coords, &mut fd);
        for (a, b) in g.iter().zip(fd.iter()) {
            assert!((a - b).abs() < 1e-6, "grad {a} vs fd {b}");
        }
    }
}
