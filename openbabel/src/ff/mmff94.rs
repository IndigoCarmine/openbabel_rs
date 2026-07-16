//! MMFF94 / MMFF94s force field energy evaluation.
//!
//! Reproduces `OBForceFieldMMFF94`'s `Compute<>` formulas
//! (`forcefieldmmff94.cpp`) over the coefficients OpenBabel precomputes. MMFF94
//! and MMFF94s share this class and these functional forms; they differ only in
//! the out-of-plane / torsion parameters chosen at setup, which are baked into
//! the exported coefficients, so one evaluator serves both.
//!
//! The energy terms are the MMFF94 set:
//! * **bond** — quartic stretch `143.9325·½·kb·δ²·(1 − 2δ + 7⁄3·δ²)`;
//! * **angle** — cubic bend `0.043844·½·ka·Δ²·(1 − 0.007·Δ)` (degrees), or a
//!   cosine form for designated *linear* angles;
//! * **stretch-bend** — `2.51210·(kbaABC·Δrab + kbaCBA·Δrbc)·Δθ` (θ in degrees);
//! * **torsion** — `½·(V1·(1+cosφ) + V2·(1−cos2φ) + V3·(1+cos3φ))`;
//! * **out-of-plane** — `0.043844·½·koop·χ²`, χ the Wilson angle in degrees;
//! * **van der Waals** — the buffered 14-7 potential `ε·erep⁷·eattr`;
//! * **electrostatic** — buffered Coulomb `qq / (r + 0.05)`.
//!
//! The per-term unit factors OpenBabel applies once to each energy sum (e.g.
//! `143.9325·½` for bonds) are folded into each term here so the total is a
//! plain sum. Energies are in kcal/mol (MMFF94's native unit).

use super::geom::{angle, dihedral, distance, wilson};
use super::{ElecTerm, EnergyModel};

const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

/// Quartic bond stretch. `kb` md/Å², `r0` Å.
#[derive(Clone, Copy, Debug)]
struct BondTerm {
    a: usize,
    b: usize,
    kb: f64,
    r0: f64,
}

/// Angle bend. `linear` angles use the cosine form; others the cubic form.
/// `theta0` is in degrees.
#[derive(Clone, Copy, Debug)]
struct AngleTerm {
    a: usize,
    b: usize,
    c: usize,
    ka: f64,
    theta0_deg: f64,
    linear: bool,
}

/// Stretch-bend coupling: bond stretches `a–b`, `b–c` coupled to the `a–b–c`
/// angle. `theta0` degrees, `rab0`/`rbc0` Å.
#[derive(Clone, Copy, Debug)]
struct StrBndTerm {
    a: usize,
    b: usize,
    c: usize,
    kba_abc: f64,
    kba_cba: f64,
    theta0_deg: f64,
    rab0: f64,
    rbc0: f64,
}

/// Threefold Fourier torsion `½·(V1·φ1 + V2·φ2 + V3·φ3)`.
#[derive(Clone, Copy, Debug)]
struct TorsionTerm {
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    v1: f64,
    v2: f64,
    v3: f64,
}

/// Out-of-plane (Wilson) bend at central atom `b`, apex `d`.
#[derive(Clone, Copy, Debug)]
struct OopTerm {
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    koop: f64,
}

/// Buffered 14-7 van der Waals. `r_ab` is the combined minimum-energy distance,
/// `r_ab7 = r_ab⁷` (precomputed by OpenBabel), `epsilon` the well depth.
#[derive(Clone, Copy, Debug)]
struct VdwTerm {
    a: usize,
    b: usize,
    epsilon: f64,
    r_ab: f64,
    r_ab7: f64,
}

/// An MMFF94 (or MMFF94s) potential-energy surface for one molecule.
pub(crate) struct Mmff94Model {
    n_atoms: usize,
    bonds: Vec<BondTerm>,
    angles: Vec<AngleTerm>,
    strbnds: Vec<StrBndTerm>,
    torsions: Vec<TorsionTerm>,
    oops: Vec<OopTerm>,
    vdws: Vec<VdwTerm>,
    elecs: Vec<ElecTerm>,
}

impl Mmff94Model {
    /// Parse the flat buffer produced by the shim's `ff_export_terms` for the
    /// MMFF94/MMFF94s force field. `None` if `format_ok` is 0 or the buffer is
    /// short.
    pub(crate) fn from_flat(f: &[f64]) -> Option<Mmff94Model> {
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
        let mut m = Mmff94Model {
            n_atoms,
            bonds: Vec::new(),
            angles: Vec::new(),
            strbnds: Vec::new(),
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
                ka: g!(),
                theta0_deg: g!(),
                linear: g!() > 0.5,
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.strbnds.push(StrBndTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                kba_abc: g!(),
                kba_cba: g!(),
                theta0_deg: g!(),
                rab0: g!(),
                rbc0: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.torsions.push(TorsionTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                v1: g!(),
                v2: g!(),
                v3: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.oops.push(OopTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                koop: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.vdws.push(VdwTerm {
                a: g!() as usize,
                b: g!() as usize,
                epsilon: g!(),
                r_ab: g!(),
                r_ab7: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.elecs.push(ElecTerm { a: g!() as usize, b: g!() as usize, qq: g!() });
        }
        Some(m)
    }
}

impl Mmff94Model {
    /// Per-component energies `[bond, angle, strbnd, torsion, oop, vdw, elec]`
    /// at `coords`, matching OpenBabel's `E_*` breakdown order. Used in tests to
    /// pinpoint which term diverges from OpenBabel.
    #[cfg(test)]
    pub(crate) fn energy_components(&self, coords: &[f64]) -> [f64; 7] {
        let mut c = [0.0f64; 7];
        for t in &self.bonds {
            let d = distance(coords, t.a, t.b) - t.r0;
            let d2 = d * d;
            c[0] += 143.9325 * 0.5 * t.kb * d2 * (1.0 - 2.0 * d + 7.0 / 3.0 * d2);
        }
        for t in &self.angles {
            let theta_rad = angle(coords, t.a, t.b, t.c);
            if t.linear {
                c[1] += 143.9325 * t.ka * (1.0 + theta_rad.cos());
            } else {
                let d = theta_rad * RAD_TO_DEG - t.theta0_deg;
                c[1] += 0.043844 * 0.5 * t.ka * d * d * (1.0 - 0.007 * d);
            }
        }
        for t in &self.strbnds {
            let theta_deg = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG;
            let d_rab = distance(coords, t.a, t.b) - t.rab0;
            let d_rbc = distance(coords, t.b, t.c) - t.rbc0;
            let d_theta = theta_deg - t.theta0_deg;
            c[2] += 2.51210 * (t.kba_abc * d_rab + t.kba_cba * d_rbc) * d_theta;
        }
        for t in &self.torsions {
            let phi = dihedral(coords, t.a, t.b, t.c, t.d);
            c[3] += 0.5
                * (t.v1 * (1.0 + phi.cos())
                    + t.v2 * (1.0 - (2.0 * phi).cos())
                    + t.v3 * (1.0 + (3.0 * phi).cos()));
        }
        for t in &self.oops {
            let chi_deg = wilson(coords, t.a, t.b, t.c, t.d) * RAD_TO_DEG;
            c[4] += 0.043844 * 0.5 * t.koop * chi_deg * chi_deg;
        }
        for t in &self.vdws {
            let rab = distance(coords, t.a, t.b);
            let rab7 = rab.powi(7);
            let erep = (1.07 * t.r_ab) / (rab + 0.07 * t.r_ab);
            let erep7 = erep.powi(7);
            let eattr = ((1.12 * t.r_ab7) / (rab7 + 0.12 * t.r_ab7)) - 2.0;
            c[5] += t.epsilon * erep7 * eattr;
        }
        for t in &self.elecs {
            let rab = distance(coords, t.a, t.b) + 0.05;
            c[6] += t.qq / rab;
        }
        c
    }
}

impl EnergyModel for Mmff94Model {
    fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    fn energy(&self, coords: &[f64]) -> f64 {
        let mut e = 0.0;

        // Bond: quartic stretch, folded 143.9325·½ unit factor.
        for t in &self.bonds {
            let d = distance(coords, t.a, t.b) - t.r0;
            let d2 = d * d;
            e += 143.9325 * 0.5 * t.kb * d2 * (1.0 - 2.0 * d + 7.0 / 3.0 * d2);
        }

        // Angle: cubic bend (degrees), or cosine form for linear angles.
        for t in &self.angles {
            let theta_rad = angle(coords, t.a, t.b, t.c);
            if t.linear {
                e += 143.9325 * t.ka * (1.0 + theta_rad.cos());
            } else {
                let d = theta_rad * RAD_TO_DEG - t.theta0_deg;
                e += 0.043844 * 0.5 * t.ka * d * d * (1.0 - 0.007 * d);
            }
        }

        // Stretch-bend: bond/angle coupling, folded 2.51210 unit factor. The
        // DEG_TO_RAD·RAD_TO_DEG of OpenBabel's expression cancels, leaving the
        // angle delta in degrees.
        for t in &self.strbnds {
            let theta_deg = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG;
            let d_rab = distance(coords, t.a, t.b) - t.rab0;
            let d_rbc = distance(coords, t.b, t.c) - t.rbc0;
            let d_theta = theta_deg - t.theta0_deg;
            e += 2.51210 * (t.kba_abc * d_rab + t.kba_cba * d_rbc) * d_theta;
        }

        // Torsion: threefold Fourier, folded ½ unit factor.
        for t in &self.torsions {
            let phi = dihedral(coords, t.a, t.b, t.c, t.d);
            let phi1 = 1.0 + phi.cos();
            let phi2 = 1.0 - (2.0 * phi).cos();
            let phi3 = 1.0 + (3.0 * phi).cos();
            e += 0.5 * (t.v1 * phi1 + t.v2 * phi2 + t.v3 * phi3);
        }

        // Out-of-plane: Wilson angle in degrees, folded 0.043844·½ unit factor.
        for t in &self.oops {
            let chi_deg = wilson(coords, t.a, t.b, t.c, t.d) * RAD_TO_DEG;
            e += 0.043844 * 0.5 * t.koop * chi_deg * chi_deg;
        }

        // Van der Waals: buffered 14-7.
        for t in &self.vdws {
            let rab = distance(coords, t.a, t.b);
            let rab7 = rab.powi(7);
            let erep = (1.07 * t.r_ab) / (rab + 0.07 * t.r_ab);
            let erep7 = erep.powi(7);
            let eattr = ((1.12 * t.r_ab7) / (rab7 + 0.12 * t.r_ab7)) - 2.0;
            e += t.epsilon * erep7 * eattr;
        }

        // Electrostatic: buffered Coulomb (qq already includes MMFF's 332.0716).
        for t in &self.elecs {
            let rab = distance(coords, t.a, t.b) + 0.05;
            e += t.qq / rab;
        }

        e
    }
}
