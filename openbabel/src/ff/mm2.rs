//! MM2 force field energy evaluation.
//!
//! **Status: present but unwired.** OpenBabel 3.2.1 does *not* build MM2 —
//! `forcefieldmm2.cpp` is excluded from its `src/CMakeLists.txt` `forcefields`
//! set — so `find_forcefield("MM2")` returns null, `export_flat("MM2")` yields
//! no terms, and [`build_model`](super::build_model) never constructs an
//! [`Mm2Model`] in practice. This evaluator (and the shim's MM2 export path and
//! header accessors) are a complete, ready port kept in the tree: reviving MM2
//! only requires compiling `forcefieldmm2.cpp` back into OpenBabel (and giving
//! its energy methods overrides of the base virtuals, which upstream never did).
//! The five force fields OpenBabel 3.2.1 actually ships — UFF, MMFF94, MMFF94s,
//! GAFF, Ghemical — are all wired and parity-verified.
//!
//! Reproduces `OBForceFieldMM2`'s energy methods (`forcefieldmm2.cpp`). MM2
//! predates OpenBabel's `OBFFCalculation` architecture: it has no precomputed
//! calc vectors and resolves each term's parameters inline, so the shim
//! replicates MM2's own iteration to emit resolved, geometry-independent
//! coefficients and this module applies the functional forms:
//!
//! * **bond** — cubic/quartic stretch `bondunit·k·δ²·(1 + c₃·δ + c₄·δ²)`;
//! * **angle** — sextic bend `angleunit·k·Δ²·(1 + c₆·Δ⁴)` (degrees);
//! * **stretch-bend** — `strbndunit·k·Δθ·(Δr_ab + Δr_bc)`;
//! * **torsion** — `torsionunit·(V1·(1+cosφ) + V2·(1−cos2φ) + V3·(1+cos3φ))`;
//! * **out-of-plane** — `oopunit·k·χ²`, χ the point-to-plane angle in degrees;
//! * **van der Waals** — Buckingham exp-6 `ε·(A·e^(−B·r/r₀) − C·(r₀/r)⁶)`;
//! * **electrostatic** — bond **dipole–dipole** interaction.
//!
//! The dipole term reproduces OpenBabel's expression verbatim, including its
//! use of bond 2's *absolute* midpoint (an origin-dependence quirk in the
//! upstream code) so that energies match at a shared coordinate origin.
//! Energies are in kcal/mol.

use super::geom::{angle, dihedral, distance, dot, point2plane_angle, pos, sub};
use super::EnergyModel;

const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
struct BondTerm {
    a: usize,
    b: usize,
    force: f64,
    l_ref: f64,
}

#[derive(Clone, Copy, Debug)]
struct AngleTerm {
    a: usize,
    b: usize,
    c: usize,
    force: f64,
    ang_ref_deg: f64,
}

#[derive(Clone, Copy, Debug)]
struct StrBndTerm {
    a: usize,
    b: usize,
    c: usize,
    force: f64,
    ang_ref_deg: f64,
    l_ref1: f64,
    l_ref2: f64,
}

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

/// Out-of-plane term; the four atoms are in OpenBabel's `Point2PlaneAngle`
/// argument order (`p0` above the plane of `p1, p2, p3`).
#[derive(Clone, Copy, Debug)]
struct OopTerm {
    p0: usize,
    p1: usize,
    p2: usize,
    p3: usize,
    force: f64,
}

/// Buckingham exp-6 van der Waals. `rr` is the summed vdW radius, `eps` the well
/// depth (both resolved per pair by the shim).
#[derive(Clone, Copy, Debug)]
struct VdwTerm {
    a: usize,
    b: usize,
    rr: f64,
    eps: f64,
}

/// Bond dipole–dipole pair: bond 1 is `a–b` (moment `dipole1`), bond 2 is `c–d`
/// (moment `dipole2`).
#[derive(Clone, Copy, Debug)]
struct DipoleTerm {
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    dipole1: f64,
    dipole2: f64,
}

/// An MM2 potential-energy surface for one molecule.
pub(crate) struct Mm2Model {
    n_atoms: usize,
    bondunit: f64,
    bond_cubic: f64,
    bond_quartic: f64,
    angleunit: f64,
    angle_sextic: f64,
    stretchbendunit: f64,
    torsionunit: f64,
    outplanebendunit: f64,
    a_expterm: f64,
    b_expterm: f64,
    c_expterm: f64,
    elec_f: f64,
    bonds: Vec<BondTerm>,
    angles: Vec<AngleTerm>,
    strbnds: Vec<StrBndTerm>,
    torsions: Vec<TorsionTerm>,
    oops: Vec<OopTerm>,
    vdws: Vec<VdwTerm>,
    dipoles: Vec<DipoleTerm>,
}

impl Mm2Model {
    /// Parse the flat buffer produced by the shim's `ff_export_terms` for MM2.
    /// `None` if `format_ok` is 0 or the buffer is short.
    pub(crate) fn from_flat(f: &[f64]) -> Option<Mm2Model> {
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
        let mut m = Mm2Model {
            n_atoms,
            bondunit: g!(),
            bond_cubic: g!(),
            bond_quartic: g!(),
            angleunit: g!(),
            angle_sextic: g!(),
            stretchbendunit: g!(),
            torsionunit: g!(),
            outplanebendunit: g!(),
            a_expterm: g!(),
            b_expterm: g!(),
            c_expterm: g!(),
            elec_f: g!(),
            bonds: Vec::new(),
            angles: Vec::new(),
            strbnds: Vec::new(),
            torsions: Vec::new(),
            oops: Vec::new(),
            vdws: Vec::new(),
            dipoles: Vec::new(),
        };
        let n = g!() as usize;
        for _ in 0..n {
            m.bonds.push(BondTerm { a: g!() as usize, b: g!() as usize, force: g!(), l_ref: g!() });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.angles.push(AngleTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                force: g!(),
                ang_ref_deg: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.strbnds.push(StrBndTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                force: g!(),
                ang_ref_deg: g!(),
                l_ref1: g!(),
                l_ref2: g!(),
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
                p0: g!() as usize,
                p1: g!() as usize,
                p2: g!() as usize,
                p3: g!() as usize,
                force: g!(),
            });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.vdws.push(VdwTerm { a: g!() as usize, b: g!() as usize, rr: g!(), eps: g!() });
        }
        let n = g!() as usize;
        for _ in 0..n {
            m.dipoles.push(DipoleTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                dipole1: g!(),
                dipole2: g!(),
            });
        }
        Some(m)
    }
}

impl Mm2Model {
    /// Per-component energies `[bond, angle, strbnd, torsion, oop, vdw, elec]`
    /// at `coords`, matching OpenBabel's MM2 `E_*` breakdown order. Used in
    /// tests to pinpoint which term diverges from OpenBabel.
    #[cfg(test)]
    pub(crate) fn energy_components(&self, coords: &[f64]) -> [f64; 7] {
        let mut c = [0.0f64; 7];
        for t in &self.bonds {
            let d = distance(coords, t.a, t.b) - t.l_ref;
            let d2 = d * d;
            c[0] +=
                self.bondunit * t.force * d2 * (1.0 + self.bond_cubic * d + self.bond_quartic * d2);
        }
        for t in &self.angles {
            let d = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG - t.ang_ref_deg;
            let d2 = d * d;
            c[1] += self.angleunit * t.force * d2 * (1.0 + self.angle_sextic * d2 * d2);
        }
        for t in &self.strbnds {
            let d_theta = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG - t.ang_ref_deg;
            let d_ab = distance(coords, t.a, t.b) - t.l_ref1;
            let d_bc = distance(coords, t.b, t.c) - t.l_ref2;
            c[2] += self.stretchbendunit * t.force * d_theta * (d_ab + d_bc);
        }
        for t in &self.torsions {
            let phi = dihedral(coords, t.a, t.b, t.c, t.d);
            c[3] += self.torsionunit
                * (t.v1 * (1.0 + phi.cos())
                    + t.v2 * (1.0 - (2.0 * phi).cos())
                    + t.v3 * (1.0 + (3.0 * phi).cos()));
        }
        for t in &self.oops {
            let chi = point2plane_angle(coords, t.p0, t.p1, t.p2, t.p3);
            c[4] += self.outplanebendunit * t.force * chi * chi;
        }
        for t in &self.vdws {
            let rab = distance(coords, t.a, t.b);
            let rrab = t.rr / rab;
            let rrab6 = rrab.powi(6);
            let abrr = rab / t.rr;
            c[5] +=
                t.eps * (self.a_expterm * (-self.b_expterm * abrr).exp() - self.c_expterm * rrab6);
        }
        for t in &self.dipoles {
            let va = pos(coords, t.a);
            let vb = pos(coords, t.b);
            let vc = pos(coords, t.c);
            let vd = pos(coords, t.d);
            let ab = sub(va, vb);
            let cd = sub(vc, vd);
            let ri2 = dot(ab, ab);
            let rk2 = dot(cd, cd);
            let r = [(vc[0] + vd[0]) / 2.0, (vc[1] + vd[1]) / 2.0, (vc[2] + vd[2]) / 2.0];
            let r2 = dot(r, r);
            let rirkr3 = (ri2 * rk2 * r2).sqrt() * r2;
            let dotp = dot(ab, cd);
            let doti = dot(ab, r);
            let dotk = dot(cd, r);
            let fik = self.elec_f * t.dipole1 * t.dipole2;
            c[6] += fik * (dotp - 3.0 * doti * dotk / r2) / rirkr3;
        }
        c
    }
}

impl EnergyModel for Mm2Model {
    fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    fn energy(&self, coords: &[f64]) -> f64 {
        let mut e = 0.0;

        for t in &self.bonds {
            let d = distance(coords, t.a, t.b) - t.l_ref;
            let d2 = d * d;
            e += self.bondunit * t.force * d2 * (1.0 + self.bond_cubic * d + self.bond_quartic * d2);
        }

        for t in &self.angles {
            let d = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG - t.ang_ref_deg;
            let d2 = d * d;
            e += self.angleunit * t.force * d2 * (1.0 + self.angle_sextic * d2 * d2);
        }

        for t in &self.strbnds {
            let d_theta = angle(coords, t.a, t.b, t.c) * RAD_TO_DEG - t.ang_ref_deg;
            let d_ab = distance(coords, t.a, t.b) - t.l_ref1;
            let d_bc = distance(coords, t.b, t.c) - t.l_ref2;
            e += self.stretchbendunit * t.force * d_theta * (d_ab + d_bc);
        }

        for t in &self.torsions {
            let phi = dihedral(coords, t.a, t.b, t.c, t.d);
            e += self.torsionunit
                * (t.v1 * (1.0 + phi.cos())
                    + t.v2 * (1.0 - (2.0 * phi).cos())
                    + t.v3 * (1.0 + (3.0 * phi).cos()));
        }

        for t in &self.oops {
            let chi = point2plane_angle(coords, t.p0, t.p1, t.p2, t.p3);
            e += self.outplanebendunit * t.force * chi * chi;
        }

        for t in &self.vdws {
            let rab = distance(coords, t.a, t.b);
            let rrab = t.rr / rab;
            let rrab6 = rrab.powi(6);
            let abrr = rab / t.rr;
            e += t.eps * (self.a_expterm * (-self.b_expterm * abrr).exp() - self.c_expterm * rrab6);
        }

        // Bond dipole–dipole. `r` is bond 2's absolute midpoint, matching
        // OpenBabel's (origin-dependent) formulation.
        for t in &self.dipoles {
            let va = pos(coords, t.a);
            let vb = pos(coords, t.b);
            let vc = pos(coords, t.c);
            let vd = pos(coords, t.d);
            let ab = sub(va, vb);
            let cd = sub(vc, vd);
            let ri2 = dot(ab, ab);
            let rk2 = dot(cd, cd);
            let r = [(vc[0] + vd[0]) / 2.0, (vc[1] + vd[1]) / 2.0, (vc[2] + vd[2]) / 2.0];
            let r2 = dot(r, r);
            let rirkr3 = (ri2 * rk2 * r2).sqrt() * r2;
            let dotp = dot(ab, cd);
            let doti = dot(ab, r);
            let dotk = dot(cd, r);
            let fik = self.elec_f * t.dipole1 * t.dipole2;
            e += fik * (dotp - 3.0 * doti * dotk / r2) / rirkr3;
        }

        e
    }
}
