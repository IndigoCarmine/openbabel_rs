//! Rust force-field numeric core.
//!
//! OpenBabel does the chemistry (perception, atom typing, and precomputing each
//! energy term's coefficients) once, under the global lock; the resulting term
//! list is exported into the [`Terms`] structs here. This module then evaluates
//! the energy and gradient and runs the minimizer **without touching OpenBabel**
//! — pure arithmetic over the term list and a coordinate buffer — so the hot
//! optimization loop carries no global state and can run in parallel.
//!
//! Each force field supplies its own functional forms in a submodule (e.g.
//! [`uff`]); the coefficients in the term structs are force-field specific (the
//! `form`/`coord` discriminants say how to read them). Coordinates are flat
//! `[x0,y0,z0, x1,…]` triples in atom-index order.
#![allow(dead_code)] // wired incrementally across P1–P5

pub(crate) mod gaff;
pub(crate) mod geom;
pub(crate) mod ghemical;
pub(crate) mod minimize;
pub(crate) mod mm2;
pub(crate) mod mmff94;
pub(crate) mod uff;

/// Harmonic bond stretch `kb·(r − r0)²`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BondTerm {
    pub a: usize,
    pub b: usize,
    pub kb: f64,
    pub r0: f64,
}

/// Angle bend. The `coord` discriminant selects the functional form (see
/// [`uff`]); `ka`, `theta0` (radians), and `c0..c2` are the precomputed
/// coefficients, `n` the periodicity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct AngleTerm {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub ka: f64,
    pub theta0: f64,
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub coord: i32,
    pub n: i32,
}

/// Torsion `V·(1 − cosNPhi0·cos(n·φ))`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TorsionTerm {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub d: usize,
    pub v: f64,
    pub n: i32,
    pub cos_nphi0: f64,
}

/// Out-of-plane bend `koop·(c0 + c1·cosψ + c2·cos2ψ)`, `b` central.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct OopTerm {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub d: usize,
    pub koop: f64,
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
}

/// Van der Waals (UFF: Lennard-Jones 12-6). `ka_squared = x_ij²` (combined
/// equilibrium distance squared), `kab` the well depth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct VdwTerm {
    pub a: usize,
    pub b: usize,
    pub kab: f64,
    pub ka_squared: f64,
}

/// Electrostatic `qq / r` (disabled by default in UFF).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ElecTerm {
    pub a: usize,
    pub b: usize,
    pub qq: f64,
}

/// A complete set of precomputed energy terms for one molecule + force field,
/// as exported from OpenBabel's force-field setup.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Terms {
    pub n_atoms: usize,
    pub bonds: Vec<BondTerm>,
    pub angles: Vec<AngleTerm>,
    pub torsions: Vec<TorsionTerm>,
    pub oops: Vec<OopTerm>,
    pub vdws: Vec<VdwTerm>,
    pub elecs: Vec<ElecTerm>,
}

impl Terms {
    /// Parse the flat, self-describing `f64` buffer produced by the shim's
    /// `ff_export_terms`. Returns `None` if the leading `format_ok` flag is 0
    /// (unknown/unsupported force field) or the buffer is truncated.
    pub(crate) fn from_flat(f: &[f64]) -> Option<Terms> {
        let mut i = 0usize;
        // Pull the next value, bailing (via `?`) if the buffer is too short.
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
        let mut t = Terms { n_atoms, ..Default::default() };

        let n = g!() as usize;
        t.bonds.reserve(n);
        for _ in 0..n {
            t.bonds.push(BondTerm { a: g!() as usize, b: g!() as usize, kb: g!(), r0: g!() });
        }
        let n = g!() as usize;
        t.angles.reserve(n);
        for _ in 0..n {
            t.angles.push(AngleTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                ka: g!(),
                theta0: g!(),
                c0: g!(),
                c1: g!(),
                c2: g!(),
                coord: g!() as i32,
                n: g!() as i32,
            });
        }
        let n = g!() as usize;
        t.torsions.reserve(n);
        for _ in 0..n {
            t.torsions.push(TorsionTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                v: g!(),
                n: g!() as i32,
                cos_nphi0: g!(),
            });
        }
        let n = g!() as usize;
        t.oops.reserve(n);
        for _ in 0..n {
            t.oops.push(OopTerm {
                a: g!() as usize,
                b: g!() as usize,
                c: g!() as usize,
                d: g!() as usize,
                koop: g!(),
                c0: g!(),
                c1: g!(),
                c2: g!(),
            });
        }
        let n = g!() as usize;
        t.vdws.reserve(n);
        for _ in 0..n {
            t.vdws.push(VdwTerm { a: g!() as usize, b: g!() as usize, kab: g!(), ka_squared: g!() });
        }
        let n = g!() as usize;
        t.elecs.reserve(n);
        for _ in 0..n {
            t.elecs.push(ElecTerm { a: g!() as usize, b: g!() as usize, qq: g!() });
        }
        Some(t)
    }
}

/// A potential-energy surface over a flat coordinate buffer.
pub(crate) trait EnergyModel {
    /// Number of atoms; the coordinate buffer must have length `3 * n_atoms()`.
    fn n_atoms(&self) -> usize;

    /// Total potential energy at `coords` (force field's own energy unit).
    fn energy(&self, coords: &[f64]) -> f64;

    /// Fill `grad` (length `3·n_atoms`) with `∂E/∂coords` and return the energy.
    ///
    /// The default uses central finite differences: correct for any `energy`
    /// implementation, but `O(n)` energy evaluations per call. Force fields
    /// override this with analytic gradients for speed.
    fn energy_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        fd_gradient(self, coords, grad);
        self.energy(coords)
    }
}

/// Build the energy model for `forcefield` from the flat term buffer the shim
/// exported, or `None` if no Rust evaluator exists for it yet (so the caller
/// falls back to OpenBabel). The dispatch point where new force fields are
/// added; each force field parses its own layout.
pub(crate) fn build_model(forcefield: &str, flat: &[f64]) -> Option<Box<dyn EnergyModel>> {
    match forcefield.to_ascii_uppercase().as_str() {
        "UFF" => Terms::from_flat(flat).map(|t| Box::new(uff::UffModel::new(t)) as Box<dyn EnergyModel>),
        "GHEMICAL" => {
            ghemical::GhemicalModel::from_flat(flat).map(|m| Box::new(m) as Box<dyn EnergyModel>)
        }
        "GAFF" => gaff::GaffModel::from_flat(flat).map(|m| Box::new(m) as Box<dyn EnergyModel>),
        "MMFF94" | "MMFF94S" => {
            mmff94::Mmff94Model::from_flat(flat).map(|m| Box::new(m) as Box<dyn EnergyModel>)
        }
        "MM2" => mm2::Mm2Model::from_flat(flat).map(|m| Box::new(m) as Box<dyn EnergyModel>),
        _ => None,
    }
}

/// Central finite-difference gradient of `m.energy` at `coords` into `grad`.
///
/// Used as the default gradient and as the reference an analytic gradient is
/// validated against in tests.
pub(crate) fn fd_gradient<M: EnergyModel + ?Sized>(m: &M, coords: &[f64], grad: &mut [f64]) {
    let mut x = coords.to_vec();
    let h = 1e-5;
    for i in 0..x.len() {
        let xi = x[i];
        x[i] = xi + h;
        let ep = m.energy(&x);
        x[i] = xi - h;
        let em = m.energy(&x);
        x[i] = xi;
        grad[i] = (ep - em) / (2.0 * h);
    }
}

#[cfg(test)]
mod parity {
    //! The Rust numeric core evaluated over OpenBabel-exported terms must
    //! reproduce OpenBabel's own force-field energy.
    use super::build_model;
    use super::minimize::minimize;
    use crate::{Algorithm, Molecule};

    fn built_3d(smiles: &str) -> Molecule {
        let mut m = Molecule::parse(smiles, "smi").unwrap();
        m.add_hydrogens();
        assert!(m.generate_3d(), "3D generation failed for {smiles}");
        m
    }

    /// Serialize the OpenBabel-reference parity tests against one another.
    ///
    /// `built_3d`, `export_flat`, `energy`, and `optimize_geometry` all drive
    /// OpenBabel's shared, non-thread-safe force-field singleton. Although every
    /// FFI call is individually locked, OpenBabel's `Setup` keeps an order-
    /// dependent cross-call cache (`_mol` plus the precomputed calc vectors), so
    /// interleaving two of these multi-call test bodies — which `cargo test`
    /// does by default — makes OpenBabel reuse stale non-bonded terms and return
    /// inconsistent energies for the *same* molecule. That non-thread-safety is
    /// exactly what this module's Rust core exists to bypass; these tests only
    /// validate the Rust core against OpenBabel's *sequential* ground truth, so
    /// they hold this guard for their full duration. (The Rust optimize path is
    /// itself lock-free and gets a real concurrency test elsewhere.)
    fn ob_guard() -> std::sync::MutexGuard<'static, ()> {
        static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
        LOCK.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// The Rust single-point energy over OpenBabel-exported terms equals
    /// OpenBabel's own force-field energy at the same geometry.
    fn assert_energy_matches(smiles: &str, ff: &str) {
        let m = built_3d(smiles);
        let model = build_model(ff, &m.export_flat(ff)).expect("Rust model");
        let coords = m.coordinates();
        let rust_e = model.energy(&coords);
        let ob_e = m.energy(ff).expect("OpenBabel energy");
        let tol = 1e-3 * (1.0 + ob_e.abs());
        assert!(
            (rust_e - ob_e).abs() < tol,
            "{smiles}/{ff}: Rust {rust_e} vs OpenBabel {ob_e} (tol {tol})"
        );
    }

    /// The Rust minimizer reduces the energy and lands near OpenBabel's own
    /// minimized energy from the same starting geometry.
    ///
    /// The built geometry is perturbed first: `generate_3d` finishes with an
    /// MMFF94 cleanup, so an un-perturbed structure can already sit at (or in) a
    /// force field's minimum with nothing to descend. A small, deterministic
    /// displacement — far too small to change perception or atom typing —
    /// restores a real basin for every force field, and is written back into the
    /// molecule so OpenBabel's minimizer starts from the identical geometry.
    fn assert_minimizer_matches(smiles: &str, ff: &str) {
        let mut m = built_3d(smiles);
        let mut start = m.coordinates();
        for (i, x) in start.iter_mut().enumerate() {
            *x += 0.1 * ((i as f64) * 1.7).sin();
        }
        assert!(m.set_coordinates(&start), "set_coordinates failed");

        let model = build_model(ff, &m.export_flat(ff)).expect("Rust model");
        let e_start = model.energy(&start);
        let out = minimize(model.as_ref(), &start, Algorithm::ConjugateGradients, 2000, 1e-8, 1);
        assert!(out.energy < e_start - 1.0, "{ff}: no reduction {e_start} -> {}", out.energy);

        let mut ob = m.clone();
        let ob_e = ob.optimize_geometry(ff, 2000).expect("OB minimize");
        // Both descend the same PES from the same start; allow a few units of
        // energy for differing line-search / convergence paths.
        assert!(
            (out.energy - ob_e).abs() < 5.0,
            "{smiles}/{ff}: Rust min {} vs OpenBabel min {ob_e}",
            out.energy
        );
    }

    #[test]
    fn uff_energy_matches_openbabel() {
        let _g = ob_guard();
        assert_energy_matches("CC", "UFF");
        assert_energy_matches("CCO", "UFF");
    }

    #[test]
    fn ghemical_energy_matches_openbabel() {
        let _g = ob_guard();
        assert_energy_matches("CC", "Ghemical");
        assert_energy_matches("CCO", "Ghemical");
    }

    #[test]
    fn gaff_energy_matches_openbabel() {
        let _g = ob_guard();
        assert_energy_matches("CC", "GAFF");
        assert_energy_matches("CCO", "GAFF");
    }

    #[test]
    fn mmff94_energy_matches_openbabel() {
        let _g = ob_guard();
        assert_energy_matches("CC", "MMFF94");
        assert_energy_matches("CCO", "MMFF94");
        assert_energy_matches("c1ccccc1", "MMFF94");
    }

    #[test]
    fn mmff94s_energy_matches_openbabel() {
        let _g = ob_guard();
        assert_energy_matches("CCO", "MMFF94s");
        assert_energy_matches("c1ccccc1", "MMFF94s");
    }

    // MM2 is excluded from OpenBabel 3.2.1's build (not in src/CMakeLists.txt's
    // `forcefields` set), so `find_forcefield("MM2")` is null and there is no
    // OpenBabel MM2 to validate against until it is compiled back in.
    #[test]
    #[ignore = "MM2 not compiled into OpenBabel 3.2.1; pending build revival"]
    fn mm2_energy_matches_openbabel() {
        let _g = ob_guard();
        assert_energy_matches("CC", "MM2");
        assert_energy_matches("CCO", "MM2");
    }

    #[test]
    fn uff_minimizer_matches_openbabel() {
        let _g = ob_guard();
        assert_minimizer_matches("CCO", "UFF");
    }

    #[test]
    fn ghemical_minimizer_matches_openbabel() {
        let _g = ob_guard();
        assert_minimizer_matches("CCO", "Ghemical");
    }

    #[test]
    fn gaff_minimizer_matches_openbabel() {
        let _g = ob_guard();
        assert_minimizer_matches("CCO", "GAFF");
    }

    #[test]
    fn mmff94_minimizer_matches_openbabel() {
        let _g = ob_guard();
        assert_minimizer_matches("CCO", "MMFF94");
    }

    #[test]
    #[ignore = "MM2 not compiled into OpenBabel 3.2.1; pending build revival"]
    fn mm2_minimizer_matches_openbabel() {
        let _g = ob_guard();
        assert_minimizer_matches("CCO", "MM2");
    }

    /// Per-component parity on benzene, which (unlike the alkanes) exercises the
    /// stretch-bend, out-of-plane, and aromatic-torsion terms. Compares every
    /// MMFF94 energy component against OpenBabel's own `E_*` breakdown so a
    /// single diverging term is pinpointed rather than hidden in the total.
    #[test]
    fn mmff94_components_match_openbabel() {
        let _g = ob_guard();
        let names = ["bond", "angle", "strbnd", "torsion", "oop", "vdw", "elec"];
        for _ in 0..4 {
            let m = built_3d("c1ccccc1");
            let coords = m.coordinates();
            let model =
                super::mmff94::Mmff94Model::from_flat(&m.export_flat("MMFF94")).expect("model");
            let rust = model.energy_components(&coords);
            let ob = m.energy_components("MMFF94");
            assert!(ob.first().copied().unwrap_or(0.0) > 0.5, "OB components unavailable");
            for (k, name) in names.iter().enumerate() {
                let (r, o) = (rust[k], ob[k + 1]);
                let tol = 1e-2 * (1.0 + o.abs());
                assert!((r - o).abs() < tol, "MMFF94 {name}: Rust {r} vs OB {o} (tol {tol})");
            }
        }
    }
}
