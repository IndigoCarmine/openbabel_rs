//! OpenBabel's analytic gradients must be the derivative of its energy.
//!
//! The minimizers step along the gradient and test for convergence on it, so a
//! gradient that disagrees with the energy sends them in a direction that does
//! not lead downhill: the line search finds nothing and the run stops where it
//! is, reporting convergence. Each test here compares
//! [`Molecule::energy_gradient`]'s analytic gradient with a central difference
//! of its own energy, over geometries kicked off their minimum so that every
//! term carries a real force.

use openbabel::{Constraints, Molecule};

/// Taxol after `generate_3d_with("med")`. Its GAFF gradient came out NaN.
const TAXOL_SDF: &str = include_str!("data/taxol_med.sdf");

const FORCE_FIELDS: [&str; 5] = ["MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical"];

fn from_smiles(smiles: &str) -> Molecule {
    let mut mol = Molecule::parse(smiles, "smi").expect("parse");
    mol.add_hydrogens();
    assert!(mol.generate_3d(), "generate 3D for {smiles}");
    mol
}

fn coordinates(mol: &Molecule) -> Vec<f64> {
    (0..mol.num_atoms())
        .flat_map(|i| {
            let (x, y, z) = mol.atom(i).unwrap().coords();
            [x, y, z]
        })
        .collect()
}

/// Displace every coordinate by up to `amount` Å, deterministically (a fixed
/// linear congruential sequence), so the geometry sits off its minimum and
/// every term has a force worth checking.
fn kick(mol: &mut Molecule, amount: f64) {
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let kicked: Vec<f64> = coordinates(mol)
        .into_iter()
        .map(|x| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let unit = (state >> 11) as f64 / (1u64 << 53) as f64; // [0, 1)
            x + amount * (2.0 * unit - 1.0)
        })
        .collect();
    assert!(mol.set_coordinates(&kicked));
}

/// The largest disagreement between the analytic gradient and a central
/// difference of the energy, relative to the gradient's size, with the atom it
/// is on. Panics on a gradient or energy that is not finite.
fn worst_mismatch(mol: &Molecule, ff: &str, constraints: &Constraints) -> (f64, usize) {
    let (energy, analytic) = mol
        .energy_gradient(ff, constraints)
        .expect("energy_gradient");
    assert!(energy.is_finite(), "{ff}: energy {energy}");
    for (atom, g) in analytic.iter().enumerate() {
        assert!(
            g.iter().all(|v| v.is_finite()),
            "{ff}: gradient on atom {atom} (element {}) is {g:?}",
            mol.atom(atom as u32).unwrap().atomic_number()
        );
    }

    const H: f64 = 1.0e-5;
    let x0 = coordinates(mol);
    let mut probe = mol.clone();
    let energy_at = |probe: &mut Molecule, x: &[f64]| {
        assert!(probe.set_coordinates(x));
        probe.energy_gradient(ff, constraints).unwrap().0
    };
    let mut worst = (0.0, 0);
    for i in 0..x0.len() {
        let mut x = x0.clone();
        x[i] = x0[i] + H;
        let up = energy_at(&mut probe, &x);
        x[i] = x0[i] - H;
        let down = energy_at(&mut probe, &x);
        let numeric = (up - down) / (2.0 * H);
        let mismatch = (numeric - analytic[i / 3][i % 3]).abs() / numeric.abs().max(1.0);
        assert!(
            mismatch.is_finite(),
            "{ff}: coordinate {i}: numeric {numeric}"
        );
        if mismatch > worst.0 {
            worst = (mismatch, i / 3);
        }
    }
    worst
}

/// Relative agreement a correct analytic gradient reaches against a central
/// difference with h = 1e-5 Å: the difference's own error is O(h²) and
/// rounding, orders of magnitude below this.
const TOLERANCE: f64 = 1.0e-4;

fn assert_gradient_matches(label: &str, mol: &Molecule, ff: &str, constraints: &Constraints) {
    let (mismatch, atom) = worst_mismatch(mol, ff, constraints);
    assert!(
        mismatch < TOLERANCE,
        "{label} / {ff}: analytic gradient is off by {mismatch:.2e} (relative) on atom {atom} \
         (element {})",
        mol.atom(atom as u32).unwrap().atomic_number()
    );
}

/// Every force field, over molecules that between them exercise every term:
/// planar centres (amide, carboxylate, aromatic, carbonyl), pyramidal ones
/// (amine, phosphine), an ion pair, rings, and taxol.
#[test]
fn analytic_gradients_match_the_energy() {
    let mut molecules: Vec<(&str, Molecule)> = [
        ("N-methylacetamide", "CC(=O)NC"),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("trimethylamine", "CN(C)C"),
        ("trimethylphosphine", "CP(C)C"),
        ("ammonium acetate", "[NH4+].CC(=O)[O-]"),
        ("cyclohexene", "C1CCC=CC1"),
    ]
    .into_iter()
    .map(|(name, smiles)| (name, from_smiles(smiles)))
    .collect();
    molecules.push((
        "taxol",
        Molecule::parse(TAXOL_SDF, "sdf").expect("taxol fixture"),
    ));

    for (name, mol) in &mut molecules {
        kick(mol, 0.05);
        for ff in FORCE_FIELDS {
            assert_gradient_matches(name, mol, ff, &Constraints::new());
        }
    }
}

/// GAFF taxol, unkicked: the geometry as generated, where GAFF's gradient was
/// NaN and L-BFGS stopped after one step.
#[test]
fn gaff_gradient_is_finite_on_generated_taxol() {
    let mol = Molecule::parse(TAXOL_SDF, "sdf").expect("taxol fixture");
    assert_gradient_matches("taxol (as generated)", &mol, "GAFF", &Constraints::new());
}

fn kicked_butane() -> Molecule {
    let mut mol = from_smiles("CCCC");
    kick(&mut mol, 0.05);
    mol
}

/// Restraints add their own energy and gradient to the surface, and the two
/// must agree like any force-field term's. The angle restraint's did not: its
/// energy is in degrees, its gradient came out per radian.
#[test]
fn restraint_gradients_match_their_energy() {
    let mol = kicked_butane();
    let mut angle = Constraints::new();
    angle.angle(0, 1, 2, 95.0);
    let mut torsion = Constraints::new();
    torsion.torsion(0, 1, 2, 3, 90.0);
    let mut distance = Constraints::new();
    distance.distance(0, 3, 3.2);

    for (label, constraints) in [
        ("angle", angle),
        ("torsion", torsion),
        ("distance", distance),
    ] {
        for ff in ["UFF", "GAFF", "MMFF94"] {
            assert_gradient_matches(
                &format!("butane, {label} restraint"),
                &mol,
                ff,
                &constraints,
            );
        }
    }
}
