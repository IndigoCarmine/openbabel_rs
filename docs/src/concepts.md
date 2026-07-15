# Key Concepts

A closer look at the parts of the API whose behavior is subtle enough to be
worth explaining in prose. Each section links to the types involved; the
[API Reference](./api-reference.md) has the exact signatures.

## SMARTS-to-SMARTS transforms {#smarts-to-smarts-transforms}

`Transform` compiles a **SMARTS→SMARTS edit** (OpenBabel's `OBChemTsfm`, the
engine behind SMIRKS-like reaction rules) and applies it to a molecule in place.

```rust,no_run
use openbabel::{Molecule, Transform};

# fn main() -> Result<(), openbabel::Error> {
// Deprotonate a carboxylic acid: O–H on a C(=O)–O becomes O⁻.
let tsfm = Transform::new("[C:1](=O)[OX2:2][H]", "[C:1](=O)[O-:2]")
    .expect("valid transform");
let mut mol = Molecule::parse("CC(=O)O", "smi")?;   // acetic acid
let n = tsfm.apply(&mut mol);                        // number of matches rewritten
# let _ = n;
# Ok(())
# }
```

The key mechanism to understand is the **class labels** (`:1`, `:2`, …). They
tie an atom in the reactant pattern to the same atom in the product pattern, so
the transform knows *which* atoms persist versus which are added or removed. A
transform can change formal charges and bond orders, and delete atoms (for
example, dropping the hydrogen when deprotonating). `apply` returns how many
matches it rewrote, and rewrites **all** of them.

This is the single densest concept in the crate — if a transform does nothing,
the usual cause is a reactant pattern that doesn't match, or class labels that
don't line up between the two sides.

## Constrained minimization {#constrained-minimization}

Force-field minimization is exposed in two layers.

The simple path is `Molecule::optimize_geometry(forcefield, steps)`, which
minimizes in one shot and returns the final energy.

For control, build a `Minimizer` — it bundles the force field, an `Algorithm`
(steepest descent / conjugate gradients / L-BFGS), a step budget, an
energy-convergence threshold, and an optional `Constraints` set:

```rust,no_run
use openbabel::{Algorithm, Minimizer, Molecule};

# fn main() {
let mut mol = Molecule::parse("CCO", "smi").unwrap();
mol.generate_3d();

let mut cfg = Minimizer::new("MMFF94");
cfg.algorithm(Algorithm::ConjugateGradients)
   .max_steps(500)
   .energy_convergence(1e-6)
   .steps_per_frame(10);

// One-shot: returns the final energy.
let _final_energy = mol.optimize_geometry_with(&cfg);
# }
```

Two things to keep in mind:

- **Convergence is energy-only.** OpenBabel 3.2.1 lets you set the *energy*
  criterion (`energy_convergence`); it combines that with a *fixed* internal
  gradient criterion for which there is no public setter. So this API does not
  expose a gradient tolerance — by design, not by omission.
- **Builder methods take `&mut self` and return `&mut Self`,** so they chain and
  can be applied conditionally.

### The `Constraints` builder

A `Constraints` set restrains a minimization. All atom indices are **0-based**
(matching `Atom::index`). You can hold atoms fixed, pin a single Cartesian axis,
restrain distances/angles/torsions to target values, exclude atoms entirely, and
tune the restraint stiffness:

```rust,no_run
use openbabel::Constraints;

let mut c = Constraints::new();
c.fix_atom(0)              // hold atom 0 in place
 .distance(1, 2, 1.54)     // pin the 1–2 bond to 1.54 Å
 .force_factor(50000.0);   // stiffen the restraints
```

Pass the set to a `Minimizer` with `Minimizer::constraints`.

## The pull-based optimization trajectory

`Molecule::minimize(&cfg)` does **not** minimize immediately. It returns an
`Optimization` — a lazy `Iterator` that advances the minimization by
`steps_per_frame` steps on each `next()` and yields an `OptStep` (cumulative
step count, energy, and coordinates) per frame. Iteration ends at convergence or
the step budget, leaving the molecule holding the final geometry.

```rust,no_run
use openbabel::{Minimizer, Molecule};

# fn main() {
let mut mol = Molecule::parse("CCO", "smi").unwrap();
mol.generate_3d();
let cfg = Minimizer::new("MMFF94");

for step in mol.minimize(&cfg) {
    println!("step {:>4}  E = {:.4}", step.step, step.energy);
    // log it, collect it, or `break` to stop early.
}
# }
```

Because the trajectory drives the parent molecule's coordinates, the iterator
**borrows the molecule mutably**. That means: read each frame's geometry from the
yielded `OptStep`, not from the molecule, while iterating — the molecule is
borrowed for the duration of the loop.

## Structure alignment {#structure-alignment}

`Molecule::align_to(&reference)` least-squares-superposes the molecule onto a
reference (the Kabsch algorithm, via OpenBabel's `OBAlign`), returns the RMSD,
and updates the aligned molecule's coordinates in place. `align_to_with` exposes
options (such as whether to include hydrogens and whether to use symmetry).

The important precondition: the default alignment assumes the **two molecules
have the same atom order**. It is ideal for comparing conformers or poses of the
same molecule (identical structures align to ~0 RMSD), not for matching two
arbitrarily-ordered molecules.

## Conformer search

`Molecule::generate_conformers(count)` runs a genetic-algorithm conformer search
and returns how many conformers it produced. It **requires an existing 3D
structure** (call `generate_3d` first). Walk the results with `num_conformers`
and `set_conformer`, and score each with `energy`. Like alignment and
distance-geometry generation, it relies on the vendored Eigen support.
