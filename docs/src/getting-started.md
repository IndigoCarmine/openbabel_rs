# Getting Started

## Prerequisites

- **Rust** (on Windows, the MSVC toolchain).
- **A C++ compiler and CMake** — `build.rs` uses them to compile OpenBabel.
- **The git submodules checked out** (OpenBabel and Eigen sources):

  ```sh
  git submodule update --init --recursive
  ```

> **Heads-up: the first build is slow.** It compiles all of OpenBabel from
> source and takes roughly **10–20 minutes**. Every build after that is
> incremental and fast. See [Building & Platform Notes](./building.md) for
> details and Windows/MSVC specifics.

## Add the dependency

This crate is not yet published to crates.io, so depend on it via git or a path:

```toml
[dependencies]
openbabel = { git = "https://github.com/IndigoCarmine/openbabel_rs" }
```

## Build, test, run

```sh
cargo build --workspace
cargo test  --workspace
cargo run   -p openbabel-cli -- "c1ccccc1"   # inspect benzene
```

## Your first molecule

```rust,no_run
use openbabel::Molecule;

# fn main() -> Result<(), openbabel::Error> {
// Parse a molecule from SMILES.
let mut mol = Molecule::parse("CCO", "smi")?;   // ethanol

// Query basic properties.
assert_eq!(mol.formula(), "C2H6O");
assert!((mol.molar_mass() - 46.07).abs() < 0.05);

// Mutate: add explicit hydrogens.
mol.add_hydrogens();
assert_eq!(mol.num_atoms(), 9);

// Serialize to another format.
println!("{}", mol.write("can")?.trim());       // canonical SMILES
# Ok(())
# }
```

A few things worth noticing, because they shape the whole API:

- **`Molecule::parse(data, format)`** takes the OpenBabel *format id* as a
  string — `"smi"`, `"mol"`, `"sdf"`, `"pdb"`, `"inchi"`, and so on. The same id
  set is used by `write`. Parsing returns [`Result`], and an unreadable input
  yields [`Error::Parse`].
- **Methods that read** (`formula`, `molar_mass`, `atoms`, …) take `&self`;
  **methods that mutate** the molecule (`add_hydrogens`, `generate_3d`,
  `optimize_geometry`, …) take `&mut self`. This mirrors OpenBabel's own model
  while giving you Rust's borrow-checker guarantees.
- **Fallible operations return `Option` or `Result`.** Many OpenBabel routines
  signal failure by returning a null/empty result rather than throwing; the safe
  API surfaces that as `None`/`Err` so you handle it explicitly.

## The demo CLI

`openbabel-demo` (the `cli` crate) parses a SMILES string and prints a broad
tour of the binding surface — formula, masses, charge, ring counts, logP, TPSA,
a SMARTS hydroxyl count, stereocenters, an optional SVG depiction, canonical
SMILES, InChI/InChIKey, a 3D structure with MMFF94 energy and optimization,
self-alignment RMSD, a conformer search, a Spectrophore descriptor, Gasteiger
partial charges, and a per-atom listing. It is a good, runnable map of what the
crate can do:

```sh
cargo run -p openbabel-cli -- "CC(=O)Oc1ccccc1C(=O)O" aspirin.svg
```
