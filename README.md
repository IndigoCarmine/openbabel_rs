# openbabel-rs

Rust bindings to the [OpenBabel](https://github.com/openbabel/openbabel)
cheminformatics toolkit, built with [`cxx`](https://cxx.rs) over a thin C++ shim.

## Workspace layout

| Crate           | Role                                                                 |
| --------------- | ------------------------------------------------------------------- |
| `openbabel-sys` | Low-level FFI: a `cxx` bridge + C++ shim; `build.rs` builds OpenBabel from source and links it. |
| `openbabel`     | Safe, idiomatic API (`Molecule`, `Atom`, `Bond`, `Error`).          |
| `cli`           | `openbabel-demo` — a small SMILES-inspection demo.                   |

OpenBabel itself is vendored as a git submodule at `vendor/openbabel-src`
(pinned to tag `openbabel-3-2-1`).

## Current scope (MVP)

- Read/write molecules in any OpenBabel format (`Molecule::parse` / `write`) —
  SMILES, MOL, SDF, PDB, canonical SMILES, …
- Properties: formula, molar mass, exact mass, total charge, atom/bond counts,
  title.
- Atom access: element, coordinates, formal charge, aromaticity, ring membership.
- Bond access: begin/end atoms, order, aromaticity, ring membership.
- `add_hydrogens` / `remove_hydrogens`.

## Building

Requirements:

- Rust (MSVC toolchain on Windows).
- A C++ compiler and CMake (used by `build.rs` to compile OpenBabel).
- The submodule checked out: `git submodule update --init --recursive`.

```sh
cargo build --workspace
cargo test  --workspace
cargo run   -p openbabel-cli -- "c1ccccc1"   # benzene
```

> The **first** build compiles all of OpenBabel from source and takes ~10–20
> minutes. Subsequent builds are incremental and fast.

### Example

```rust
use openbabel::Molecule;

let mut mol = Molecule::parse("CCO", "smi")?;      // ethanol
assert_eq!(mol.formula(), "C2H6O");
assert!((mol.molar_mass() - 46.07).abs() < 0.05);
mol.add_hydrogens();
assert_eq!(mol.num_atoms(), 9);
println!("{}", mol.write("can")?.trim());          // canonical SMILES
```

## Platform notes

`build.rs` handles two Windows/MSVC specifics that are easy to trip over:

- The `cmake` crate replaces `CMAKE_CXX_FLAGS`, dropping the `/DWIN32
  /D_WINDOWS /EHsc /GR` that OpenBabel relies on (e.g. its `strcasecmp` shim is
  guarded by `#if defined(WIN32)`); we add them back.
- OpenBabel discovers its format plugins (`.obf`) in the directory of
  `openbabel-3.dll` (not via `BABEL_LIBDIR` on Windows), so `build.rs` copies
  the DLL **and** the plugins next to the test/exe binaries. Data files are
  located via `BABEL_DATADIR`, which the safe wrapper sets through the C
  runtime so OpenBabel's `getenv` observes it.

## License

Bindings © their authors. OpenBabel is GPL-2.0; code that links it inherits
that license. See `vendor/openbabel-src/COPYING`.
