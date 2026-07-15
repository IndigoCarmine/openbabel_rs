# openbabel-rs

Rust bindings to the [OpenBabel](https://github.com/openbabel/openbabel)
cheminformatics toolkit, built with [`cxx`](https://cxx.rs) over a thin C++ shim.

## Workspace layout

| Crate           | Role                                                                 |
| --------------- | ------------------------------------------------------------------- |
| `openbabel-sys` | Low-level FFI: a `cxx` bridge + C++ shim; `build.rs` builds OpenBabel from source and links it. |
| `openbabel`     | Safe, idiomatic API (`Molecule`, `Atom`, `Bond`, `Error`).          |
| `cli`           | `openbabel-demo` — a small SMILES-inspection demo.                   |

Dependencies are vendored as git submodules: OpenBabel at `vendor/openbabel-src`
(tag `openbabel-3-2-1`) and Eigen (header-only, for structure alignment) at
`vendor/eigen` (tag `3.4.0`). Fetch both with
`git submodule update --init --recursive`.

## Current scope

Core (`Molecule`, `Atom`, `Bond`):

- Read/write molecules in any OpenBabel format (`Molecule::parse` / `write`) —
  SMILES, MOL, SDF, PDB, canonical SMILES, …
- Properties: formula, molar mass, exact mass, total charge, atom/bond counts,
  title.
- Atom access: element, coordinates, formal charge, aromaticity, ring membership.
- Bond access: begin/end atoms, order, aromaticity, ring membership.
- `add_hydrogens` / `remove_hydrogens`.

Analysis:

- SMARTS substructure matching (`SmartsPattern`): `matches`, `num_matches`,
  `match_indices`.
- Fingerprints & Tanimoto similarity (`Fingerprint`): `FP2`, `FP3`, `FP4`,
  `MACCS`, …
- Numeric descriptors (`Molecule::descriptor`, plus `logp` / `tpsa` /
  `molar_refractivity`).

3D structures & force fields:

- 3D coordinate generation (`Molecule::generate_3d` / `generate_3d_with`,
  `dimension`, `has_3d`).
- Force-field single-point energy and geometry optimization
  (`Molecule::energy` / `optimize_geometry`) with `MMFF94`, `MMFF94s`, `UFF`,
  `GAFF`, `Ghemical`; unit via `forcefield_energy_unit`.

> `generate_3d` uses the fast fragment-builder path. Distance-geometry
> generation, which needs Eigen, is also available now that Eigen is enabled.

Charges, atom chemistry & identifiers:

- Partial atomic charges (`Molecule::compute_charges`) via the `gasteiger`,
  `mmff94`, `eem`, `eqeq`, `qeq`, and `qtpie` models; per-atom values through
  `Atom::partial_charge`.
- Richer atom chemistry: `degree`, `total_valence`, `implicit_hydrogens`,
  `hybridization`, `is_hbond_donor`, `is_hbond_acceptor`.
- InChI / InChIKey identifiers (`Molecule::inchi` / `inchikey`), backed by the
  InChI library bundled with OpenBabel.

Structure alignment:

- Least-squares superposition of one molecule onto another
  (`Molecule::align_to` / `align_to_with`, Kabsch algorithm): returns the RMSD
  and updates the aligned molecule's coordinates in place. Backed by OpenBabel's
  `OBAlign`, enabled by the vendored Eigen submodule.

## Thread safety

OpenBabel is not thread-safe — it keeps global mutable state (shared plugin
singletons, aromaticity/ring perception caches). This crate therefore
serializes every call into OpenBabel behind a global lock, so the safe API
cannot be used to trigger data races. Calls from multiple threads are correct
but do not run concurrently; for throughput, use multiple processes.

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
- Building the bundled InChI from source on MSVC needs two fixes, both applied
  by `build.rs`: OpenBabel force-sets `OPENBABEL_USE_SYSTEM_INCHI=ON` when
  `OB_USE_PREBUILT_BINARIES` is on (its MSVC default), which makes configuration
  demand a system InChI we don't have — so we turn that flag off; and the
  vendored InChI omits four AuxInfo functions its ABI wrappers reference, so
  `build.rs` drops inert stubs into the InChI tree so `inchi.dll` links. Both
  are runtime-inert (OpenBabel never calls those functions for InChI output).
  `inchi.dll` is copied alongside the other runtime DLLs.

## License

Bindings © their authors. OpenBabel is GPL-2.0; code that links it inherits
that license. See `vendor/openbabel-src/COPYING`.
