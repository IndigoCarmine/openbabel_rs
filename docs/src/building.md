# Building & Platform Notes

## Requirements

- **Rust** (MSVC toolchain on Windows).
- **A C++ compiler and CMake** — `build.rs` uses them to compile OpenBabel.
- **Submodules checked out**: `git submodule update --init --recursive`.

```sh
cargo build --workspace
cargo test  --workspace
cargo run   -p openbabel-cli -- "c1ccccc1"   # benzene
```

> The **first** build compiles all of OpenBabel from source and takes ~10–20
> minutes. Subsequent builds are incremental and fast.

## Why building can be slow the first time

`openbabel-sys/build.rs` compiles OpenBabel from the `vendor/openbabel-src`
submodule (via CMake) and links it, rather than relying on a system install.
This keeps the binding self-contained and version-matched, at the cost of a long
initial build. CMake caches its work, so later `cargo build`s are no-ops for the
C++ side unless the shim or the cxx bridge changes. See
[Architecture → Building OpenBabel from source](./architecture.md#building-openbabel-from-source).

## Windows / MSVC specifics

`build.rs` handles three Windows/MSVC details that are easy to trip over:

- **C++ flags.** The `cmake` crate replaces `CMAKE_CXX_FLAGS`, dropping the
  `/DWIN32 /D_WINDOWS /EHsc /GR` that OpenBabel relies on (for example, its
  `strcasecmp` shim is guarded by `#if defined(WIN32)`); `build.rs` adds them
  back.
- **Plugin & DLL discovery.** OpenBabel discovers its format plugins (`.obf`) in
  the directory of `openbabel-3.dll` (not via `BABEL_LIBDIR` on Windows), so
  `build.rs` copies the DLL **and** the plugins next to the test/exe binaries.
  Data files are located via `BABEL_DATADIR`, which the safe wrapper sets
  through the C runtime so OpenBabel's `getenv` observes it (see
  [Architecture → Runtime initialization](./architecture.md#runtime-initialization-finding-plugins-and-data)).
- **Bundled InChI on MSVC.** Building the bundled InChI from source needs two
  fixes, both applied by `build.rs`: OpenBabel force-sets
  `OPENBABEL_USE_SYSTEM_INCHI=ON` when `OB_USE_PREBUILT_BINARIES` is on (its MSVC
  default), which would demand a system InChI we don't have — so `build.rs`
  turns that off; and the vendored InChI omits four AuxInfo functions its ABI
  wrappers reference, so `build.rs` drops inert stubs into the InChI tree so
  `inchi.dll` links. Both are runtime-inert (OpenBabel never calls those
  functions for InChI output). `inchi.dll` is copied alongside the other runtime
  DLLs.

## Runtime data files

At runtime, OpenBabel needs its plugin (`.obf`) and data directories. `build.rs`
bakes their paths into a generated `paths.rs`, and the safe API's `init()`
points `BABEL_LIBDIR` / `BABEL_DATADIR` at the bundled copies through the C
runtime. You do not need to set any environment variables yourself.
