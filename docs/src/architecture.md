# Architecture

This chapter explains how the bindings are put together — the layers, how Rust
talks to OpenBabel's C++, how OpenBabel gets built, and the two cross-cutting
concerns that shape the whole safe API: **thread safety** and **runtime
discovery of OpenBabel's data files**.

## Three layers

```text
   your code
      │
      ▼
┌──────────────┐   safe, idiomatic Rust: Molecule, Atom, Bond, …
│  openbabel   │   every call runs under a global lock (with_ob)
└──────┬───────┘
       │  calls
       ▼
┌──────────────┐   unsafe FFI: a #[cxx::bridge] mod that mirrors the C++ shim
│ openbabel-sys│   build.rs compiles OpenBabel from source and links it
└──────┬───────┘
       │  FFI
       ▼
┌──────────────┐   ob_shim::* C++ functions wrapping OpenBabel's C++ API
│  C++ shim    │   shim/shim.h + shim/shim.cc
└──────┬───────┘
       │
       ▼
   OpenBabel (libopenbabel, built from vendor/openbabel-src)
```

- **`openbabel-sys`** is the low-level plumbing. It is `unsafe` and mirrors the
  shim one-to-one. You normally don't use it directly.
- **`openbabel`** is the crate you use. It turns raw pointers and out-parameters
  into `Molecule`, `Atom`, `Result`, `Option`, iterators, and RAII, and it
  enforces the thread-safety rule described below.

## The cxx bridge and the C++ shim

OpenBabel's public interface is C++ (templates, references, STL types), which
Rust cannot call directly. So `openbabel-sys` defines a small C++ **shim**
(`shim/shim.h` and `shim/shim.cc`) that exposes a flat, C-friendly surface —
free functions in the `ob_shim` namespace that take and return simple types, and
opaque owner types (`Molecule`, `Smarts`, `Transform`, `Constraints`,
`Optimizer`) that wrap the corresponding OpenBabel objects.

The Rust side declares the *same* surface inside a
[`#[cxx::bridge]`](https://cxx.rs) module (`openbabel-sys/src/lib.rs`). `cxx`
checks — **at compile time** — that every signature on the Rust side matches the
C++ header. A mismatch is a build error, not undefined behavior at runtime. The
type mapping follows cxx's conventions:

| Rust                    | C++                 | Used for            |
| ----------------------- | ------------------- | ------------------- |
| `&Molecule`             | `const Molecule&`   | read-only calls     |
| `Pin<&mut Molecule>`    | `Molecule&`         | mutating calls      |
| `&str`                  | `rust::Str`         | passing text in     |
| `String`                | `rust::String`      | returning text out  |
| `UniquePtr<Molecule>`   | `std::unique_ptr`   | ownership transfer  |

Because the shim mirrors the bridge exactly, adding a capability is a matter of:
write the C++ wrapper in `shim.*`, declare it in the bridge, then expose it
safely from `openbabel`.

## Building OpenBabel from source {#building-openbabel-from-source}

`openbabel-sys/build.rs` compiles OpenBabel as part of the crate build, so there
is no runtime dependency on a system OpenBabel install. Its pipeline:

1. **Build + install OpenBabel** from the `vendor/openbabel-src` submodule into
   `OUT_DIR` using the [`cmake`](https://docs.rs/cmake) crate. Eigen (from
   `vendor/eigen`) is pointed at OpenBabel's `find_package(Eigen3)`, which
   enables `HAVE_EIGEN3` — this compiles `OBAlign` (structure superposition) and
   unlocks distance-geometry 3D generation.
2. **Compile the cxx bridge + shim**, pointing the C++ compiler at the
   freshly installed OpenBabel headers.
3. **Link** against the OpenBabel import library.
4. **Make the runtime discoverable**: bake the plugin/data directory paths into
   a generated `paths.rs`, and (on Windows, which has no rpath) copy
   `openbabel-3.dll` and its plugins next to the test/exe binaries.

The first build compiles all of OpenBabel and is slow (~10–20 minutes);
subsequent builds are incremental unless the shim or bridge changes. See
[Building & Platform Notes](./building.md) for the MSVC-specific fixes this
script applies.

## Thread safety: one global lock

**OpenBabel is not thread-safe.** It keeps global mutable state — shared plugin
singletons, aromaticity/ring-perception caches, and so on — so concurrent calls
from multiple threads can corrupt memory.

Because `openbabel` exposes a *safe* API, it must make that corruption
impossible. It does so by serializing **every** entry into OpenBabel behind a
single global lock. The mechanism is a small helper in `openbabel/src/lib.rs`:

```rust,ignore
static OB_LOCK: Mutex<()> = Mutex::new(());

pub(crate) fn with_ob<R>(f: impl FnOnce() -> R) -> R {
    init();
    let _guard = OB_LOCK.lock().unwrap_or_else(|p| p.into_inner());
    f()
}
```

Every FFI call goes through `with_ob`. Two consequences worth internalizing:

- **Calls from multiple threads are correct but do not run concurrently.** For
  throughput, use multiple *processes* rather than threads.
- **The lock is not reentrant.** Inside a `with_ob` closure, keep the raw
  `ffi::…` calls only — never call back into a public method that would itself
  take the lock, or you will deadlock. This is why the internal code passes
  around `as_inner()` / `as_inner_pin_mut()` handles inside the closure.

## Runtime initialization: finding plugins and data {#runtime-initialization-finding-plugins-and-data}

OpenBabel loads its format plugins (`.obf`) from `BABEL_LIBDIR` and reads
element/forcefield data from `BABEL_DATADIR`, both lazily on first use. The safe
API sets these for you, exactly once, via `init()` (called from `with_ob`):

- The paths point at the plugin/data directories **bundled with the linked
  library** (baked in by `build.rs`), so they always match the OpenBabel version
  you built.
- They are set through the **C runtime** (`ffi::set_env`), not `std::env::set_var`.
  On Windows, OpenBabel reads them with `getenv()`, which does not observe
  variables set through the Win32 environment block that `std::env` uses.
- `init()` **overrides** any pre-existing values. A stale system OpenBabel
  install often leaves `BABEL_DATADIR` pointing at its own, version-mismatched
  data directory, which silently breaks data-driven plugins (descriptors, some
  fingerprints). The bundled data must match the bundled library.

You rarely call `init()` yourself — the public API does it — but it is `pub` and
idempotent if you want to force it early.
