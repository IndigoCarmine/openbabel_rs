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

## Prebuilt OpenBabel for CI {#prebuilt-openbabel}

Every tagged release has OpenBabel, already built, attached for each platform
the release was tested on:

| Archive                                         | Built on         |
| ----------------------------------------------- | ---------------- |
| `openbabel-sys-x86_64-unknown-linux-gnu.tar.gz` | `ubuntu-latest`  |
| `openbabel-sys-x86_64-pc-windows-msvc.tar.gz`   | `windows-latest` |
| `openbabel-sys-aarch64-apple-darwin.tar.gz`     | `macos-latest`   |

Point `OPENBABEL_SYS_PREBUILT_DIR` at an unpacked archive and `build.rs` uses it
instead of compiling OpenBabel — no CMake, no Perl, no 10–20 minutes. The cxx
bridge and the shim are still compiled, which is quick.

**The archive must come from the same openbabel-sys.** The shim is compiled
against the archive's headers and links into its library, so an archive from
another version fails to link or, worse, misbehaves. Each archive records what
it was built from in `openbabel-sys-prebuilt.txt` — the target, the OpenBabel
and Eigen submodule commits, and a hash of `build.rs` — and `build.rs` fails the
build on any difference instead of quietly falling back to a source build.
Reading the commits needs `git` on `PATH`.

In practice: depend on openbabel_rs by tag, and download the release with that
tag. In the application's workflow, after installing Rust and before any cargo
command:

```yaml
- name: Fetch prebuilt OpenBabel
  shell: bash
  env:
    GH_TOKEN: ${{ github.token }}
  run: |
    tag=$(sed -n 's/.*openbabel_rs[^?]*?tag=\([^#]*\)#.*/\1/p' Cargo.lock | head -n 1)
    triple=$(rustc -vV | sed -n 's/^host: //p')
    archive="openbabel-sys-$triple.tar.gz"
    if [ -z "$tag" ]; then
      echo "::warning::openbabel_rs is not pinned to a tag; building OpenBabel from source"
      exit 0
    fi
    cd "$RUNNER_TEMP"
    if ! gh release download "$tag" -R IndigoCarmine/openbabel_rs -p "$archive"; then
      echo "::warning::openbabel_rs $tag has no $archive; building OpenBabel from source"
      exit 0
    fi
    mkdir -p openbabel-sys
    tar -xzf "$archive" -C openbabel-sys
    echo "OPENBABEL_SYS_PREBUILT_DIR=$RUNNER_TEMP/openbabel-sys" >> "$GITHUB_ENV"
```

The tag is read from `Cargo.lock`, so it cannot drift from the version cargo
actually builds. A few things to know:

- The archives are built on GitHub's `*-latest` runners and meant for them. An
  older Linux distribution (older glibc) or macOS release may refuse to load
  them.
- When cross-compiling, use the target triple rather than the host's.
- cargo still clones the submodules of a git dependency; caching `~/.cargo/git`
  (for example with `Swatinem/rust-cache`) saves that download too.

The archives come from `.github/workflows/tests.yml`: on a tag push each
platform builds OpenBabel with `OPENBABEL_SYS_EXPORT_DIR` set, packs the export,
throws the source build away and runs the tests again against the unpacked
archive. Only when every platform passes are the archives attached to the
release — so tag every version bump. A manual run of the workflow does all of
it short of publishing.

## Linux / macOS: a relocatable install

OpenBabel normally bakes its install prefix into the libraries it installs — as
the rpath on Linux and the install name on macOS — so the tree only loads where
it was built. `build.rs` turns on OpenBabel's wheel-build mode (`BUILD_BY_PIP`),
which uses loader-relative paths (`$ORIGIN`, `@rpath`) instead; that is what lets
a prebuilt archive be unpacked anywhere. With the Python bindings off, it also
leaves out Cairo (PNG depiction), which would otherwise tie the tree to the
build machine's Cairo.

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
