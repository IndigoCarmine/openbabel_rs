//! Shared helpers for the ported OpenBabel test suite (the `ob_*` integration
//! tests).
//!
//! These tests are faithful Rust ports of OpenBabel's own C++ tests under
//! `vendor/openbabel-src/test/`. They load the *exact same* input and expected
//! result files that OpenBabel's `test_runner` uses, resolved from the vendored
//! submodule (OpenBabel's `TESTDATADIR`, i.e. `test/files/`).
//!
//! This module is included with `mod common;` from each `ob_*.rs` test file; it
//! is not itself a test binary (Cargo only treats top-level `tests/*.rs` files
//! as targets, not files inside subdirectories).

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Absolute path to OpenBabel's `test/files/` directory in the vendored
/// submodule — the analogue of the C++ suite's `TESTDATADIR`.
///
/// The `openbabel` crate lives at `<workspace>/openbabel`; the submodule is at
/// `<workspace>/vendor/openbabel-src`. The submodule must be checked out
/// (`git submodule update --init`) — which it must be anyway for the crate to
/// build, since `openbabel-sys` compiles OpenBabel from this same source.
pub fn test_data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("vendor")
        .join("openbabel-src")
        .join("test")
        .join("files")
}

/// Resolve `name` (e.g. `"attype.00.smi"`) inside the OpenBabel test-data
/// directory. Mirrors `OBTestUtil::GetFilename` / the C++ `TESTDATADIR + name`.
pub fn ob_test_file(name: &str) -> PathBuf {
    test_data_dir().join(name)
}

/// Read an OpenBabel test-data file to a string, panicking with a clear message
/// (including the resolved path) if it cannot be read — this usually means the
/// git submodule has not been checked out.
pub fn read_test_file(name: &str) -> String {
    let path = ob_test_file(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "cannot read OpenBabel test data file {}: {e}\n\
             (did you run `git submodule update --init`?)",
            path.display()
        )
    })
}

/// Split a line into whitespace-separated tokens, matching OpenBabel's default
/// `tokenize()` (which splits on space, tab, carriage return and newline).
///
/// Empty tokens are dropped, so a blank or all-whitespace line yields an empty
/// vector — exactly like the C++ reference-file readers.
pub fn tokenize(line: &str) -> Vec<&str> {
    line.split([' ', '\t', '\r', '\n'])
        .filter(|t| !t.is_empty())
        .collect()
}
