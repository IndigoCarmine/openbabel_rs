//! Safe, idiomatic Rust bindings to the [OpenBabel](https://openbabel.org)
//! cheminformatics toolkit.
//!
//! This crate wraps the unsafe [`openbabel_sys`] FFI layer (a cxx bridge over a
//! C++ shim) in a safe API. It currently covers a core MVP surface — molecule
//! I/O, basic properties, and atom/bond access — and will grow over time.
//!
//! ```no_run
//! let mut mol = openbabel::Molecule::parse("CCO", "smi").unwrap();
//! assert_eq!(mol.formula(), "C2H6O");
//! mol.add_hydrogens();
//! println!("{}", mol.write("can").unwrap()); // canonical SMILES
//! ```

mod atom;
mod bond;
mod error;
mod mol;

pub use atom::Atom;
pub use bond::Bond;
pub use error::Error;
pub use mol::Molecule;

use std::sync::Once;

static INIT: Once = Once::new();

/// Point OpenBabel at the format plugins and data files installed alongside the
/// linked library.
///
/// OpenBabel loads its format plugins (`.obf`) from `BABEL_LIBDIR` and reads
/// element/forcefield data from `BABEL_DATADIR`, both lazily on first use. We
/// set them (unless already set in the environment) before any conversion
/// happens. This runs at most once and is safe to call repeatedly; the public
/// API calls it for you.
pub fn init() {
    INIT.call_once(|| {
        // Set via the C runtime (not std::env::set_var): on Windows OpenBabel
        // reads these with getenv(), which does not observe variables set
        // through the Win32 environment block that std::env uses.
        if std::env::var_os("BABEL_LIBDIR").is_none() {
            openbabel_sys::ffi::set_env("BABEL_LIBDIR", openbabel_sys::paths::BABEL_LIBDIR);
        }
        if std::env::var_os("BABEL_DATADIR").is_none() {
            openbabel_sys::ffi::set_env("BABEL_DATADIR", openbabel_sys::paths::BABEL_DATADIR);
        }
    });
}

/// OpenBabel release version string, e.g. `"3.2.1"`.
pub fn version() -> String {
    openbabel_sys::ffi::release_version()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_reported() {
        let v = version();
        assert!(!v.is_empty(), "version string should not be empty");
        assert!(v.starts_with('3'), "unexpected OpenBabel version: {v:?}");
    }
}
