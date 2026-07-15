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
mod fingerprint;
mod mol;
mod smarts;

pub use atom::Atom;
pub use bond::Bond;
pub use error::Error;
pub use fingerprint::Fingerprint;
pub use mol::{Molecule, SvgOptions};
pub use smarts::SmartsPattern;

use std::sync::{Mutex, Once};

static INIT: Once = Once::new();

/// Serializes all access to OpenBabel.
///
/// OpenBabel is not thread-safe: it keeps global mutable state (shared plugin
/// singletons, aromaticity/ring perception caches, …), so concurrent calls
/// corrupt memory. Because this crate exposes a *safe* API, it must make that
/// impossible — every entry point runs under this lock via [`with_ob`].
static OB_LOCK: Mutex<()> = Mutex::new(());

/// Run `f` with OpenBabel initialized and the global lock held.
///
/// All FFI calls into OpenBabel go through here. The closure must not call
/// another `with_ob` (the lock is not reentrant); keep the raw `ffi::…` calls
/// inside it rather than calling back into public methods.
pub(crate) fn with_ob<R>(f: impl FnOnce() -> R) -> R {
    init();
    let _guard = OB_LOCK.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    f()
}

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
        //
        // We *override* any pre-existing values: a stale system OpenBabel
        // install (e.g. an older release) commonly leaves BABEL_DATADIR
        // pointing at its own, version-mismatched data directory, which
        // silently breaks data-driven plugins (descriptors, some fingerprints).
        // Our bundled data must match our bundled library.
        openbabel_sys::ffi::set_env("BABEL_LIBDIR", openbabel_sys::paths::BABEL_LIBDIR);
        openbabel_sys::ffi::set_env("BABEL_DATADIR", openbabel_sys::paths::BABEL_DATADIR);
    });
}

/// OpenBabel release version string, e.g. `"3.2.1"`.
pub fn version() -> String {
    openbabel_sys::ffi::release_version()
}

/// The energy unit reported by a force field (e.g. `"kcal/mol"`), or `None` if
/// OpenBabel has no such force field.
pub fn forcefield_energy_unit(forcefield: &str) -> Option<String> {
    let unit = with_ob(|| openbabel_sys::ffi::forcefield_unit(forcefield));
    if unit.is_empty() {
        None
    } else {
        Some(unit)
    }
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
