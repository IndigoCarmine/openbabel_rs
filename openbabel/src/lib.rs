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
//!
//! A prose guide (architecture, feature tour, key concepts) — available in both
//! English and Japanese — lives at
//! <https://indigocarmine.github.io/openbabel_rs/>.
#![doc(html_root_url = "https://indigocarmine.github.io/openbabel_rs/api/")]
#![deny(missing_docs)]

mod atom;
mod bond;
mod constraints;
pub mod elements;
mod error;
mod ff;
mod fingerprint;
mod minimize;
mod mol;
mod reaction;
mod residue;
mod ring;
mod smarts;
mod transform;
mod unitcell;

pub use atom::{Atom, AtomMut, Winding};
pub use bond::{Bond, BondMut};
pub use constraints::{Axis, Constraints};
pub use error::Error;
pub use fingerprint::Fingerprint;
pub use minimize::{Algorithm, Minimizer, OptStep, Optimization, StopReason};
pub use mol::{Molecule, SvgOptions};
pub use reaction::Reaction;
pub use residue::Residue;
pub use ring::Ring;
pub use smarts::SmartsPattern;
pub use transform::Transform;
pub use unitcell::{LatticeType, UnitCell};

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

/// The plugin and data directories to hand OpenBabel, as
/// `(BABEL_LIBDIR, BABEL_DATADIR)`.
///
/// A shipped application does not have the build machine's directories, so an
/// application-relative layout wins when it is there: the `.obf` plugins next to
/// the executable and the data directory at `<exe_dir>/data`. That mirrors what
/// `openbabel-sys`'s build script already lays out on Windows, where OpenBabel
/// finds plugins beside `openbabel-3.dll` regardless of `BABEL_LIBDIR`.
///
/// Both directories are decided together, from one signal, so a bundled data
/// directory can never end up paired with the build tree's library or the other
/// way round — data and library have to match.
fn runtime_dirs(exe_dir: Option<&std::path::Path>) -> (String, String) {
    let bundled = exe_dir.and_then(|dir| {
        // The plugins sit beside the executable, but they are not the signal: a
        // cargo target directory has them there too and yet wants the baked
        // paths. Only a packaged application has the data directory beside it.
        let data = dir.join("data");
        if !data.is_dir() {
            return None;
        }
        Some((dir.to_str()?.to_owned(), data.to_str()?.to_owned()))
    });

    bundled.unwrap_or_else(|| {
        (
            openbabel_sys::paths::BABEL_LIBDIR.to_owned(),
            openbabel_sys::paths::BABEL_DATADIR.to_owned(),
        )
    })
}

/// Point OpenBabel at its format plugins and data files.
///
/// OpenBabel loads its format plugins (`.obf`) from `BABEL_LIBDIR` and reads
/// element/forcefield data from `BABEL_DATADIR`, both lazily on first use, so
/// they have to be set before any conversion happens. This runs at most once and
/// is safe to call repeatedly; the public API calls it for you.
///
/// See [`runtime_dirs`] for where they end up pointing. To ship an application,
/// put the `.obf` plugins (and, on Windows, `openbabel-3.dll`) next to the
/// executable and the data directory at `<exe_dir>/data` — otherwise the binary
/// depends on directories that only exist on the machine that built it.
pub fn init() {
    INIT.call_once(|| {
        let exe = std::env::current_exe().ok();
        let (libdir, datadir) = runtime_dirs(exe.as_deref().and_then(|p| p.parent()));

        // Set via the C runtime (not std::env::set_var): on Windows OpenBabel
        // reads these with getenv(), which does not observe variables set
        // through the Win32 environment block that std::env uses.
        //
        // We *override* any pre-existing values: a stale system OpenBabel
        // install (e.g. an older release) commonly leaves BABEL_DATADIR
        // pointing at its own, version-mismatched data directory, which
        // silently breaks data-driven plugins (descriptors, some fingerprints).
        // Our data must match our library.
        openbabel_sys::ffi::set_env("BABEL_LIBDIR", &libdir);
        openbabel_sys::ffi::set_env("BABEL_DATADIR", &datadir);
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

/// Serialize several molecules into one multi-record document in `format`.
///
/// Each molecule is written with [`Molecule::write`] and the records are
/// concatenated — for `"sdf"` this yields a valid multi-molecule SDF, for
/// `"smi"` one SMILES per line, and so on. The inverse of
/// [`Molecule::parse_many`]. Returns [`Error::UnknownFormat`] if `format` is
/// unknown.
pub fn write_many(molecules: &[Molecule], format: &str) -> Result<String, Error> {
    let mut out = String::new();
    for mol in molecules {
        out.push_str(&mol.write(format)?);
    }
    Ok(out)
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

    /// A packaged application: `data` sits beside the executable, so neither
    /// directory may point back at the machine that built it.
    #[test]
    fn a_bundled_layout_wins_over_the_baked_paths() {
        let tmp = std::env::temp_dir().join("openbabel_rs_bundled_layout");
        let data = tmp.join("data");
        std::fs::create_dir_all(&data).expect("create bundle");

        let (libdir, datadir) = runtime_dirs(Some(&tmp));
        assert_eq!(libdir, tmp.to_str().unwrap());
        assert_eq!(datadir, data.to_str().unwrap());

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// A cargo target directory has the plugins next to the test binary but no
    /// `data`, so it has to keep using the build-time paths.
    #[test]
    fn a_build_tree_falls_back_to_the_baked_paths() {
        let tmp = std::env::temp_dir().join("openbabel_rs_build_tree");
        std::fs::create_dir_all(&tmp).expect("create dir");
        std::fs::write(tmp.join("formats_common.obf"), b"").expect("plugin");

        let (libdir, datadir) = runtime_dirs(Some(&tmp));
        assert_eq!(libdir, openbabel_sys::paths::BABEL_LIBDIR);
        assert_eq!(datadir, openbabel_sys::paths::BABEL_DATADIR);

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn an_unknown_exe_location_falls_back_to_the_baked_paths() {
        let (libdir, datadir) = runtime_dirs(None);
        assert_eq!(libdir, openbabel_sys::paths::BABEL_LIBDIR);
        assert_eq!(datadir, openbabel_sys::paths::BABEL_DATADIR);
    }

    /// The bug this resolution exists for: a data directory only the build
    /// machine has means a shipped binary silently loses every data-driven
    /// plugin — force fields included — while still starting up fine. A force
    /// field's energy unit is read from that data, so it stands in for "the data
    /// directory actually resolved".
    #[test]
    fn forcefield_data_resolves_after_init() {
        init();
        assert_eq!(
            forcefield_energy_unit("MMFF94").as_deref(),
            Some("kcal/mol"),
            "force-field data did not load — BABEL_DATADIR is wrong"
        );
    }
}
