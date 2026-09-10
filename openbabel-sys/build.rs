//! Build script for `openbabel-sys`.
//!
//! Pipeline:
//!   1. Build + install OpenBabel (from the `vendor/openbabel-src` submodule)
//!      into `OUT_DIR` via the `cmake` crate -- or, with
//!      `OPENBABEL_SYS_PREBUILT_DIR` set, take an install tree built elsewhere
//!      (see `use_prebuilt`).
//!   2. Compile the cxx bridge (`src/lib.rs`) + C++ shim (`shim/shim.cc`),
//!      pointing the C++ compiler at the installed OpenBabel headers.
//!   3. Link against the OpenBabel import library.
//!   4. Make the runtime discoverable: bake the plugin/data directories into a
//!      generated `paths.rs`, and copy `openbabel-3.dll` next to the eventual
//!      test/exe binaries (Windows has no rpath).
//!
//! The first build compiles all of OpenBabel and is slow (~10-20 min); later
//! builds are incremental (cmake no-ops) unless the shim or bridge changes.
//! A prebuilt skips that compile altogether, which is what CI wants.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Names an OpenBabel install tree built elsewhere, to use instead of compiling
/// one. See [`use_prebuilt`].
const PREBUILT_DIR_ENV: &str = "OPENBABEL_SYS_PREBUILT_DIR";

/// Names a directory to copy the install tree into, laid out as a prebuilt.
/// See [`export_prebuilt`].
const EXPORT_DIR_ENV: &str = "OPENBABEL_SYS_EXPORT_DIR";

/// The file at the root of an install tree recording what it was built from.
/// See [`prebuilt_key`].
const KEY_FILE: &str = "openbabel-sys-prebuilt.txt";

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let workspace_root = manifest_dir.parent().unwrap();
    let ob_src = workspace_root.join("vendor").join("openbabel-src");
    // Eigen (header-only) is vendored as a submodule. Pointing OpenBabel's
    // `find_package(Eigen3)` at it enables HAVE_EIGEN3, which compiles OBAlign
    // (structure superposition) and unlocks distance-geometry 3D generation.
    let eigen_dir = workspace_root.join("vendor").join("eigen");

    println!("cargo:rerun-if-changed=shim/shim.cc");
    println!("cargo:rerun-if-changed=shim/shim.h");
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed={PREBUILT_DIR_ENV}");
    println!("cargo:rerun-if-env-changed={EXPORT_DIR_ENV}");

    assert!(
        ob_src.join("CMakeLists.txt").exists(),
        "OpenBabel source not found at {}.\n\
         Run: git submodule update --init --recursive",
        ob_src.display()
    );
    assert!(
        eigen_dir.join("Eigen").join("Core").exists(),
        "Eigen headers not found at {}.\n\
         Run: git submodule update --init --recursive",
        eigen_dir.display()
    );

    let key = prebuilt_key(&manifest_dir, &ob_src, &eigen_dir);
    let prebuilt_dir = env::var_os(PREBUILT_DIR_ENV).filter(|dir| !dir.is_empty());

    // 1. Build + install OpenBabel into OUT_DIR, or take a prebuilt install.
    let dst = match &prebuilt_dir {
        Some(dir) => use_prebuilt(PathBuf::from(dir), key.as_deref()),
        None => build_openbabel(&ob_src, &eigen_dir, key.as_deref()),
    };

    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let is_windows = target_os == "windows";

    // Installed layout (Windows):
    //   <dst>/include/openbabel3/openbabel/*.h   headers
    //   <dst>/bin/openbabel-3.lib                import library
    //   <dst>/bin/openbabel-3.dll                runtime library
    //   <dst>/bin/*.obf                          format/plugin modules (BABEL_LIBDIR)
    //   <dst>/bin/data/                          runtime data (BABEL_DATADIR)
    //
    // Installed layout (Unix / macOS):
    //   <dst>/include/openbabel3/openbabel/*.h        headers
    //   <dst>/lib/libopenbabel.dylib|.so              runtime + link library
    //   <dst>/lib/openbabel/<version>/*.so            plugin modules (BABEL_LIBDIR)
    //   <dst>/share/openbabel/<version>/              runtime data (BABEL_DATADIR)
    //
    // Both, added by this script:
    //   <dst>/include/openbabel-sys/forcefields/*.h   patched private headers
    //   <dst>/openbabel-sys-prebuilt.txt              what the tree was built from
    let include_dir = dst.join("include").join("openbabel3");

    assert!(
        include_dir.join("openbabel").join("mol.h").exists(),
        "expected OpenBabel headers at {}",
        include_dir.display()
    );

    // Resolve the platform-specific link directory + library name, and the
    // runtime plugin (BABEL_LIBDIR) / data (BABEL_DATADIR) directories.
    let (link_search_dir, link_lib_name, babel_libdir, babel_datadir) = if is_windows {
        let bin_dir = dst.join("bin");
        (
            bin_dir.clone(),
            "openbabel-3".to_string(),
            bin_dir.clone(),
            bin_dir.join("data"),
        )
    } else {
        let lib_dir = dst.join("lib");
        // OpenBabel installs its plugins/data under a versioned subdirectory
        // (e.g. `lib/openbabel/3.2.1`, `share/openbabel/3.2.1`). Discover the
        // version rather than hardcoding it so a submodule bump keeps working.
        let libdir = versioned_subdir(&lib_dir.join("openbabel"))
            .unwrap_or_else(|| lib_dir.join("openbabel"));
        let datadir = versioned_subdir(&dst.join("share").join("openbabel"))
            .unwrap_or_else(|| dst.join("share").join("openbabel"));
        (lib_dir, "openbabel".to_string(), libdir, datadir)
    };

    // The vendored data files are checked out through git, and on Windows
    // `core.autocrlf=true` — the default there — rewrites every one of them to
    // CRLF. Undo that before anything reads them.
    normalize_data_line_endings(&babel_datadir);

    // Lay the install tree out for packaging as a prebuilt, if asked -- after
    // the normalization above, so an archive never carries CRLF data.
    if let Some(export_dir) = env::var_os(EXPORT_DIR_ENV).filter(|dir| !dir.is_empty()) {
        if prebuilt_dir.is_some() {
            println!(
                "cargo:warning={EXPORT_DIR_ENV} is ignored while {PREBUILT_DIR_ENV} is set: \
                 there is no freshly built tree to export"
            );
        } else {
            export_prebuilt(&dst, Path::new(&export_dir));
        }
    }

    // 2. Compile the cxx bridge + C++ shim against the installed headers.
    //
    // The shim pulls in OpenBabel headers, so it needs the same MSVC settings
    // the library itself was built with: `WIN32`/`_WINDOWS` (enable OpenBabel's
    // Windows shims such as `strcasecmp`), `/EHsc` (the shim uses try/catch),
    // `/GR` (RTTI), and `/utf-8` (our sources contain non-ASCII comments).
    //
    // It also includes `<openbabel/math/align.h>` for OBAlign, and that header
    // includes `<Eigen/Core>` unconditionally, so Eigen must be on the shim's
    // include path too. (No `-DHAVE_EIGEN3` is needed: the header's HAVE_EIGEN3
    // regions are ABI-neutral member declarations, so the shim and library
    // agree on layout without it.)
    let mut build = cxx_build::bridge("src/lib.rs");
    build
        .file("shim/shim.cc")
        .include("shim")
        .include(&include_dir)
        .include(&eigen_dir)
        // OpenBabel's concrete force-field headers (e.g. forcefielduff.h) are
        // private to its source tree, but the Rust term exporter
        // (ff_export_terms) includes one to read a force field's precomputed
        // calculation vectors. Their patched copies live in the install tree
        // (see `install_private_headers`), so a prebuilt carries them too.
        .include(private_include_dir(&dst))
        .std("c++17");
    if env::var("CARGO_CFG_TARGET_ENV").as_deref() == Ok("msvc") {
        build
            .define("WIN32", None)
            .define("_WINDOWS", None)
            .flag("/EHsc")
            .flag("/GR")
            .flag("/utf-8");
    }
    build.compile("obshim");

    // 3. Link the OpenBabel library.
    //    Windows: the import library `openbabel-3.lib` lives in `bin/`.
    //    Unix/macOS: `libopenbabel.{dylib,so}` lives in `lib/`.
    println!("cargo:rustc-link-search=native={}", link_search_dir.display());
    println!("cargo:rustc-link-lib=dylib={}", link_lib_name);

    // On Unix there is no DLL-next-to-exe fallback: the dynamic loader must be
    // able to find `libopenbabel` at runtime. Bake the install directory into
    // the binaries' rpath so `cargo test`/`cargo run` work without the caller
    // having to set DYLD_LIBRARY_PATH / LD_LIBRARY_PATH.
    if !is_windows {
        println!(
            "cargo:rustc-link-arg=-Wl,-rpath,{}",
            link_search_dir.display()
        );
    }

    // 4a. Bake the runtime directories into a generated module so the safe
    //     wrapper can point BABEL_LIBDIR / BABEL_DATADIR at them.
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let generated = format!(
        "// @generated by build.rs — absolute paths into the installed OpenBabel.\n\
         pub const BABEL_LIBDIR: &str = r\"{}\";\n\
         pub const BABEL_DATADIR: &str = r\"{}\";\n",
        babel_libdir.display(),
        babel_datadir.display(),
    );
    fs::write(out_dir.join("paths.rs"), generated).expect("write paths.rs");

    // 4b. Also expose them to dependent build scripts via `links` metadata.
    println!("cargo:babel_libdir={}", babel_libdir.display());
    println!("cargo:babel_datadir={}", babel_datadir.display());

    // 4c. Copy the runtime next to the test/exe binaries. Cargo places those in
    //     target/<profile>/ and target/<profile>/deps/, and Windows resolves
    //     DLLs from the executable's directory.
    //
    //     Crucially, on Windows OpenBabel discovers its format plugins (`.obf`)
    //     in the *directory of openbabel-3.dll* (dlhandler_win32.cpp uses
    //     GetModuleFileName; it does NOT consult BABEL_LIBDIR). So the plugins
    //     must sit beside the copied DLL, not merely be pointed at by an env
    //     var. Data files, by contrast, are found via getenv(BABEL_DATADIR),
    //     which the safe wrapper sets through the C runtime.
    //
    //     We copy every `.dll` (openbabel-3.dll plus its runtime dependencies
    //     such as the bundled `inchi` library, which `inchiformat.obf` links
    //     against) and every `.obf` plugin.
    //
    //     Unix/macOS need none of this: the main library is found via the rpath
    //     baked in above, and the plugins are found through the BABEL_LIBDIR env
    //     var (consulted by dlhandler_unix.cpp), so we skip the copy entirely.
    if is_windows {
        if let Some(profile_dir) = out_dir.ancestors().nth(3) {
            let mut runtime = Vec::new();
            if let Ok(entries) = fs::read_dir(&link_search_dir) {
                for entry in entries.flatten() {
                    let p = entry.path();
                    match p.extension().and_then(|e| e.to_str()) {
                        Some("dll") | Some("obf") => runtime.push(p),
                        _ => {}
                    }
                }
            }
            for dest_dir in [profile_dir.to_path_buf(), profile_dir.join("deps")] {
                let _ = fs::create_dir_all(&dest_dir);
                for src in &runtime {
                    if let Some(name) = src.file_name() {
                        copy_if_changed(src, &dest_dir.join(name));
                    }
                }
            }
        }
    }
}

/// Build OpenBabel from the vendored sources and install it into `OUT_DIR`,
/// returning the install root.
fn build_openbabel(ob_src: &Path, eigen_dir: &Path, key: Option<&str>) -> PathBuf {
    // Fail early with an actionable message if cmake is missing. The `cmake`
    // crate shells out to the `cmake` binary and, when it is absent, panics
    // deep inside the crate with a generic "is `cmake` not installed?" — which
    // buries the one thing the user needs: the command to install it. Cargo has
    // no way to install system build tools itself, so we can only guide.
    ensure_cmake_present();

    // Make the bundled InChI library linkable on MSVC (see fn docs).
    ensure_inchi_auxinfo_stubs(ob_src);

    // Add read-only accessors to the force-field headers so the Rust term
    // exporter can read their precomputed calculation vectors (see fn docs).
    ensure_ff_accessors(ob_src);

    // We force a Release build regardless of the cargo profile so the C++ side
    // always uses the release CRT (/MD on MSVC), matching what the cc/cxx-build
    // compiled shim uses — mixing CRTs would be an ABI hazard. It also avoids
    // needing debug builds of any optional dependency.
    //
    // Optional features that pull in heavy external deps (Boost/Cairo/
    // RapidJSON) are disabled; the formats we need do not require them.
    //
    // Eigen IS enabled (via EIGEN3_INCLUDE_DIR below): it defines HAVE_EIGEN3,
    // which compiles OBAlign (least-squares structure superposition) and the
    // `align`/`conformer` ops, and enables distance-geometry 3D generation.
    //
    // InChI IS enabled: OpenBabel bundles the InChI library source and builds
    // it into a separate `inchi` shared library, so this needs no external
    // dependency and unlocks InChI / InChIKey output.
    //
    // Two InChI-specific knobs matter on MSVC. OpenBabel's CMakeLists, when
    // `OB_USE_PREBUILT_BINARIES` is ON (its default on MSVC), *force*-sets
    // `OPENBABEL_USE_SYSTEM_INCHI` ON in the cache — which then makes
    // `find_package(Inchi REQUIRED)` fail (there is no system InChI here) and
    // aborts configuration. Since we build everything from source (no prebuilt
    // binaries), we turn that flag OFF, which lets `OPENBABEL_USE_SYSTEM_INCHI`
    // stay OFF so the bundled InChI source is compiled instead. We also set
    // `OPENBABEL_USE_SYSTEM_INCHI=OFF` explicitly to overwrite any stale ON
    // value a previous forced configure may have baked into the CMake cache.
    // CMake prefers forward slashes in `-D` path values on Windows.
    let eigen_include = eigen_dir.to_string_lossy().replace('\\', "/");
    let mut cfg = cmake::Config::new(ob_src);
    cfg.profile("Release")
        .define("BUILD_GUI", "OFF")
        .define("ENABLE_TESTS", "OFF")
        .define("BUILD_SHARED", "ON")
        .define("WITH_MAEPARSER", "OFF")
        .define("WITH_COORDGEN", "OFF")
        .define("WITH_JSON", "OFF")
        .define("WITH_INCHI", "ON")
        .define("OB_USE_PREBUILT_BINARIES", "OFF")
        .define("OPENBABEL_USE_SYSTEM_INCHI", "OFF")
        .define("EIGEN3_INCLUDE_DIR", &eigen_include);

    // On Linux and macOS OpenBabel bakes its install prefix into what it
    // installs -- the rpath on Linux, the install name on macOS -- so the tree
    // only loads where it was built: a prebuilt unpacked anywhere else fails,
    // on Linux as `inchiformat` quietly not finding `libinchi`, on macOS as the
    // binary not starting at all. OpenBabel's wheel build switches to
    // loader-relative paths (`$ORIGIN`, `@rpath`), and `BUILD_BY_PIP` is the
    // only way in: CMakeLists sets CMAKE_INSTALL_RPATH as a normal variable,
    // which shadows any `-D`.
    //
    // With the Python bindings off the flag does two other things. It skips
    // Cairo, which only PNG depiction uses and which would tie the tree to the
    // build machine's Cairo. And it runs Python to ask `cmeel` for a CMake
    // prefix -- disabled here, so a stray cmeel install cannot redirect
    // `find_package(Eigen3)` away from the vendored copy.
    if env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        cfg.define("BUILD_BY_PIP", "ON")
            .define("CMAKE_DISABLE_FIND_PACKAGE_Python", "ON");
    }

    // OpenBabel's `data/CMakeLists.txt` shells out to `bin2hex.pl` to compile
    // its data tables into headers, so Perl is a hard build requirement. On
    // Windows the only Perl most machines have is the one inside Git, which Git
    // deliberately keeps off the system PATH -- so the same machine builds fine
    // from Git Bash and fails from PowerShell with "Could NOT find Perl". Hand
    // CMake the path when we can find one it would not.
    if let Some(perl) = find_perl() {
        cfg.define(
            "PERL_EXECUTABLE",
            perl.to_string_lossy().replace(char::from(92u8), "/"),
        );
    }

    // The `cmake` crate overrides CMAKE_CXX_FLAGS/CMAKE_C_FLAGS with its own
    // minimal set (`-nologo -MD -Brepro -W0`), which drops the `/DWIN32
    // /D_WINDOWS /EHsc /GR` that CMake's MSVC platform module normally adds.
    // OpenBabel relies on `WIN32` being defined to enable its Windows shims
    // (e.g. `#define strcasecmp _stricmp` in babelconfig.h), and needs C++
    // exceptions (/EHsc) and RTTI (/GR). Restore them here so the build
    // matches a stock CMake MSVC configuration.
    if env::var("CARGO_CFG_TARGET_ENV").as_deref() == Ok("msvc") {
        for flag in ["/DWIN32", "/D_WINDOWS", "/EHsc", "/GR"] {
            cfg.cxxflag(flag);
        }
        for flag in ["/DWIN32", "/D_WINDOWS"] {
            cfg.cflag(flag);
        }
    }

    let dst = cfg.build();

    install_private_headers(ob_src, &dst);

    // Record what the tree was built from, for whoever uses an export of it as
    // a prebuilt. Without a key there is nothing to record, and a file left by
    // an earlier build would vouch for sources it never saw.
    let key_path = dst.join(KEY_FILE);
    match key {
        Some(key) => {
            if fs::read_to_string(&key_path).ok().as_deref() != Some(key) {
                fs::write(&key_path, key).expect("write the prebuilt key file");
            }
        }
        None => {
            let _ = fs::remove_file(&key_path);
        }
    }

    dst
}

/// Take OpenBabel from an install tree built elsewhere instead of compiling it,
/// returning a copy of the tree inside `OUT_DIR`.
///
/// The tree has to be exactly what this checkout would build: the shim is still
/// compiled here, against the tree's headers, and links into its library, so a
/// tree from another version fails to link at best and at worst links and
/// misbehaves. [`prebuilt_key`] pins down "what this checkout would build", and
/// anything short of an exact match is refused. Falling back to a source build
/// instead would turn a mismatched download into a 20-minute build nobody asked
/// for, with nothing in the log to say why.
///
/// The copy is not for tidiness. The rpath this script bakes in reaches only
/// openbabel-sys's own binaries; a dependent's tests and binaries find
/// `libopenbabel` through the loader path cargo sets for `cargo test` and
/// `cargo run`, and cargo puts a link-search directory on it only when that
/// directory lies under `target/<profile>/`. A tree used where it was unpacked
/// links fine and then fails to load on Linux and macOS. The copy also leaves
/// nothing pointing at the unpacked archive, which CI tends to delete.
///
/// None of the source build's requirements apply: no CMake, no Perl, and the
/// vendored OpenBabel tree is left unpatched.
fn use_prebuilt(dir: PathBuf, key: Option<&str>) -> PathBuf {
    let key_path = dir.join(KEY_FILE);
    println!("cargo:rerun-if-changed={}", key_path.display());

    let found = fs::read_to_string(&key_path).unwrap_or_else(|err| {
        panic!(
            "\n{PREBUILT_DIR_ENV} is {}, but {KEY_FILE} cannot be read there ({err}).\n\
             Point it at an unpacked openbabel-sys prebuilt archive, or unset it to \
             build OpenBabel from source.\n",
            dir.display(),
        )
    });
    let Some(key) = key else {
        panic!(
            "\n{PREBUILT_DIR_ENV} is set, but the commits checked out in vendor/openbabel-src \
             and vendor/eigen cannot be determined, so the prebuilt cannot be checked against \
             them. That needs `git` on PATH and submodules checked out by git.\n\
             Unset {PREBUILT_DIR_ENV} to build OpenBabel from source.\n"
        )
    };
    if found != key {
        panic!(
            "\nThe OpenBabel prebuilt at {} was not built from this openbabel-sys.\n\n\
             The prebuilt was built from:\n{found}\n\
             This checkout needs:\n{key}\n\
             Use the archive attached to the release this openbabel-sys comes from, \
             or unset {PREBUILT_DIR_ENV} to build OpenBabel from source.\n",
            dir.display(),
        )
    }

    // A copy that finished carries the key file (it goes last), so one cut
    // short is redone rather than trusted.
    let dst = PathBuf::from(env::var("OUT_DIR").unwrap()).join("prebuilt");
    if fs::read_to_string(dst.join(KEY_FILE)).ok().as_deref() != Some(key) {
        let _ = fs::remove_dir_all(&dst);
        copy_install_tree(&dir, &dst);
    }

    dst
}

/// Describe what an OpenBabel install tree is built from, so a prebuilt can be
/// checked against the checkout about to use it.
///
/// Two checkouts produce interchangeable trees when they agree on the target,
/// on the OpenBabel and Eigen commits, and on this build script -- which holds
/// the CMake options and every patch applied to the OpenBabel tree, so a change
/// to either shows up in its hash. The script is hashed without CR bytes: a
/// Windows checkout and cargo's own checkout of the same commit can disagree on
/// line endings.
///
/// `None` when a submodule's commit cannot be determined.
fn prebuilt_key(manifest_dir: &Path, ob_src: &Path, eigen_dir: &Path) -> Option<String> {
    const CR: u8 = 13;

    let mut script = fs::read(manifest_dir.join("build.rs")).ok()?;
    script.retain(|&b| b != CR);

    Some(format!(
        "openbabel-sys prebuilt 1\n\
         target {}\n\
         openbabel-src {}\n\
         eigen {}\n\
         build.rs {:016x}\n",
        env::var("TARGET").unwrap(),
        submodule_commit(ob_src)?,
        submodule_commit(eigen_dir)?,
        fnv1a64(&script),
    ))
}

/// The commit checked out in `dir`, which must be a git checkout of its own.
///
/// Only asked when `dir` has its own `.git` -- a file in a submodule, a
/// directory in cargo's checkouts of a git dependency. Without one, git would
/// walk up and answer for whichever repository happens to contain `dir`.
fn submodule_commit(dir: &Path) -> Option<String> {
    if !dir.join(".git").exists() {
        return None;
    }
    let out = Command::new("git")
        .arg("-C")
        .arg(dir)
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())?;
    let commit = String::from_utf8(out.stdout).ok()?.trim().to_owned();
    (!commit.is_empty() && commit.bytes().all(|b| b.is_ascii_hexdigit())).then_some(commit)
}

/// 64-bit FNV-1a. Spelled out rather than taken from std, whose hashers may
/// change between Rust releases -- and a prebuilt is made by one toolchain and
/// checked by another.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &b| {
        (hash ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

/// Copy the install tree at `dst` into `dir` in the layout [`use_prebuilt`]
/// takes: CMake's install directories plus the key file. CMake's intermediate
/// `build` directory -- nearly all of `OUT_DIR` -- stays behind.
///
/// `dir` is cleared first when it holds an earlier export, and refused when it
/// holds anything else: it comes from an environment variable, and wiping
/// whatever that happens to name is not a build script's call.
fn export_prebuilt(dst: &Path, dir: &Path) {
    assert!(
        dst.join(KEY_FILE).is_file(),
        "{EXPORT_DIR_ENV} is set, but this build has no {KEY_FILE}: the commits checked out \
         in vendor/openbabel-src and vendor/eigen could not be determined (is git on PATH?)"
    );
    if fs::read_dir(dir).is_ok_and(|mut entries| entries.next().is_some()) {
        assert!(
            dir.join(KEY_FILE).is_file(),
            "{EXPORT_DIR_ENV} is {}, which is not empty and holds no earlier export; \
             refusing to overwrite it",
            dir.display()
        );
        fs::remove_dir_all(dir).expect("clear the previous prebuilt export");
    }
    copy_install_tree(dst, dir);
}

/// Copy what makes up a prebuilt from the install tree at `from` into `to`:
/// CMake's install directories, then the key file -- last, so that only a
/// finished copy carries it.
fn copy_install_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).expect("create a prebuilt directory");
    for name in ["bin", "include", "lib", "share", KEY_FILE] {
        let path = from.join(name);
        if path.exists() {
            copy_tree(&path, &to.join(name));
        }
    }
}

/// Recursively copy `from` to `to`, recreating symlinks as symlinks where the
/// host has them -- a Unix install links `libopenbabel.so` to its versioned
/// name, and following those links would ship the library three times.
fn copy_tree(from: &Path, to: &Path) {
    let meta = fs::symlink_metadata(from).expect("stat an install tree entry");
    if meta.is_dir() {
        fs::create_dir_all(to).expect("create a prebuilt export directory");
        for entry in fs::read_dir(from)
            .expect("list an install tree directory")
            .flatten()
        {
            copy_tree(&entry.path(), &to.join(entry.file_name()));
        }
        return;
    }
    #[cfg(unix)]
    if meta.file_type().is_symlink() {
        let target = fs::read_link(from).expect("read an install tree symlink");
        std::os::unix::fs::symlink(target, to).expect("recreate a symlink in the export");
        return;
    }
    fs::copy(from, to).expect("copy an install tree file");
}

/// Where the patched private force-field headers sit in an install tree.
fn private_include_dir(root: &Path) -> PathBuf {
    root.join("include")
        .join("openbabel-sys")
        .join("forcefields")
}

/// Copy OpenBabel's force-field headers into the install tree at `dst`.
///
/// OpenBabel never installs them -- they are private to its source tree -- but
/// the shim includes them, patched by [`ensure_ff_accessors`]. Keeping the
/// patched copies in the install tree gives both build modes one include path,
/// and a prebuilt needs no OpenBabel sources beside it.
fn install_private_headers(ob_src: &Path, dst: &Path) {
    let from = ob_src.join("src").join("forcefields");
    let to = private_include_dir(dst);
    fs::create_dir_all(&to).expect("create the private header directory");
    for entry in fs::read_dir(&from)
        .expect("list OpenBabel's force-field sources")
        .flatten()
    {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("h") {
            copy_if_changed(&path, &to.join(entry.file_name()));
        }
    }
}

/// Return the single versioned subdirectory of `parent` (e.g. the `3.2.1` in
/// `lib/openbabel/3.2.1`), or `None` if `parent` has no subdirectories.
///
/// OpenBabel installs its plugins and data under a directory named after the
/// library version. We discover it at build time rather than hardcoding a
/// version so bumping the `vendor/openbabel-src` submodule needs no edit here.
/// If several subdirectories exist, the lexicographically greatest is chosen.
fn versioned_subdir(parent: &Path) -> Option<PathBuf> {
    fs::read_dir(parent)
        .ok()?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_dir())
        .max()
}

/// Verify the `cmake` binary is on `PATH`, aborting with an OS-specific install
/// hint if not.
///
/// Building OpenBabel from source requires cmake, but Cargo cannot install
/// system build tools. Without this check the build panics deep inside the
/// `cmake` crate with a message that omits *how* to get cmake; here we surface
/// the exact command for the host platform instead. Honors the `CMAKE`
/// environment variable, which the `cmake` crate itself consults to locate a
/// non-`PATH` binary.
fn ensure_cmake_present() {
    let cmake = env::var_os("CMAKE").unwrap_or_else(|| "cmake".into());
    let found = Command::new(&cmake)
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if found {
        return;
    }

    let hint = match env::var("CARGO_CFG_TARGET_OS").as_deref() {
        Ok("macos") => "Install it with:  brew install cmake",
        Ok("windows") => {
            "Install it with:  winget install Kitware.CMake\n\
             (or download from https://cmake.org/download/ and add it to PATH)"
        }
        Ok("linux") => {
            "Install it with your package manager, e.g.:\n\
             \x20 Debian/Ubuntu:  sudo apt-get install cmake\n\
             \x20 Fedora:         sudo dnf install cmake\n\
             \x20 Arch:           sudo pacman -S cmake"
        }
        _ => "Install cmake from https://cmake.org/download/ and ensure it is on PATH",
    };

    panic!(
        "\n\
         `cmake` was not found on PATH, but it is required to build OpenBabel \
         from source.\n\
         {hint}\n\
         (Or set the CMAKE environment variable to the full path of a cmake binary.)\n"
    );
}

/// Write inert implementations of four InChI AuxInfo entry points into the
/// bundled InChI source tree so the `inchi` shared library links on MSVC.
///
/// OpenBabel vendors a trimmed IUPAC InChI library under
/// `src/formats/libinchi/`. Its `inchi_dll.c` compiles the `cdecl_`/`pasc_`
/// ABI wrappers, which reference four AuxInfo-parsing functions —
/// `Get_inchi_Input_FromAuxInfo`, `Get_std_inchi_Input_FromAuxInfo`,
/// `Free_inchi_Input`, `Free_std_inchi_Input` — whose *implementations* were
/// never vendored. On Unix the shared library links regardless (unresolved
/// symbols are permitted); on MSVC the DLL link is fatal (LNK2019/LNK1120).
///
/// OpenBabel never calls these (they turn an InChI AuxInfo string back into a
/// structure, unrelated to the InChI/InChIKey *output* this crate exposes), so
/// inert stubs are safe. The file lands in the `libinchi` directory where the
/// target's `file(GLOB *.c)` will pick it up. It is written only when missing
/// or stale to avoid needless reconfigure/rebuild churn.
fn ensure_inchi_auxinfo_stubs(ob_src: &Path) {
    let libinchi = ob_src.join("src").join("formats").join("libinchi");
    if !libinchi.join("inchi_dll.c").exists() {
        // InChI source not present (e.g. WITH_INCHI disabled downstream) —
        // nothing to patch.
        return;
    }
    let stub_path = libinchi.join("ob_rs_auxinfo_stubs.c");
    let stub = "\
/* @generated by openbabel-sys/build.rs — do not edit.\n\
 *\n\
 * Inert implementations of four InChI AuxInfo entry points that OpenBabel's\n\
 * vendored InChI references (via the cdecl_/pasc_ wrappers in inchi_dll.c) but\n\
 * does not ship. Without them the `inchi` DLL fails to link on MSVC. OpenBabel\n\
 * never calls these, so returning failure is safe. See build.rs for details.\n\
 */\n\
#include \"inchi_api.h\"\n\
\n\
INCHI_API int INCHI_DECL Get_inchi_Input_FromAuxInfo(\n\
    char *szInchiAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo,\n\
    InchiInpData *pInchiInp)\n\
{\n\
    (void)szInchiAuxInfo; (void)bDoNotAddH; (void)bDiffUnkUndfStereo;\n\
    (void)pInchiInp;\n\
    return -1; /* failure: AuxInfo parsing is unsupported in this build */\n\
}\n\
\n\
INCHI_API int INCHI_DECL Get_std_inchi_Input_FromAuxInfo(\n\
    char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp)\n\
{\n\
    (void)szInchiAuxInfo; (void)bDoNotAddH; (void)pInchiInp;\n\
    return -1;\n\
}\n\
\n\
INCHI_API void INCHI_DECL Free_inchi_Input(inchi_Input *pInp)\n\
{\n\
    (void)pInp;\n\
}\n\
\n\
INCHI_API void INCHI_DECL Free_std_inchi_Input(inchi_Input *pInp)\n\
{\n\
    (void)pInp;\n\
}\n";
    let up_to_date = fs::read_to_string(&stub_path)
        .map(|existing| existing == stub)
        .unwrap_or(false);
    if !up_to_date {
        fs::write(&stub_path, stub).expect("write InChI AuxInfo stub source");
    }
}

/// Inject read-only accessors for each force field's precomputed calculation
/// vectors into its header, so the shim's `ff_export_terms` can read them.
///
/// The `_bondcalculations` / … vectors are `protected` members of the concrete
/// `OBForceFieldXXX` classes. Rather than reimplement OpenBabel's parameter
/// setup in Rust, we let OpenBabel precompute the terms and read them out — but
/// that needs access. We add `public` **inline** accessors returning const
/// references; being inline, they add no out-of-line symbols and do not change
/// the class layout or vtable, so the already-built OpenBabel library's ABI is
/// unaffected (the accessors are compiled only into the shim translation unit).
///
/// Idempotent (keyed on a marker), matching the write-if-changed style of
/// [`ensure_inchi_auxinfo_stubs`]. A fresh `git submodule update` resets the
/// headers; the next build re-applies them. New force fields add an entry here.
fn ensure_ff_accessors(ob_src: &Path) {
    let dir = ob_src.join("src").join("forcefields");

    patch_ff_header(
        &dir.join("forcefielduff.h"),
        "  }; // class OBForceFieldUFF",
        "\
    const std::vector<OBFFBondCalculationUFF>&          RsBondCalcs()    const { return _bondcalculations; }\n\
    const std::vector<OBFFAngleCalculationUFF>&         RsAngleCalcs()   const { return _anglecalculations; }\n\
    const std::vector<OBFFTorsionCalculationUFF>&       RsTorsionCalcs() const { return _torsioncalculations; }\n\
    const std::vector<OBFFOOPCalculationUFF>&           RsOopCalcs()     const { return _oopcalculations; }\n\
    const std::vector<OBFFVDWCalculationUFF>&           RsVdwCalcs()     const { return _vdwcalculations; }\n\
    const std::vector<OBFFElectrostaticCalculationUFF>& RsElecCalcs()    const { return _electrostaticcalculations; }\n",
    );

    patch_ff_header(
        &dir.join("forcefieldghemical.h"),
        "  }; // class OBForceFieldGhemical",
        "\
    const std::vector<OBFFBondCalculationGhemical>&          RsBondCalcs()    const { return _bondcalculations; }\n\
    const std::vector<OBFFAngleCalculationGhemical>&         RsAngleCalcs()   const { return _anglecalculations; }\n\
    const std::vector<OBFFTorsionCalculationGhemical>&       RsTorsionCalcs() const { return _torsioncalculations; }\n\
    const std::vector<OBFFVDWCalculationGhemical>&           RsVdwCalcs()     const { return _vdwcalculations; }\n\
    const std::vector<OBFFElectrostaticCalculationGhemical>& RsElecCalcs()    const { return _electrostaticcalculations; }\n",
    );

    patch_ff_header(
        &dir.join("forcefieldgaff.h"),
        "  }; // class OBForceFieldGaff",
        "\
    const std::vector<OBFFBondCalculationGaff>&          RsBondCalcs()    const { return _bondcalculations; }\n\
    const std::vector<OBFFAngleCalculationGaff>&         RsAngleCalcs()   const { return _anglecalculations; }\n\
    const std::vector<OBFFTorsionCalculationGaff>&       RsTorsionCalcs() const { return _torsioncalculations; }\n\
    const std::vector<OBFFOOPCalculationGaff>&           RsOopCalcs()     const { return _oopcalculations; }\n\
    const std::vector<OBFFVDWCalculationGaff>&           RsVdwCalcs()     const { return _vdwcalculations; }\n\
    const std::vector<OBFFElectrostaticCalculationGaff>& RsElecCalcs()    const { return _electrostaticcalculations; }\n",
    );

    // forcefieldmmff94.h holds both MMFF94 and MMFF94s (one class, `mmff94s`
    // flag); its class-closing marker is an upstream copy-paste of MM2's, but is
    // unique within this file. MMFF94 adds a stretch-bend calculation vector.
    patch_ff_header(
        &dir.join("forcefieldmmff94.h"),
        "  }; // class OBForceFieldMM2",
        "\
    const std::vector<OBFFBondCalculationMMFF94>&          RsBondCalcs()    const { return _bondcalculations; }\n\
    const std::vector<OBFFAngleCalculationMMFF94>&         RsAngleCalcs()   const { return _anglecalculations; }\n\
    const std::vector<OBFFStrBndCalculationMMFF94>&        RsStrBndCalcs()  const { return _strbndcalculations; }\n\
    const std::vector<OBFFTorsionCalculationMMFF94>&       RsTorsionCalcs() const { return _torsioncalculations; }\n\
    const std::vector<OBFFOOPCalculationMMFF94>&           RsOopCalcs()     const { return _oopcalculations; }\n\
    const std::vector<OBFFVDWCalculationMMFF94>&           RsVdwCalcs()     const { return _vdwcalculations; }\n\
    const std::vector<OBFFElectrostaticCalculationMMFF94>& RsElecCalcs()    const { return _electrostaticcalculations; }\n",
    );

    // MM2 predates the OBFFCalculation architecture: it has no calc vectors and
    // resolves each term's parameters inline during energy evaluation. So the
    // exporter needs the raw parameter tables, the global unit constants, and a
    // forwarder to OpenBabel's (protected) parameter lookup — the shim then
    // replicates MM2's own iteration to emit resolved terms.
    patch_ff_header(
        &dir.join("forcefieldmm2.h"),
        "  }; // class OBForceFieldMM2",
        "\
    OBMol& RsMol() { return _mol; }\n\
    double RsBondUnit() const { return bondunit; }\n\
    double RsBondCubic() const { return bond_cubic; }\n\
    double RsBondQuartic() const { return bond_quartic; }\n\
    double RsAngleUnit() const { return angleunit; }\n\
    double RsAngleSextic() const { return angle_sextic; }\n\
    double RsStretchBendUnit() const { return stretchbendunit; }\n\
    double RsTorsionUnit() const { return torsionunit; }\n\
    double RsOutPlaneBendUnit() const { return outplanebendunit; }\n\
    double RsAExpterm() const { return a_expterm; }\n\
    double RsBExpterm() const { return b_expterm; }\n\
    double RsCExpterm() const { return c_expterm; }\n\
    double RsDielectric() const { return dielectric; }\n\
    std::vector<OBFFParameter>& RsVec(int which) {\n\
      switch (which) {\n\
      case 0: return _ffbondparams;\n\
      case 1: return _ffangleparams;\n\
      case 2: return _ffstretchbendparams;\n\
      case 3: return _fftorsionparams;\n\
      case 4: return _ffoutplanebendparams;\n\
      case 5: return _ffvdwprparams;\n\
      case 6: return _ffvdwparams;\n\
      default: return _ffdipoleparams;\n\
      }\n\
    }\n\
    OBFFParameter* RsParam(int a, int b, int c, int d, int which) { return GetParameter(a, b, c, d, RsVec(which)); }\n\
    int RsParamIdx(int a, int b, int c, int d, int which) { return GetParameterIdx(a, b, c, d, RsVec(which)); }\n",
    );
}

/// Insert a `public:` accessor block just before `marker` (a class's closing
/// brace) in `header`, unless already present. See [`ensure_ff_accessors`].
fn patch_ff_header(header: &Path, marker: &str, accessors: &str) {
    let Ok(src) = fs::read_to_string(header) else {
        return; // source not present — nothing to patch.
    };
    if src.contains("@generated by openbabel-sys/build.rs") {
        return; // already patched (marker present in every generated block)
    }
    let Some(pos) = src.find(marker) else {
        return; // upstream layout changed; skip rather than corrupt the header.
    };
    let block = format!(
        "  public:\n    // @generated by openbabel-sys/build.rs — read-only accessors for the\n    // Rust force-field term exporter (ff_export_terms). Inline: no ABI change.\n{accessors}\n"
    );
    let mut patched = String::with_capacity(src.len() + block.len());
    patched.push_str(&src[..pos]);
    patched.push_str(&block);
    patched.push_str(&src[pos..]);
    fs::write(header, patched).expect("patch force-field header with Rust accessors");
}

/// Locate a Perl interpreter for OpenBabel's data-header generation.
///
/// `None` when Perl is already on `PATH`, since CMake will find that one
/// itself, and also when nothing turns up -- in which case CMake fails with its
/// own message, which says what is missing.
fn find_perl() -> Option<PathBuf> {
    if Command::new("perl")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        return None;
    }

    // Derive Git's install root from where its helpers live, rather than
    // guessing at Program Files: `git --exec-path` prints something like
    // `C:/Program Files/Git/mingw64/libexec/git-core`, and Perl sits at
    // `<root>/usr/bin/perl.exe`.
    let from_git = Command::new("git")
        .arg("--exec-path")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| {
            let exec_path = PathBuf::from(String::from_utf8_lossy(&o.stdout).trim());
            exec_path
                .ancestors()
                .map(|root| root.join("usr").join("bin").join("perl.exe"))
                .find(|p| p.is_file())
        });
    if from_git.is_some() {
        return from_git;
    }

    ["C:/Program Files/Git/usr/bin/perl.exe"]
        .iter()
        .map(PathBuf::from)
        .find(|p| p.is_file())
}

/// Rewrite CRLF back to LF across the installed data directory.
///
/// `rigid-fragments-index.txt` and `ring-fragments-index.txt` are *byte-offset*
/// indexes into their fragment files, generated upstream against files with LF
/// endings. Checked out with CRLF every offset is short by the number of
/// preceding lines — measured at 46 190 bytes for a fragment on line 46 210 —
/// so `OBBuilder::GetFragmentCoord` seeks into the middle of an unrelated
/// fragment and parses whatever is there. It reports that as
///
/// ```text
/// Rigid fragment O=C1CC(=O)NC(=O)N1 in rigid-fragments.txt has all zero coordinates.
/// ```
///
/// and then builds the molecule anyway, piling the atoms it could not place on
/// the origin or writing NaN coordinates — while `generate_3d()` still returns
/// true. Measured on 3-phenyl-o-benzoquinone: 14 of 40 builds usable with the
/// CRLF data, 40 of 40 after this normalization.
///
/// Every data file OpenBabel ships is text bar one PNG, and the parsers are
/// written for LF, so the whole directory is normalized rather than just the
/// two fragment files — any other offset- or column-sensitive table would fail
/// the same way. Files holding a NUL byte are left alone: that is the PNG, and
/// anything else binary that might be added later.
fn normalize_data_line_endings(dir: &Path) {
    // Spelled as byte values rather than escapes: this function is about CR
    // and LF, and a source escape for them is the one thing an editor or a
    // checkout can quietly rewrite.
    const CR: u8 = 13;
    const LF: u8 = 10;

    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Ok(bytes) = fs::read(&path) else {
            continue;
        };
        // A NUL means binary — the shipped PNG, and whatever else may join it.
        let has_crlf = bytes.windows(2).any(|w| w[0] == CR && w[1] == LF);
        if bytes.contains(&0) || !has_crlf {
            continue;
        }
        let mut out = Vec::with_capacity(bytes.len());
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] == CR && bytes.get(i + 1) == Some(&LF) {
                i += 1;
                continue;
            }
            out.push(bytes[i]);
            i += 1;
        }
        let _ = fs::write(&path, out);
    }
}

/// Copy `src` to `dst` unless `dst` already holds the same bytes.
///
/// Compared by content, not modification time: an unpacked prebuilt keeps the
/// times it was packed with, easily older than a DLL an earlier source build
/// left in `target/` -- and a time check would then keep that stale DLL beside
/// the new plugins.
fn copy_if_changed(src: &Path, dst: &Path) {
    let same = match (fs::metadata(src), fs::metadata(dst)) {
        (Ok(s), Ok(d)) if s.len() == d.len() => fs::read(src)
            .ok()
            .zip(fs::read(dst).ok())
            .is_some_and(|(s, d)| s == d),
        _ => false,
    };
    if !same {
        let _ = fs::copy(src, dst);
    }
}
