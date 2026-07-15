//! Build script for `openbabel-sys`.
//!
//! Pipeline:
//!   1. Build + install OpenBabel (from the `vendor/openbabel-src` submodule)
//!      into `OUT_DIR` via the `cmake` crate.
//!   2. Compile the cxx bridge (`src/lib.rs`) + C++ shim (`shim/shim.cc`),
//!      pointing the C++ compiler at the freshly installed OpenBabel headers.
//!   3. Link against the OpenBabel import library.
//!   4. Make the runtime discoverable: bake the plugin/data directories into a
//!      generated `paths.rs`, and copy `openbabel-3.dll` next to the eventual
//!      test/exe binaries (Windows has no rpath).
//!
//! The first build compiles all of OpenBabel and is slow (~10-20 min); later
//! builds are incremental (cmake no-ops) unless the shim or bridge changes.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let workspace_root = manifest_dir.parent().unwrap();
    let ob_src = workspace_root.join("vendor").join("openbabel-src");

    println!("cargo:rerun-if-changed=shim/shim.cc");
    println!("cargo:rerun-if-changed=shim/shim.h");
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=build.rs");

    assert!(
        ob_src.join("CMakeLists.txt").exists(),
        "OpenBabel source not found at {}.\n\
         Run: git submodule update --init --recursive",
        ob_src.display()
    );

    // 1. Build + install OpenBabel into OUT_DIR.
    //
    // We force a Release build regardless of the cargo profile so the C++ side
    // always uses the release CRT (/MD on MSVC), matching what the cc/cxx-build
    // compiled shim uses — mixing CRTs would be an ABI hazard. It also avoids
    // needing debug builds of any optional dependency.
    //
    // Optional features that pull in heavy external deps (Boost/Eigen/Cairo/
    // InChI/RapidJSON) are disabled for the MVP; the core formats we need
    // (SMILES/MOL/SDF/PDB, in formats_common) do not require them.
    let mut cfg = cmake::Config::new(&ob_src);
    cfg.profile("Release")
        .define("BUILD_GUI", "OFF")
        .define("ENABLE_TESTS", "OFF")
        .define("BUILD_SHARED", "ON")
        .define("WITH_MAEPARSER", "OFF")
        .define("WITH_COORDGEN", "OFF")
        .define("WITH_JSON", "OFF")
        .define("WITH_INCHI", "OFF");

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

    // Installed layout (Windows):
    //   <dst>/include/openbabel3/openbabel/*.h   headers
    //   <dst>/bin/openbabel-3.lib                import library
    //   <dst>/bin/openbabel-3.dll                runtime library
    //   <dst>/bin/*.obf                          format/plugin modules (BABEL_LIBDIR)
    //   <dst>/bin/data/                          runtime data (BABEL_DATADIR)
    let include_dir = dst.join("include").join("openbabel3");
    let bin_dir = dst.join("bin");
    let data_dir = bin_dir.join("data");

    assert!(
        include_dir.join("openbabel").join("mol.h").exists(),
        "expected OpenBabel headers at {}",
        include_dir.display()
    );

    // 2. Compile the cxx bridge + C++ shim against the installed headers.
    //
    // The shim pulls in OpenBabel headers, so it needs the same MSVC settings
    // the library itself was built with: `WIN32`/`_WINDOWS` (enable OpenBabel's
    // Windows shims such as `strcasecmp`), `/EHsc` (the shim uses try/catch),
    // `/GR` (RTTI), and `/utf-8` (our sources contain non-ASCII comments).
    let mut build = cxx_build::bridge("src/lib.rs");
    build
        .file("shim/shim.cc")
        .include("shim")
        .include(&include_dir)
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

    // 3. Link the OpenBabel import library (openbabel-3.lib lives in bin/).
    println!("cargo:rustc-link-search=native={}", bin_dir.display());
    println!("cargo:rustc-link-lib=dylib=openbabel-3");

    // 4a. Bake the runtime directories into a generated module so the safe
    //     wrapper can point BABEL_LIBDIR / BABEL_DATADIR at them.
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let generated = format!(
        "// @generated by build.rs — absolute paths into the installed OpenBabel.\n\
         pub const BABEL_LIBDIR: &str = r\"{}\";\n\
         pub const BABEL_DATADIR: &str = r\"{}\";\n",
        bin_dir.display(),
        data_dir.display(),
    );
    fs::write(out_dir.join("paths.rs"), generated).expect("write paths.rs");

    // 4b. Also expose them to dependent build scripts via `links` metadata.
    println!("cargo:babel_libdir={}", bin_dir.display());
    println!("cargo:babel_datadir={}", data_dir.display());

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
    if let Some(profile_dir) = out_dir.ancestors().nth(3) {
        let mut runtime = vec![bin_dir.join("openbabel-3.dll")];
        if let Ok(entries) = fs::read_dir(&bin_dir) {
            for entry in entries.flatten() {
                let p = entry.path();
                if p.extension().and_then(|e| e.to_str()) == Some("obf") {
                    runtime.push(p);
                }
            }
        }
        for dest_dir in [profile_dir.to_path_buf(), profile_dir.join("deps")] {
            let _ = fs::create_dir_all(&dest_dir);
            for src in &runtime {
                if let Some(name) = src.file_name() {
                    copy_if_newer(src, &dest_dir.join(name));
                }
            }
        }
    }
}

/// Copy `src` to `dst`, skipping the write if `dst` is already up to date.
fn copy_if_newer(src: &Path, dst: &Path) {
    if !src.exists() {
        return;
    }
    let up_to_date = fs::metadata(dst)
        .ok()
        .zip(fs::metadata(src).ok())
        .and_then(|(d, s)| Some((d.modified().ok()?, s.modified().ok()?)))
        .map(|(dm, sm)| dm >= sm)
        .unwrap_or(false);
    if !up_to_date {
        let _ = fs::copy(src, dst);
    }
}
