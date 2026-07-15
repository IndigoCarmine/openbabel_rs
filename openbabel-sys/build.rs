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

    // Make the bundled InChI library linkable on MSVC (see fn docs).
    ensure_inchi_auxinfo_stubs(&ob_src);

    // 1. Build + install OpenBabel into OUT_DIR.
    //
    // We force a Release build regardless of the cargo profile so the C++ side
    // always uses the release CRT (/MD on MSVC), matching what the cc/cxx-build
    // compiled shim uses — mixing CRTs would be an ABI hazard. It also avoids
    // needing debug builds of any optional dependency.
    //
    // Optional features that pull in heavy external deps (Boost/Eigen/Cairo/
    // RapidJSON) are disabled; the formats we need do not require them.
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
    let mut cfg = cmake::Config::new(&ob_src);
    cfg.profile("Release")
        .define("BUILD_GUI", "OFF")
        .define("ENABLE_TESTS", "OFF")
        .define("BUILD_SHARED", "ON")
        .define("WITH_MAEPARSER", "OFF")
        .define("WITH_COORDGEN", "OFF")
        .define("WITH_JSON", "OFF")
        .define("WITH_INCHI", "ON")
        .define("OB_USE_PREBUILT_BINARIES", "OFF")
        .define("OPENBABEL_USE_SYSTEM_INCHI", "OFF");

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
    //
    //     We copy every `.dll` (openbabel-3.dll plus its runtime dependencies
    //     such as the bundled `inchi` library, which `inchiformat.obf` links
    //     against) and every `.obf` plugin.
    if let Some(profile_dir) = out_dir.ancestors().nth(3) {
        let mut runtime = Vec::new();
        if let Ok(entries) = fs::read_dir(&bin_dir) {
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
                    copy_if_newer(src, &dest_dir.join(name));
                }
            }
        }
    }
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
