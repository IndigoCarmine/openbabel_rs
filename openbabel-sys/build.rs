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
    // Eigen (header-only) is vendored as a submodule. Pointing OpenBabel's
    // `find_package(Eigen3)` at it enables HAVE_EIGEN3, which compiles OBAlign
    // (structure superposition) and unlocks distance-geometry 3D generation.
    let eigen_dir = workspace_root.join("vendor").join("eigen");

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
    assert!(
        eigen_dir.join("Eigen").join("Core").exists(),
        "Eigen headers not found at {}.\n\
         Run: git submodule update --init --recursive",
        eigen_dir.display()
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
        .define("OPENBABEL_USE_SYSTEM_INCHI", "OFF")
        .define("EIGEN3_INCLUDE_DIR", &eigen_include);

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
                        copy_if_newer(src, &dest_dir.join(name));
                    }
                }
            }
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
