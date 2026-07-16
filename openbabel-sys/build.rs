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
use std::process::Command;

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

    // Fail early with an actionable message if cmake is missing. The `cmake`
    // crate shells out to the `cmake` binary and, when it is absent, panics
    // deep inside the crate with a generic "is `cmake` not installed?" — which
    // buries the one thing the user needs: the command to install it. Cargo has
    // no way to install system build tools itself, so we can only guide.
    ensure_cmake_present();

    // Make the bundled InChI library linkable on MSVC (see fn docs).
    ensure_inchi_auxinfo_stubs(&ob_src);

    // Add read-only accessors to the force-field headers so the Rust term
    // exporter can read their precomputed calculation vectors (see fn docs).
    ensure_ff_accessors(&ob_src);

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
        // OpenBabel's concrete force-field headers (e.g. forcefielduff.h) live
        // in the source tree, not the installed headers. The Rust term exporter
        // (ff_export_terms) includes one to read a force field's precomputed
        // calculation vectors; put that directory on the shim's include path.
        .include(ob_src.join("src").join("forcefields"))
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
