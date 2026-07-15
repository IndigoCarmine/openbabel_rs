//! Low-level FFI bridge to OpenBabel via [`cxx`].
//!
//! Everything in the bridge mirrors, one-to-one, the C++ shim in
//! `shim/shim.h` / `shim/shim.cc`. cxx checks the signatures match at compile
//! time, so a mismatch is a build error rather than UB.
//!
//! cxx type mapping reminders:
//!   * `&Molecule`         <-> `const Molecule&`
//!   * `Pin<&mut Molecule>` <-> `Molecule&`   (needed for mutating calls)
//!   * `&str`              <-> `rust::Str`
//!   * `String`            <-> `rust::String`
//!
//! Prefer the safe `openbabel` crate; this crate is the unsafe plumbing.

#[cxx::bridge(namespace = "ob_shim")]
pub mod ffi {
    unsafe extern "C++" {
        include!("shim.h");

        /// Opaque owner of an OpenBabel `OBMol`.
        type Molecule;

        /// OpenBabel release version, e.g. `"3.2.1"`.
        fn release_version() -> String;

        /// Set an environment variable through the C runtime so OpenBabel's
        /// `getenv` lookups (e.g. `BABEL_LIBDIR`) observe it.
        fn set_env(key: &str, value: &str);

        /// Create an empty molecule.
        fn mol_new() -> UniquePtr<Molecule>;

        /// Parse `data` as `format` ("smi", "mol", "sdf", "pdb", ...).
        /// Returns a null pointer on unknown format or parse failure.
        fn mol_read(format: &str, data: &str) -> UniquePtr<Molecule>;

        /// Serialize `mol` to `format`. Sets `ok` to false on unknown format.
        fn mol_write(mol: &Molecule, format: &str, ok: &mut bool) -> String;

        // Whole-molecule properties.
        fn mol_formula(mol: &Molecule) -> String;
        fn mol_mol_wt(mol: &Molecule) -> f64;
        fn mol_exact_mass(mol: &Molecule) -> f64;
        fn mol_total_charge(mol: &Molecule) -> i32;
        fn mol_num_atoms(mol: &Molecule) -> u32;
        fn mol_num_bonds(mol: &Molecule) -> u32;
        fn mol_title(mol: &Molecule) -> String;
        fn mol_set_title(mol: Pin<&mut Molecule>, title: &str);
        fn mol_add_hydrogens(mol: Pin<&mut Molecule>);
        fn mol_delete_hydrogens(mol: Pin<&mut Molecule>);

        // Atom accessors (idx is 1-based, 1..=num_atoms).
        fn atom_atomic_num(mol: &Molecule, idx: u32) -> u32;
        fn atom_x(mol: &Molecule, idx: u32) -> f64;
        fn atom_y(mol: &Molecule, idx: u32) -> f64;
        fn atom_z(mol: &Molecule, idx: u32) -> f64;
        fn atom_formal_charge(mol: &Molecule, idx: u32) -> i32;
        fn atom_is_aromatic(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_in_ring(mol: &Molecule, idx: u32) -> bool;

        // Bond accessors (idx is 0-based, 0..num_bonds). Atom indices returned
        // are 1-based.
        fn bond_begin_idx(mol: &Molecule, idx: u32) -> u32;
        fn bond_end_idx(mol: &Molecule, idx: u32) -> u32;
        fn bond_order(mol: &Molecule, idx: u32) -> u32;
        fn bond_is_aromatic(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_in_ring(mol: &Molecule, idx: u32) -> bool;
    }
}

/// Absolute paths, baked at build time, to the OpenBabel runtime that this
/// crate was linked against.
pub mod paths {
    include!(concat!(env!("OUT_DIR"), "/paths.rs"));
}
