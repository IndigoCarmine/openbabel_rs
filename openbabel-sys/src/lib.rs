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

        /// Opaque owner of a compiled OpenBabel `OBSmartsPattern`.
        type Smarts;

        /// Opaque owner of a compiled `OBChemTsfm` (SMARTS→SMARTS transform).
        type Transform;

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
        fn atom_partial_charge(mol: &Molecule, idx: u32) -> f64;
        fn atom_is_aromatic(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_in_ring(mol: &Molecule, idx: u32) -> bool;
        fn atom_degree(mol: &Molecule, idx: u32) -> u32;
        fn atom_total_valence(mol: &Molecule, idx: u32) -> u32;
        fn atom_implicit_h_count(mol: &Molecule, idx: u32) -> u32;
        fn atom_hybridization(mol: &Molecule, idx: u32) -> u32;
        fn atom_is_hbond_donor(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_hbond_acceptor(mol: &Molecule, idx: u32) -> bool;

        // Bond accessors (idx is 0-based, 0..num_bonds). Atom indices returned
        // are 1-based.
        fn bond_begin_idx(mol: &Molecule, idx: u32) -> u32;
        fn bond_end_idx(mol: &Molecule, idx: u32) -> u32;
        fn bond_order(mol: &Molecule, idx: u32) -> u32;
        fn bond_is_aromatic(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_in_ring(mol: &Molecule, idx: u32) -> bool;

        // SMARTS substructure matching.
        /// Compile a SMARTS pattern; null on invalid syntax.
        fn smarts_new(pattern: &str) -> UniquePtr<Smarts>;
        /// Number of atoms in the pattern (= length of each match).
        fn smarts_atom_count(smarts: &Smarts) -> u32;
        /// Whether the pattern matches `mol` at least once.
        fn smarts_matches(smarts: &Smarts, mol: &Molecule) -> bool;
        /// Unique matches flattened to 1-based atom indices; reshape with
        /// `smarts_atom_count`.
        fn smarts_match_atoms(smarts: &Smarts, mol: &Molecule) -> Vec<u32>;

        // Fingerprints & similarity.
        /// Fingerprint of `mol` via plugin `id` ("FP2"/"FP3"/"FP4"/"MACCS");
        /// empty on unknown id.
        fn fingerprint(mol: &Molecule, id: &str) -> Vec<u32>;
        /// Tanimoto coefficient between two fingerprints.
        fn tanimoto(a: &[u32], b: &[u32]) -> f64;

        // Descriptors.
        /// Numeric descriptor `id` ("logP"/"TPSA"/"MR"/"MW"/...) of `mol`;
        /// sets `ok` to false for an unknown id.
        fn descriptor(mol: &Molecule, id: &str, ok: &mut bool) -> f64;

        // Force fields.
        /// Single-point energy under force field `ff_id`; sets `ok` to false on
        /// unknown field or setup failure.
        fn mol_energy(mol: &Molecule, ff_id: &str, ok: &mut bool) -> f64;
        /// Energy unit of `ff_id` (e.g. "kcal/mol"); empty if unknown.
        fn forcefield_unit(ff_id: &str) -> String;
        /// Minimize `mol` in place (`steps` conjugate-gradient steps); returns
        /// final energy, sets `ok` to false on failure.
        fn mol_optimize(mol: Pin<&mut Molecule>, ff_id: &str, steps: u32, ok: &mut bool) -> f64;

        // 3D generation.
        /// Generate 3D coordinates in place; `speed` ∈ {fastest,fast,med,slow,best}.
        fn mol_make_3d(mol: Pin<&mut Molecule>, speed: &str) -> bool;
        /// Coordinate dimension: 0, 2, or 3.
        fn mol_dimension(mol: &Molecule) -> u32;

        // Partial charges.
        /// Assign partial charges with model `model` ("gasteiger"/"mmff94"/
        /// "eem"/...); false on unknown model or failure.
        fn mol_compute_charges(mol: Pin<&mut Molecule>, model: &str) -> bool;

        // Structure alignment.
        /// Least-squares superpose `mol` onto `reference`, updating `mol`'s
        /// coordinates in place and returning the RMSD. `include_h` counts H in
        /// the fit; `symmetry` allows symmetry-equivalent remapping. Sets `ok`
        /// to false on atom-count mismatch or alignment failure.
        fn mol_align(
            mol: Pin<&mut Molecule>,
            reference: &Molecule,
            include_h: bool,
            symmetry: bool,
            ok: &mut bool,
        ) -> f64;

        // 2D depiction.
        /// Generate 2D coordinates in place via the `gen2D` op; false on failure.
        fn mol_make_2d(mol: Pin<&mut Molecule>) -> bool;
        /// Render `mol` to an SVG document. `all_carbons` labels every carbon;
        /// `atom_indices` annotates atoms with their index. 2D coordinates are
        /// generated automatically if absent. Sets `ok` to false on failure.
        fn mol_to_svg(
            mol: &Molecule,
            all_carbons: bool,
            atom_indices: bool,
            ok: &mut bool,
        ) -> String;

        // Stereochemistry.
        /// Force (re)perception of stereochemistry from structure.
        fn mol_perceive_stereo(mol: Pin<&mut Molecule>);
        /// Counts of perceived tetrahedral / cis-trans stereo units.
        fn mol_num_tetrahedral_stereo(mol: &Molecule) -> u32;
        fn mol_num_cistrans_stereo(mol: &Molecule) -> u32;
        /// Whether atom `idx` (1-based) is a tetrahedral stereocenter.
        fn atom_is_tetrahedral_stereo(mol: &Molecule, idx: u32) -> bool;
        /// Winding at atom `idx`: 1 = clockwise, 2 = anticlockwise, 0 = none.
        fn atom_tetrahedral_winding(mol: &Molecule, idx: u32) -> i32;
        /// Whether bond `idx` (0-based) is a cis/trans stereo unit.
        fn bond_is_cistrans_stereo(mol: &Molecule, idx: u32) -> bool;

        // Reaction / SMIRKS-like transforms.
        /// Compile a reactant→product SMARTS transform; null if either is invalid.
        fn transform_new(reactant: &str, product: &str) -> UniquePtr<Transform>;
        /// Apply the transform to every match in `mol`, editing it in place;
        /// false if nothing matched or it failed.
        fn transform_apply(t: &Transform, mol: Pin<&mut Molecule>) -> bool;

        // Conformer search.
        /// Genetic-algorithm conformer search targeting `count` conformers,
        /// stored in `mol` (which must already have a 3D structure). Returns the
        /// number of conformers now stored.
        fn mol_generate_conformers(mol: Pin<&mut Molecule>, count: u32) -> u32;
        /// Number of stored conformers.
        fn mol_num_conformers(mol: &Molecule) -> u32;
        /// Make conformer `index` the active coordinates (no-op if out of range).
        fn mol_set_conformer(mol: Pin<&mut Molecule>, index: u32);
    }
}

/// Absolute paths, baked at build time, to the OpenBabel runtime that this
/// crate was linked against.
pub mod paths {
    include!(concat!(env!("OUT_DIR"), "/paths.rs"));
}
