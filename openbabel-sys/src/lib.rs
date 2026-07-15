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

        /// Opaque owner of an `OBFFConstraints` set (fixed atoms, distance /
        /// angle / torsion restraints, ignored atoms).
        type Constraints;

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

        // Element data (by atomic number).
        fn element_symbol(atomic_number: u32) -> String;
        fn element_name(atomic_number: u32) -> String;
        fn element_atomic_number(symbol: &str) -> u32;
        fn element_mass(atomic_number: u32) -> f64;
        fn element_exact_mass(atomic_number: u32) -> f64;
        fn element_electronegativity(atomic_number: u32) -> f64;
        fn element_covalent_radius(atomic_number: u32) -> f64;
        fn element_vdw_radius(atomic_number: u32) -> f64;
        fn element_max_bonds(atomic_number: u32) -> u32;

        // More atom accessors.
        fn atom_type(mol: &Molecule, idx: u32) -> String;
        fn atom_isotope(mol: &Molecule, idx: u32) -> u32;
        fn atom_atomic_mass(mol: &Molecule, idx: u32) -> f64;
        fn atom_exact_mass(mol: &Molecule, idx: u32) -> f64;
        fn atom_spin_multiplicity(mol: &Molecule, idx: u32) -> i32;
        fn atom_heavy_degree(mol: &Molecule, idx: u32) -> u32;
        fn atom_hetero_degree(mol: &Molecule, idx: u32) -> u32;
        fn atom_is_chiral(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_heteroatom(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_metal(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_polar_hydrogen(mol: &Molecule, idx: u32) -> bool;
        fn atom_member_of_ring_count(mol: &Molecule, idx: u32) -> u32;
        fn atom_member_of_ring_size(mol: &Molecule, idx: u32) -> u32;
        fn atom_is_in_ring_size(mol: &Molecule, idx: u32, size: u32) -> bool;

        // More bond accessors.
        fn bond_length(mol: &Molecule, idx: u32) -> f64;
        fn bond_equilibrium_length(mol: &Molecule, idx: u32) -> f64;
        fn bond_is_rotor(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_amide(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_ester(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_carbonyl(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_closure(mol: &Molecule, idx: u32) -> bool;

        // More whole-molecule methods.
        fn mol_num_heavy_atoms(mol: &Molecule) -> u32;
        fn mol_num_rotors(mol: &Molecule) -> u32;
        fn mol_num_rings(mol: &Molecule) -> u32;
        fn mol_spaced_formula(mol: &Molecule) -> String;
        fn mol_spin_multiplicity(mol: &Molecule) -> u32;
        /// Translate the molecule so its centroid is at the origin.
        fn mol_center(mol: Pin<&mut Molecule>);
        /// Valence angle i-j-k in degrees (1-based atom indices).
        fn mol_angle(mol: &Molecule, i: u32, j: u32, k: u32) -> f64;
        /// Torsion angle i-j-k-l in degrees (1-based atom indices).
        fn mol_torsion(mol: &Molecule, i: u32, j: u32, k: u32, l: u32) -> f64;
        /// Deep-copy the molecule.
        fn mol_clone(mol: &Molecule) -> UniquePtr<Molecule>;
        /// Remove fragments smaller than `threshold` atoms (0 = keep largest).
        fn mol_strip_salts(mol: Pin<&mut Molecule>, threshold: u32) -> bool;
        /// Split into disconnected fragments (each a fresh molecule).
        fn mol_separate(mol: &Molecule) -> UniquePtr<CxxVector<Molecule>>;
        /// Set/get a string property (OBPairData) by key.
        fn mol_set_property(mol: Pin<&mut Molecule>, key: &str, value: &str);
        fn mol_get_property(mol: &Molecule, key: &str, ok: &mut bool) -> String;

        // Residues (biopolymer/PDB substructure; `res_idx` is 0-based).
        /// Number of perceived/loaded residues.
        fn mol_num_residues(mol: &Molecule) -> u32;
        /// Residue name (e.g. "GLY", "HOH"); empty if `res_idx` out of range.
        fn residue_name(mol: &Molecule, res_idx: u32) -> String;
        /// Residue sequence number (PDB `resSeq`).
        fn residue_number(mol: &Molecule, res_idx: u32) -> i32;
        /// Residue sequence number as text (may include insertion code).
        fn residue_number_string(mol: &Molecule, res_idx: u32) -> String;
        /// Chain identifier as a string ("A", "B", …); empty if unset.
        fn residue_chain(mol: &Molecule, res_idx: u32) -> String;
        /// Insertion code as a string; empty if unset.
        fn residue_insertion_code(mol: &Molecule, res_idx: u32) -> String;
        /// Atom counts for the residue.
        fn residue_num_atoms(mol: &Molecule, res_idx: u32) -> u32;
        fn residue_num_heavy_atoms(mol: &Molecule, res_idx: u32) -> u32;
        /// 0-based atom indices belonging to the residue.
        fn residue_atom_indices(mol: &Molecule, res_idx: u32) -> Vec<u32>;
        /// 0-based residue index for atom `idx` (1-based); -1 if none.
        fn atom_residue_index(mol: &Molecule, idx: u32) -> i32;
        /// PDB atom name for atom `idx` within its residue (e.g. " CA "); empty
        /// if the atom has no residue.
        fn atom_residue_atom_id(mol: &Molecule, idx: u32) -> String;
        /// Whether atom `idx` is a HETATM in its residue.
        fn atom_is_hetatm(mol: &Molecule, idx: u32) -> bool;
        /// PDB serial number of atom `idx` within its residue; 0 if none.
        fn atom_serial_number(mol: &Molecule, idx: u32) -> u32;

        // Spectra.
        /// Spectrophore descriptor (default 48 values; needs 3D coordinates).
        /// Empty on failure.
        fn mol_spectrophore(mol: &Molecule) -> Vec<f64>;
        /// Vibrational frequencies (cm⁻¹); empty unless the molecule carries
        /// `OBVibrationData` (e.g. read from a comp-chem output).
        fn mol_vibration_frequencies(mol: &Molecule) -> Vec<f64>;
        /// Vibrational IR intensities (km/mol); empty unless present.
        fn mol_vibration_intensities(mol: &Molecule) -> Vec<f64>;

        /// All atom coordinates flattened as `[x0,y0,z0, x1,y1,z1, …]`.
        fn mol_coordinates(mol: &Molecule) -> Vec<f64>;

        // Force-field constraints (atom indices 0-based).
        /// Create an empty constraint set.
        fn constraints_new() -> UniquePtr<Constraints>;
        fn constraints_add_ignore(c: Pin<&mut Constraints>, atom: u32);
        fn constraints_add_atom(c: Pin<&mut Constraints>, atom: u32);
        fn constraints_add_atom_x(c: Pin<&mut Constraints>, atom: u32);
        fn constraints_add_atom_y(c: Pin<&mut Constraints>, atom: u32);
        fn constraints_add_atom_z(c: Pin<&mut Constraints>, atom: u32);
        fn constraints_add_distance(c: Pin<&mut Constraints>, a: u32, b: u32, length: f64);
        fn constraints_add_angle(c: Pin<&mut Constraints>, a: u32, b: u32, d: u32, angle: f64);
        fn constraints_add_torsion(
            c: Pin<&mut Constraints>,
            a: u32,
            b: u32,
            d: u32,
            e: u32,
            torsion: f64,
        );
        fn constraints_set_factor(c: Pin<&mut Constraints>, factor: f64);

        // Geometry optimization. Each call runs the whole minimization atomically
        // (the safe wrapper holds the global lock for its full duration) so the
        // shared static constraint state can't be corrupted by a concurrent run.
        // `algorithm` is 0 = steepest descent, 1 = conjugate gradients, 2 = L-BFGS.
        /// Run to completion; write final coordinates to `mol` and return the
        /// final energy. Sets `ok` false if the force field is unknown / setup fails.
        fn optimizer_run_to_end(
            mol: Pin<&mut Molecule>,
            ff_id: &str,
            algorithm: u32,
            steps: u32,
            econv: f64,
            constraints: &Constraints,
            ok: &mut bool,
        ) -> f64;
        /// Run to completion recording a frame every `frame_interval` steps;
        /// flattened as repeated `[energy, x0,y0,z0, …]` (`1 + 3·num_atoms` each).
        /// `mol` ends at the final geometry. Empty on unknown FF / setup failure.
        fn optimizer_run_trajectory(
            mol: Pin<&mut Molecule>,
            ff_id: &str,
            algorithm: u32,
            steps: u32,
            econv: f64,
            constraints: &Constraints,
            frame_interval: u32,
        ) -> Vec<f64>;

        // Molecule construction & editing.
        /// Add an atom of `atomic_num`; returns its 0-based index.
        fn mol_add_atom(mol: Pin<&mut Molecule>, atomic_num: u32) -> u32;
        /// Bond 0-based atoms `begin`–`end` with `order`; false if out of range.
        fn mol_add_bond(mol: Pin<&mut Molecule>, begin: u32, end: u32, order: u32) -> bool;
        fn mol_delete_atom(mol: Pin<&mut Molecule>, idx: u32) -> bool;
        fn mol_delete_bond(mol: Pin<&mut Molecule>, idx: u32) -> bool;
        fn mol_begin_modify(mol: Pin<&mut Molecule>);
        fn mol_end_modify(mol: Pin<&mut Molecule>);
        fn mol_clear(mol: Pin<&mut Molecule>);
        fn mol_translate(mol: Pin<&mut Molecule>, x: f64, y: f64, z: f64);
        /// Overwrite all coordinates; false unless `coords.len() == 3·num_atoms`.
        fn mol_set_coordinates(mol: Pin<&mut Molecule>, coords: &[f64]) -> bool;
        fn mol_set_dimension(mol: Pin<&mut Molecule>, dim: u32);
        fn mol_connect_the_dots(mol: Pin<&mut Molecule>);
        fn mol_perceive_bond_orders(mol: Pin<&mut Molecule>);
        fn mol_add_polar_hydrogens(mol: Pin<&mut Molecule>) -> bool;
        fn mol_convert_dative_bonds(mol: Pin<&mut Molecule>) -> bool;
        fn mol_assign_spin_multiplicity(mol: Pin<&mut Molecule>) -> bool;
        fn mol_add_hydrogens_ph(mol: Pin<&mut Molecule>, ph: f64) -> bool;

        // Atom setters (idx is 1-based).
        fn atom_set_atomic_num(mol: Pin<&mut Molecule>, idx: u32, atomic_num: u32);
        fn atom_set_formal_charge(mol: Pin<&mut Molecule>, idx: u32, charge: i32);
        fn atom_set_position(mol: Pin<&mut Molecule>, idx: u32, x: f64, y: f64, z: f64);
        fn atom_set_isotope(mol: Pin<&mut Molecule>, idx: u32, isotope: u32);
        fn atom_set_spin_multiplicity(mol: Pin<&mut Molecule>, idx: u32, spin: i32);
        fn atom_set_partial_charge(mol: Pin<&mut Molecule>, idx: u32, charge: f64);
        fn atom_set_type(mol: Pin<&mut Molecule>, idx: u32, type_name: &str);
        fn atom_set_implicit_h(mol: Pin<&mut Molecule>, idx: u32, count: u32);

        // Bond setters (idx is 0-based).
        fn bond_set_order(mol: Pin<&mut Molecule>, idx: u32, order: u32);
        fn bond_set_length(mol: Pin<&mut Molecule>, idx: u32, length: f64) -> bool;

        /// Read EVERY record from `data` in `format` (multi-record SDF, one
        /// SMILES per line, …).
        fn mol_read_many(format: &str, data: &str) -> UniquePtr<CxxVector<Molecule>>;

        // Ring access (SSSR; ring_idx is 0-based, 0..mol_num_rings).
        fn ring_size(mol: &Molecule, ring_idx: u32) -> u32;
        fn ring_atom_indices(mol: &Molecule, ring_idx: u32) -> Vec<u32>;
        fn ring_is_aromatic(mol: &Molecule, ring_idx: u32) -> bool;
    }
}

/// Absolute paths, baked at build time, to the OpenBabel runtime that this
/// crate was linked against.
pub mod paths {
    include!(concat!(env!("OUT_DIR"), "/paths.rs"));
}
