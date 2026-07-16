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

        /// Opaque owner of an `OBReaction` (reactant / product / agent lists).
        type Reaction;

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
        /// Flat `[x,y,z,…]` of conformer `index` without changing the active one.
        fn mol_conformer_coordinates(mol: &Molecule, index: u32) -> Vec<f64>;
        /// Energy of each conformer under `ff_id`; restores the active conformer.
        fn mol_conformer_energies(mol: &Molecule, ff_id: &str) -> Vec<f64>;

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

        // File I/O (empty `format` = auto-detect from the file extension).
        fn mol_read_file(path: &str, format: &str) -> UniquePtr<Molecule>;
        fn mol_read_file_many(path: &str, format: &str) -> UniquePtr<CxxVector<Molecule>>;
        fn mol_write_file(mol: &Molecule, path: &str, format: &str, ok: &mut bool);

        // Ring access (SSSR; ring_idx is 0-based, 0..mol_num_rings).
        fn ring_size(mol: &Molecule, ring_idx: u32) -> u32;
        fn ring_atom_indices(mol: &Molecule, ring_idx: u32) -> Vec<u32>;
        fn ring_is_aromatic(mol: &Molecule, ring_idx: u32) -> bool;

        // Graph navigation (atom idx 1-based; returned indices 0-based).
        fn atom_neighbor_indices(mol: &Molecule, idx: u32) -> Vec<u32>;
        fn atom_bond_indices(mol: &Molecule, idx: u32) -> Vec<u32>;
        fn atom_count_bonds_of_order(mol: &Molecule, idx: u32, order: u32) -> u32;
        fn atom_explicit_h_count(mol: &Molecule, idx: u32) -> u32;
        /// 0-based bond index joining 0-based atoms `a`/`b`; -1 if not bonded.
        fn mol_bond_between(mol: &Molecule, a: u32, b: u32) -> i32;
        /// 0-based atom across bond `bond_idx` from 0-based `atom_idx`; -1 if not on it.
        fn bond_other_atom(mol: &Molecule, bond_idx: u32, atom_idx: u32) -> i32;

        // Crystallography (unit cell).
        fn mol_has_unit_cell(mol: &Molecule) -> bool;
        /// `[a, b, c, alpha, beta, gamma]`; empty if no unit cell.
        fn mol_cell_parameters(mol: &Molecule) -> Vec<f64>;
        fn mol_cell_volume(mol: &Molecule) -> f64;
        fn mol_cell_spacegroup(mol: &Molecule) -> String;
        fn mol_cell_lattice_type(mol: &Molecule) -> u32;
        fn mol_cell_to_fractional(mol: &Molecule, x: f64, y: f64, z: f64) -> Vec<f64>;
        fn mol_cell_to_cartesian(mol: &Molecule, x: f64, y: f64, z: f64) -> Vec<f64>;

        // Symmetry & canonical ordering (one value per atom, atom order).
        fn mol_symmetry_classes(mol: &Molecule) -> Vec<u32>;
        fn mol_canonical_ranks(mol: &Molecule) -> Vec<u32>;

        // Reactions (formats "rxn" = MDL RXN, "rsmi" = reaction SMILES).
        fn reaction_new() -> UniquePtr<Reaction>;
        fn reaction_read(format: &str, data: &str) -> UniquePtr<Reaction>;
        fn reaction_write(r: &Reaction, format: &str, ok: &mut bool) -> String;
        fn reaction_num_reactants(r: &Reaction) -> u32;
        fn reaction_num_products(r: &Reaction) -> u32;
        fn reaction_num_agents(r: &Reaction) -> u32;
        fn reaction_reactant(r: &Reaction, i: u32) -> UniquePtr<Molecule>;
        fn reaction_product(r: &Reaction, i: u32) -> UniquePtr<Molecule>;
        fn reaction_agent(r: &Reaction, i: u32) -> UniquePtr<Molecule>;
        fn reaction_add_reactant(r: Pin<&mut Reaction>, mol: &Molecule);
        fn reaction_add_product(r: Pin<&mut Reaction>, mol: &Molecule);
        fn reaction_add_agent(r: Pin<&mut Reaction>, mol: &Molecule);
        fn reaction_title(r: &Reaction) -> String;
        fn reaction_set_title(r: Pin<&mut Reaction>, title: &str);
        fn reaction_comment(r: &Reaction) -> String;
        fn reaction_set_comment(r: Pin<&mut Reaction>, comment: &str);
        fn reaction_is_reversible(r: &Reaction) -> bool;
        fn reaction_set_reversible(r: Pin<&mut Reaction>, value: bool);

        // Subgraph isomorphism & automorphisms (flat, `width` atoms per row).
        fn mol_substructure_mappings(
            query: &Molecule,
            target: &Molecule,
            width: &mut u32,
        ) -> Vec<u32>;
        fn mol_automorphisms(mol: &Molecule, width: &mut u32) -> Vec<u32>;

        // Geometry & topology (niche; atom indices 0-based).
        fn mol_set_torsion(mol: Pin<&mut Molecule>, a: u32, b: u32, c: u32, d: u32, radians: f64);
        fn mol_find_children(mol: &Molecule, from: u32, to: u32) -> Vec<u32>;
        fn mol_largest_fragment(mol: &Molecule) -> Vec<u32>;
        fn mol_set_total_charge(mol: Pin<&mut Molecule>, charge: i32);
        fn mol_set_total_spin(mol: Pin<&mut Molecule>, spin: u32);

        // Per-atom / per-bond string data (atom idx 1-based, bond idx 0-based).
        fn atom_set_data(mol: Pin<&mut Molecule>, idx: u32, key: &str, value: &str);
        fn atom_get_data(mol: &Molecule, idx: u32, key: &str, ok: &mut bool) -> String;
        fn bond_set_data(mol: Pin<&mut Molecule>, idx: u32, key: &str, value: &str);
        fn bond_get_data(mol: &Molecule, idx: u32, key: &str, ok: &mut bool) -> String;

        // Inter-atom distance (atoms 1-based) & 2D wedge/hash bond stereo.
        fn mol_distance(mol: &Molecule, i: u32, j: u32) -> f64;
        fn bond_is_wedge(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_hash(mol: &Molecule, idx: u32) -> bool;
        fn bond_set_wedge(mol: Pin<&mut Molecule>, idx: u32, value: bool);
        fn bond_set_hash(mol: Pin<&mut Molecule>, idx: u32, value: bool);

        // Persistent atom ids, connectivity relations, LSSR (atoms 1-based).
        fn atom_id(mol: &Molecule, idx: u32) -> u64;
        fn atom_set_id(mol: Pin<&mut Molecule>, idx: u32, id: u64);
        fn atom_is_connected(mol: &Molecule, a: u32, b: u32) -> bool;
        fn atom_is_one_three(mol: &Molecule, a: u32, b: u32) -> bool;
        fn atom_is_one_four(mol: &Molecule, a: u32, b: u32) -> bool;
        fn mol_lssr_sizes(mol: &Molecule) -> Vec<u32>;

        // Perception state flags & targeted hydrogen editing (atom idx 1-based).
        fn mol_has_aromatic_perceived(mol: &Molecule) -> bool;
        fn mol_has_sssr_perceived(mol: &Molecule) -> bool;
        fn mol_has_ring_atoms_perceived(mol: &Molecule) -> bool;
        fn mol_has_chains_perceived(mol: &Molecule) -> bool;
        fn mol_has_hydrogens_added(mol: &Molecule) -> bool;
        fn mol_has_nonzero_coords(mol: &Molecule) -> bool;
        fn mol_add_hydrogens_to_atom(mol: Pin<&mut Molecule>, idx: u32) -> bool;
        fn mol_delete_hydrogens_of_atom(mol: Pin<&mut Molecule>, idx: u32) -> bool;

        // Structured torsion / angle data: flat 0-based atom indices (3 per
        // angle [vertex,a,b]; 4 per torsion [a,b,c,d]).
        fn mol_find_angles(mol: &Molecule) -> Vec<u32>;
        fn mol_find_torsions(mol: &Molecule) -> Vec<u32>;

        // Perception-state flag setters (the setter side of the Has* queries).
        fn mol_set_aromatic_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_sssr_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_ring_atoms_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_chains_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_hydrogens_added(mol: Pin<&mut Molecule>, value: bool);

        // Axial / equatorial ring position (atom idx 1-based).
        fn atom_is_axial(mol: &Molecule, idx: u32) -> bool;

        // Remaining perception-state flags (readers + setters).
        fn mol_has_lssr_perceived(mol: &Molecule) -> bool;
        fn mol_has_atom_types_perceived(mol: &Molecule) -> bool;
        fn mol_has_ring_types_perceived(mol: &Molecule) -> bool;
        fn mol_has_chirality_perceived(mol: &Molecule) -> bool;
        fn mol_has_partial_charges_perceived(mol: &Molecule) -> bool;
        fn mol_has_hybridization_perceived(mol: &Molecule) -> bool;
        fn mol_has_closure_bonds_perceived(mol: &Molecule) -> bool;
        fn mol_is_corrected_for_ph(mol: &Molecule) -> bool;
        fn mol_has_spin_multiplicity_assigned(mol: &Molecule) -> bool;
        fn mol_set_lssr_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_atom_types_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_ring_types_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_chirality_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_partial_charges_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_hybridization_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_closure_bonds_perceived(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_corrected_for_ph(mol: Pin<&mut Molecule>, value: bool);
        fn mol_set_spin_multiplicity_assigned(mol: Pin<&mut Molecule>, value: bool);

        // Atom functional-group / environment predicates (atom idx 1-based).
        fn atom_is_carboxyl_oxygen(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_phosphate_oxygen(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_sulfate_oxygen(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_nitro_oxygen(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_amide_nitrogen(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_aromatic_noxide(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_nonpolar_hydrogen(mol: &Molecule, idx: u32) -> bool;
        fn atom_is_hbond_donor_h(mol: &Molecule, idx: u32) -> bool;
        fn atom_count_free_oxygens(mol: &Molecule, idx: u32) -> u32;
        fn atom_count_free_sulfurs(mol: &Molecule, idx: u32) -> u32;
        fn atom_count_ring_bonds(mol: &Molecule, idx: u32) -> u32;
        fn atom_smallest_bond_angle(mol: &Molecule, idx: u32) -> f64;
        fn atom_average_bond_angle(mol: &Molecule, idx: u32) -> f64;
        fn atom_lewis_acid_base_counts(mol: &Molecule, idx: u32, acid: &mut i32, base: &mut i32);

        // Bond classification predicates (bond idx 0-based).
        fn bond_is_primary_amide(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_secondary_amide(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_tertiary_amide(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_wedge_or_hash(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_cis_or_trans(mol: &Molecule, idx: u32) -> bool;
        fn bond_is_double_bond_geometry(mol: &Molecule, idx: u32) -> bool;
    }
}

/// Absolute paths, baked at build time, to the OpenBabel runtime that this
/// crate was linked against.
pub mod paths {
    include!(concat!(env!("OUT_DIR"), "/paths.rs"));
}
