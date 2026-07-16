// Thin C++ shim over OpenBabel, exposing a cxx-friendly surface.
//
// Rules for this file:
//   * Every function mirrors a declaration in `../src/lib.rs`. Keep them in
//     sync — cxx verifies the signatures match at build time.
//   * No C++ exception may escape across the FFI boundary. Every body is
//     wrapped so OpenBabel errors become the documented failure value.
//   * `Molecule` is defined completely here (not just forward-declared): cxx's
//     generated glue for `UniquePtr<Molecule>` instantiates the destructor and
//     therefore needs the complete type visible from this header.
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/phmodel.h>
#include <openbabel/reaction.h>

#include "rust/cxx.h"

namespace ob_shim {

// Opaque (to Rust) wrapper owning an OpenBabel OBMol. Owning by value keeps
// lifetimes trivial: dropping the Rust UniquePtr frees everything.
struct Molecule {
  OpenBabel::OBMol mol;
};

// Opaque (to Rust) wrapper owning a compiled SMARTS pattern. Defined here (not
// forward-declared) for the same reason as Molecule: cxx's UniquePtr glue needs
// the complete type.
struct Smarts {
  OpenBabel::OBSmartsPattern pat;
};

// Opaque (to Rust) wrapper owning a compiled SMARTS→SMARTS transformation
// (OBChemTsfm), applied to a molecule to edit it in place (SMIRKS-like).
struct Transform {
  OpenBabel::OBChemTsfm tsfm;
};

// Opaque (to Rust) set of force-field constraints (fixed atoms, distance /
// angle / torsion restraints, ignored atoms). Built incrementally, then handed
// to an optimizer_run_* call. Complete type (not forward-declared) for the
// UniquePtr glue.
struct Constraints {
  OpenBabel::OBFFConstraints c;
};

// Opaque (to Rust) wrapper owning a chemical reaction (OBReaction): lists of
// reactant / product / agent molecules plus title/comment/reversible flags.
// Complete type (not forward-declared) for the UniquePtr glue.
struct Reaction {
  OpenBabel::OBReaction rxn;
};

// OpenBabel release version string, e.g. "3.2.1".
rust::String release_version();

// Set an environment variable via the C runtime, so that OpenBabel's
// getenv()-based lookups (BABEL_LIBDIR / BABEL_DATADIR) see it. Rust's
// std::env::set_var uses the Win32 environment block, which the MSVC CRT's
// getenv does not read; going through the CRT here keeps them in sync.
void set_env(rust::Str key, rust::Str value);

// Create an empty molecule.
std::unique_ptr<Molecule> mol_new();

// Parse `data` as `format` (an OpenBabel format id such as "smi", "mol",
// "sdf", "pdb"). Returns a null pointer if the format is unknown or the data
// fails to parse.
std::unique_ptr<Molecule> mol_read(rust::Str format, rust::Str data);

// Serialize `mol` to `format`. `ok` is set to false if the format is unknown
// (returned string then empty); a valid empty result keeps `ok` true.
rust::String mol_write(const Molecule &mol, rust::Str format, bool &ok);

// --- Whole-molecule properties -------------------------------------------
rust::String mol_formula(const Molecule &mol);
double mol_mol_wt(const Molecule &mol);
double mol_exact_mass(const Molecule &mol);
int mol_total_charge(const Molecule &mol);
uint32_t mol_num_atoms(const Molecule &mol);
uint32_t mol_num_bonds(const Molecule &mol);
rust::String mol_title(const Molecule &mol);
void mol_set_title(Molecule &mol, rust::Str title);
void mol_add_hydrogens(Molecule &mol);
void mol_delete_hydrogens(Molecule &mol);

// --- Atom accessors (idx is 1-based, 1..=num_atoms, as in OpenBabel) ------
uint32_t atom_atomic_num(const Molecule &mol, uint32_t idx);
double atom_x(const Molecule &mol, uint32_t idx);
double atom_y(const Molecule &mol, uint32_t idx);
double atom_z(const Molecule &mol, uint32_t idx);
int atom_formal_charge(const Molecule &mol, uint32_t idx);
double atom_partial_charge(const Molecule &mol, uint32_t idx);
bool atom_is_aromatic(const Molecule &mol, uint32_t idx);
bool atom_is_in_ring(const Molecule &mol, uint32_t idx);
uint32_t atom_degree(const Molecule &mol, uint32_t idx);       // explicit connections
uint32_t atom_total_valence(const Molecule &mol, uint32_t idx);
uint32_t atom_implicit_h_count(const Molecule &mol, uint32_t idx);
uint32_t atom_hybridization(const Molecule &mol, uint32_t idx); // 1=sp,2=sp2,3=sp3,...
bool atom_is_hbond_donor(const Molecule &mol, uint32_t idx);
bool atom_is_hbond_acceptor(const Molecule &mol, uint32_t idx);

// --- Bond accessors (idx is 0-based, 0..num_bonds, as in OpenBabel) -------
// Returned atom indices are 1-based (OpenBabel's GetBeginAtomIdx/EndAtomIdx).
uint32_t bond_begin_idx(const Molecule &mol, uint32_t idx);
uint32_t bond_end_idx(const Molecule &mol, uint32_t idx);
uint32_t bond_order(const Molecule &mol, uint32_t idx);
bool bond_is_aromatic(const Molecule &mol, uint32_t idx);
bool bond_is_in_ring(const Molecule &mol, uint32_t idx);

// --- SMARTS substructure matching ----------------------------------------
// Compile a SMARTS pattern. Returns null if the pattern is invalid.
std::unique_ptr<Smarts> smarts_new(rust::Str pattern);
// Number of atoms in the pattern (the length of each match).
uint32_t smarts_atom_count(const Smarts &smarts);
// Whether the pattern matches `mol` at least once.
bool smarts_matches(const Smarts &smarts, const Molecule &mol);
// Unique matches, flattened: length is num_matches * atom_count, and atom
// indices are 1-based. Reshape on the Rust side using atom_count.
rust::Vec<uint32_t> smarts_match_atoms(const Smarts &smarts, const Molecule &mol);

// --- Fingerprints & similarity -------------------------------------------
// Compute the fingerprint of `mol` using the fingerprint plugin `id`
// ("FP2", "FP3", "FP4", "MACCS"). Returns an empty vector on unknown id.
rust::Vec<uint32_t> fingerprint(const Molecule &mol, rust::Str id);
// Tanimoto coefficient between two fingerprints.
double tanimoto(rust::Slice<const uint32_t> a, rust::Slice<const uint32_t> b);

// --- Descriptors ----------------------------------------------------------
// Evaluate a numeric descriptor plugin `id` ("logP", "TPSA", "MR", "MW", ...)
// on `mol`. `ok` is set to false for an unknown id.
double descriptor(const Molecule &mol, rust::Str id, bool &ok);

// --- Force fields ---------------------------------------------------------
// Single-point energy of `mol` under force field `ff_id` ("MMFF94", "UFF",
// "GAFF", "Ghemical", ...). NaN with ok=false if the field is unknown or setup
// fails. Meaningful only when `mol` has 3D coordinates.
double mol_energy(const Molecule &mol, rust::Str ff_id, bool &ok);
// The energy unit reported by force field `ff_id` (e.g. "kcal/mol"); empty for
// an unknown field.
rust::String forcefield_unit(rust::Str ff_id);
// Minimize `mol`'s geometry in place: `steps` conjugate-gradient steps under
// `ff_id`. Returns the final energy (NaN with ok=false on failure).
double mol_optimize(Molecule &mol, rust::Str ff_id, uint32_t steps, bool &ok);

// --- 3D structure generation ----------------------------------------------
// Generate 3D coordinates in place, like `obabel --gen3d`. `speed` is one of
// "fastest"/"fast"/"med"/"slow"/"best" (or a digit "1".."5"). Returns false on
// failure.
bool mol_make_3d(Molecule &mol, rust::Str speed);
// Coordinate dimension of `mol`: 0 (none), 2, or 3.
uint32_t mol_dimension(const Molecule &mol);

// --- Partial charges ------------------------------------------------------
// Assign partial atomic charges using charge model `model` ("gasteiger",
// "mmff94", "eem", "eqeq", "qeq", "qtpie"). Returns false on unknown model or
// failure. After this, atom_partial_charge() reflects the assigned charges.
bool mol_compute_charges(Molecule &mol, rust::Str model);

// --- Structure alignment (OBAlign) ----------------------------------------
// Least-squares superpose `mol` onto `reference` (Kabsch algorithm). Updates
// `mol`'s coordinates in place to the aligned pose and returns the RMSD.
// `include_h` counts hydrogens in the fit (otherwise heavy-atom RMSD);
// `symmetry` lets symmetry-equivalent atoms be remapped for the best fit.
// `ok` is false if the two molecules differ in atom count or alignment fails.
// The molecules must have the same atom ordering for the result to be sensible.
double mol_align(Molecule &mol, const Molecule &reference, bool include_h,
                 bool symmetry, bool &ok);

// --- 2D depiction ---------------------------------------------------------
// Generate 2D coordinates in place via the `gen2D` op. Returns false on failure.
bool mol_make_2d(Molecule &mol);
// Render `mol` to an SVG document. `all_carbons` draws a label on every carbon
// (the default labels only terminal carbons); `atom_indices` annotates each
// atom with its index. 2D coordinates are generated automatically if absent.
// `ok` is false on failure (the returned string is then empty).
rust::String mol_to_svg(const Molecule &mol, bool all_carbons, bool atom_indices,
                        bool &ok);

// --- Stereochemistry ------------------------------------------------------
// Force (re)perception of stereochemistry from the molecule's structure
// (SMILES @/@@ and /\, or 2D/3D coordinates).
void mol_perceive_stereo(Molecule &mol);
// Count of tetrahedral / cis-trans stereo units perceived in `mol`.
uint32_t mol_num_tetrahedral_stereo(const Molecule &mol);
uint32_t mol_num_cistrans_stereo(const Molecule &mol);
// Whether atom `idx` (1-based) is a tetrahedral stereocenter.
bool atom_is_tetrahedral_stereo(const Molecule &mol, uint32_t idx);
// Winding at tetrahedral stereocenter `idx`: 1 = clockwise, 2 = anticlockwise,
// 0 = not a (specified) stereocenter.
int atom_tetrahedral_winding(const Molecule &mol, uint32_t idx);
// Whether bond `idx` (0-based) is a cis/trans (double-bond) stereo unit.
bool bond_is_cistrans_stereo(const Molecule &mol, uint32_t idx);

// --- Reaction / SMIRKS-like transforms (OBChemTsfm) -----------------------
// Compile a transformation from a reactant SMARTS to a product SMARTS. Returns
// null if either pattern is invalid.
std::unique_ptr<Transform> transform_new(rust::Str reactant, rust::Str product);
// Apply the transformation to every match in `mol`, editing it in place.
// Returns false if nothing matched or the transform failed.
bool transform_apply(const Transform &t, Molecule &mol);

// --- Conformer search (OBConformerSearch) ---------------------------------
// Run a genetic-algorithm conformer search targeting `count` diverse
// conformers, storing them in `mol`. Requires `mol` to already have a 3D
// structure. Returns the number of conformers now stored.
uint32_t mol_generate_conformers(Molecule &mol, uint32_t count);
// Number of stored conformers.
uint32_t mol_num_conformers(const Molecule &mol);
// Make conformer `index` the active coordinates (no-op if out of range).
void mol_set_conformer(Molecule &mol, uint32_t index);
// Flat [x,y,z,...] of conformer `index`, without changing the active conformer.
rust::Vec<double> mol_conformer_coordinates(const Molecule &mol, uint32_t index);
// Energy of each conformer under `ff_id`, restoring the active conformer.
rust::Vec<double> mol_conformer_energies(const Molecule &mol, rust::Str ff_id);

// --- Element data (OBElements, keyed by atomic number) --------------------
rust::String element_symbol(uint32_t atomic_number);
rust::String element_name(uint32_t atomic_number);
uint32_t element_atomic_number(rust::Str symbol);
double element_mass(uint32_t atomic_number);         // standard atomic weight
double element_exact_mass(uint32_t atomic_number);   // most abundant isotope
double element_electronegativity(uint32_t atomic_number);  // Pauling
double element_covalent_radius(uint32_t atomic_number);    // Angstrom
double element_vdw_radius(uint32_t atomic_number);         // Angstrom
uint32_t element_max_bonds(uint32_t atomic_number);

// --- More atom accessors (idx is 1-based) ---------------------------------
rust::String atom_type(const Molecule &mol, uint32_t idx);  // OB atom type, e.g. "C3"
uint32_t atom_isotope(const Molecule &mol, uint32_t idx);
double atom_atomic_mass(const Molecule &mol, uint32_t idx);
double atom_exact_mass(const Molecule &mol, uint32_t idx);
int atom_spin_multiplicity(const Molecule &mol, uint32_t idx);
uint32_t atom_heavy_degree(const Molecule &mol, uint32_t idx);   // heavy-atom neighbours
uint32_t atom_hetero_degree(const Molecule &mol, uint32_t idx);  // heteroatom neighbours
bool atom_is_chiral(const Molecule &mol, uint32_t idx);
bool atom_is_heteroatom(const Molecule &mol, uint32_t idx);
bool atom_is_metal(const Molecule &mol, uint32_t idx);
bool atom_is_polar_hydrogen(const Molecule &mol, uint32_t idx);
uint32_t atom_member_of_ring_count(const Molecule &mol, uint32_t idx);
uint32_t atom_member_of_ring_size(const Molecule &mol, uint32_t idx);  // smallest ring, 0 if none
bool atom_is_in_ring_size(const Molecule &mol, uint32_t idx, uint32_t size);

// --- More bond accessors (idx is 0-based) ---------------------------------
double bond_length(const Molecule &mol, uint32_t idx);
double bond_equilibrium_length(const Molecule &mol, uint32_t idx);
bool bond_is_rotor(const Molecule &mol, uint32_t idx);
bool bond_is_amide(const Molecule &mol, uint32_t idx);
bool bond_is_ester(const Molecule &mol, uint32_t idx);
bool bond_is_carbonyl(const Molecule &mol, uint32_t idx);
bool bond_is_closure(const Molecule &mol, uint32_t idx);  // ring-closure bond

// --- More whole-molecule methods ------------------------------------------
uint32_t mol_num_heavy_atoms(const Molecule &mol);
uint32_t mol_num_rotors(const Molecule &mol);   // rotatable bonds
uint32_t mol_num_rings(const Molecule &mol);    // SSSR ring count
rust::String mol_spaced_formula(const Molecule &mol);
uint32_t mol_spin_multiplicity(const Molecule &mol);
// Translate the molecule so its centroid is at the origin.
void mol_center(Molecule &mol);
// Bond/valence angle i-j-k in degrees (atom indices 1-based); 0 if invalid.
double mol_angle(const Molecule &mol, uint32_t i, uint32_t j, uint32_t k);
// Torsion angle i-j-k-l in degrees (atom indices 1-based); 0 if invalid.
double mol_torsion(const Molecule &mol, uint32_t i, uint32_t j, uint32_t k, uint32_t l);
// Deep-copy the molecule.
std::unique_ptr<Molecule> mol_clone(const Molecule &mol);
// Remove disconnected fragments with fewer than `threshold` atoms (0 keeps
// only the largest fragment). Returns true if anything was removed.
bool mol_strip_salts(Molecule &mol, uint32_t threshold);
// Split into disconnected fragments (each a fresh molecule).
std::unique_ptr<std::vector<Molecule>> mol_separate(const Molecule &mol);
// String property (OBPairData) access by key.
void mol_set_property(Molecule &mol, rust::Str key, rust::Str value);
rust::String mol_get_property(const Molecule &mol, rust::Str key, bool &ok);

// --- Residues (OBResidue; res_idx is 0-based) -----------------------------
uint32_t mol_num_residues(const Molecule &mol);
rust::String residue_name(const Molecule &mol, uint32_t res_idx);
int residue_number(const Molecule &mol, uint32_t res_idx);
rust::String residue_number_string(const Molecule &mol, uint32_t res_idx);
rust::String residue_chain(const Molecule &mol, uint32_t res_idx);  // 1 char, "" if none
rust::String residue_insertion_code(const Molecule &mol, uint32_t res_idx);
uint32_t residue_num_atoms(const Molecule &mol, uint32_t res_idx);
uint32_t residue_num_heavy_atoms(const Molecule &mol, uint32_t res_idx);
// 0-based indices of the atoms that belong to this residue.
rust::Vec<uint32_t> residue_atom_indices(const Molecule &mol, uint32_t res_idx);
// Per-atom residue info (atom idx is 1-based). residue index is -1 if the atom
// has no residue.
int atom_residue_index(const Molecule &mol, uint32_t idx);
rust::String atom_residue_atom_id(const Molecule &mol, uint32_t idx);  // e.g. "CA"
bool atom_is_hetatm(const Molecule &mol, uint32_t idx);
uint32_t atom_serial_number(const Molecule &mol, uint32_t idx);

// --- Spectra --------------------------------------------------------------
// Spectrophore descriptor (48 values by default). Needs 3D coordinates;
// returns an empty vector otherwise.
rust::Vec<double> mol_spectrophore(const Molecule &mol);
// Vibrational frequencies / IR intensities, if the molecule carries
// OBVibrationData (e.g. read from a computational-chemistry output). Empty
// otherwise.
rust::Vec<double> mol_vibration_frequencies(const Molecule &mol);
rust::Vec<double> mol_vibration_intensities(const Molecule &mol);

// --- Bulk coordinates -----------------------------------------------------
// All atom coordinates flattened as [x0,y0,z0, x1,y1,z1, ...] (3 per atom).
rust::Vec<double> mol_coordinates(const Molecule &mol);

// --- Force-field constraints (atom indices are 0-based here; the shim adds 1
// for OpenBabel's 1-based OBFFConstraints) --------------------------------
std::unique_ptr<Constraints> constraints_new();
void constraints_add_ignore(Constraints &c, uint32_t atom);
void constraints_add_atom(Constraints &c, uint32_t atom);        // fix position
void constraints_add_atom_x(Constraints &c, uint32_t atom);      // fix x only
void constraints_add_atom_y(Constraints &c, uint32_t atom);      // fix y only
void constraints_add_atom_z(Constraints &c, uint32_t atom);      // fix z only
void constraints_add_distance(Constraints &c, uint32_t a, uint32_t b, double length);
void constraints_add_angle(Constraints &c, uint32_t a, uint32_t b, uint32_t d, double angle);
void constraints_add_torsion(Constraints &c, uint32_t a, uint32_t b, uint32_t d,
                             uint32_t e, double torsion);
void constraints_set_factor(Constraints &c, double factor);

// --- Geometry optimization -----------------------------------------------
// Both functions run the WHOLE minimization inside a single call, so the safe
// wrapper can hold OpenBabel's global lock for the entire run. That atomicity
// matters: OpenBabel's force-field constraint state (`_constraints`) is a static
// (per-class) member shared by every instance, so interleaving two optimizations
// would corrupt it. Each call also clears that static constraint set before
// returning, restoring the "no constraints" state the other force-field calls
// (mol_energy / mol_optimize) rely on.
//
// `algorithm` is 0 = steepest descent, 1 = conjugate gradients, 2 = L-BFGS;
// `steps` is the step budget, `econv` the energy-convergence criterion, and
// `constraints` the restraint set (empty for none).

// Run to completion; write final coordinates back to `mol` and return the final
// energy. Sets `ok` to false (energy NaN) if the force field is unknown or setup
// fails.
double optimizer_run_to_end(Molecule &mol, rust::Str ff_id, uint32_t algorithm,
                            uint32_t steps, double econv,
                            const Constraints &constraints, bool &ok);

// Run to completion recording a frame every `frame_interval` steps. The result
// is flattened as repeated frames, each `1 + 3 * num_atoms` doubles laid out as
// [energy, x0, y0, z0, x1, y1, z1, ...]. `mol` ends at the final geometry. Empty
// on unknown force field or setup failure.
rust::Vec<double> optimizer_run_trajectory(Molecule &mol, rust::Str ff_id,
                                           uint32_t algorithm, uint32_t steps,
                                           double econv, const Constraints &constraints,
                                           uint32_t frame_interval);

// --- Molecule construction & editing --------------------------------------
// Add a new atom of `atomic_num`; returns its 0-based index.
uint32_t mol_add_atom(Molecule &mol, uint32_t atomic_num);
// Add a bond between 0-based atoms `begin` and `end` with `order` (1/2/3);
// false if the indices are out of range.
bool mol_add_bond(Molecule &mol, uint32_t begin, uint32_t end, uint32_t order);
// Delete atom / bond at a 0-based index; false if out of range.
bool mol_delete_atom(Molecule &mol, uint32_t idx);
bool mol_delete_bond(Molecule &mol, uint32_t idx);
// Suspend / resume perception around a batch of edits (nestable).
void mol_begin_modify(Molecule &mol);
void mol_end_modify(Molecule &mol);
// Remove every atom and bond, leaving an empty molecule.
void mol_clear(Molecule &mol);
// Translate every atom by (x, y, z).
void mol_translate(Molecule &mol, double x, double y, double z);
// Overwrite all coordinates from a flat [x0,y0,z0, ...] slice; false unless its
// length is exactly 3 * num_atoms.
bool mol_set_coordinates(Molecule &mol, rust::Slice<const double> coords);
// Set the coordinate dimension (0, 2, or 3). Needed to mark a hand-built
// structure as 3D before mol_connect_the_dots (which only runs when dim == 3).
void mol_set_dimension(Molecule &mol, uint32_t dim);
// Infer connectivity from 3D coordinates (covalent radii), then bond orders.
void mol_connect_the_dots(Molecule &mol);
void mol_perceive_bond_orders(Molecule &mol);
// Add only polar hydrogens (on N/O/P/S). false on failure.
bool mol_add_polar_hydrogens(Molecule &mol);
// Convert dative bonds (e.g. -[N+]([O-])= nitro to neutral dative form). false
// if nothing changed.
bool mol_convert_dative_bonds(Molecule &mol);
// (Re)assign radical spin multiplicities from valence. false on failure.
bool mol_assign_spin_multiplicity(Molecule &mol);
// Add hydrogens with pH-based (de)protonation correction at `ph`. false on
// failure.
bool mol_add_hydrogens_ph(Molecule &mol, double ph);

// --- Atom setters (idx is 1-based, as elsewhere) --------------------------
void atom_set_atomic_num(Molecule &mol, uint32_t idx, uint32_t atomic_num);
void atom_set_formal_charge(Molecule &mol, uint32_t idx, int charge);
void atom_set_position(Molecule &mol, uint32_t idx, double x, double y, double z);
void atom_set_isotope(Molecule &mol, uint32_t idx, uint32_t isotope);
void atom_set_spin_multiplicity(Molecule &mol, uint32_t idx, int spin);
void atom_set_partial_charge(Molecule &mol, uint32_t idx, double charge);
void atom_set_type(Molecule &mol, uint32_t idx, rust::Str type_name);
void atom_set_implicit_h(Molecule &mol, uint32_t idx, uint32_t count);

// --- Bond setters (idx is 0-based) ----------------------------------------
void bond_set_order(Molecule &mol, uint32_t idx, uint32_t order);
// Move the end atom so the bond has `length` (keeps the begin atom fixed).
bool bond_set_length(Molecule &mol, uint32_t idx, double length);

// --- Multi-molecule input -------------------------------------------------
// Read EVERY record from `data` in `format` (e.g. multi-record SDF, one SMILES
// per line). Empty on unknown format; stops at the first record that fails.
std::unique_ptr<std::vector<Molecule>> mol_read_many(rust::Str format, rust::Str data);

// --- Ring access (SSSR; ring_idx is 0-based, 0..mol_num_rings) -------------
uint32_t ring_size(const Molecule &mol, uint32_t ring_idx);
// 0-based atom indices forming the ring.
rust::Vec<uint32_t> ring_atom_indices(const Molecule &mol, uint32_t ring_idx);
bool ring_is_aromatic(const Molecule &mol, uint32_t ring_idx);

// --- Graph navigation (atom idx 1-based; returned atom/bond idx 0-based) ---
// 0-based indices of the atoms bonded to atom `idx`.
rust::Vec<uint32_t> atom_neighbor_indices(const Molecule &mol, uint32_t idx);
// 0-based indices of the bonds incident to atom `idx`.
rust::Vec<uint32_t> atom_bond_indices(const Molecule &mol, uint32_t idx);
// Number of bonds from atom `idx` with the given order.
uint32_t atom_count_bonds_of_order(const Molecule &mol, uint32_t idx, uint32_t order);
// Number of explicit (present-in-graph) hydrogens on atom `idx`.
uint32_t atom_explicit_h_count(const Molecule &mol, uint32_t idx);
// 0-based bond index joining 0-based atoms `a` and `b`; -1 if not bonded.
int mol_bond_between(const Molecule &mol, uint32_t a, uint32_t b);
// 0-based index of the atom at the other end of bond `bond_idx` from 0-based
// atom `atom_idx`; -1 if that atom is not on the bond.
int bond_other_atom(const Molecule &mol, uint32_t bond_idx, uint32_t atom_idx);

// --- Crystallography (unit cell) ------------------------------------------
bool mol_has_unit_cell(const Molecule &mol);
// [a, b, c, alpha, beta, gamma] (lengths in Å, angles in degrees); empty if no
// unit cell.
rust::Vec<double> mol_cell_parameters(const Molecule &mol);
double mol_cell_volume(const Molecule &mol);
rust::String mol_cell_spacegroup(const Molecule &mol);
// LatticeType: 0 undefined,1 triclinic,2 monoclinic,3 orthorhombic,4 tetragonal,
// 5 rhombohedral,6 hexagonal,7 cubic.
uint32_t mol_cell_lattice_type(const Molecule &mol);
// Convert Cartesian <-> fractional coordinates; each returns [x,y,z] (empty if
// no unit cell).
rust::Vec<double> mol_cell_to_fractional(const Molecule &mol, double x, double y, double z);
rust::Vec<double> mol_cell_to_cartesian(const Molecule &mol, double x, double y, double z);

// --- Symmetry & canonical ordering (one value per atom, atom order) --------
// Topological symmetry classes; atoms sharing a value are graph-equivalent.
rust::Vec<uint32_t> mol_symmetry_classes(const Molecule &mol);
// Canonical rank (1-based) per atom — a repeatable canonical labelling.
rust::Vec<uint32_t> mol_canonical_ranks(const Molecule &mol);

// --- Reactions (OBReaction; formats "rxn" MDL RXN, "rsmi" reaction SMILES) --
std::unique_ptr<Reaction> reaction_new();
// Parse a reaction document; null on unknown format or parse failure.
std::unique_ptr<Reaction> reaction_read(rust::Str format, rust::Str data);
// Serialize; sets ok=false (and returns empty) on unknown format / failure.
rust::String reaction_write(const Reaction &r, rust::Str format, bool &ok);
uint32_t reaction_num_reactants(const Reaction &r);
uint32_t reaction_num_products(const Reaction &r);
uint32_t reaction_num_agents(const Reaction &r);
// Deep copy of the i-th component as a standalone Molecule; null if out of range.
std::unique_ptr<Molecule> reaction_reactant(const Reaction &r, uint32_t i);
std::unique_ptr<Molecule> reaction_product(const Reaction &r, uint32_t i);
std::unique_ptr<Molecule> reaction_agent(const Reaction &r, uint32_t i);
// Append a deep copy of mol to the respective component list.
void reaction_add_reactant(Reaction &r, const Molecule &mol);
void reaction_add_product(Reaction &r, const Molecule &mol);
void reaction_add_agent(Reaction &r, const Molecule &mol);
rust::String reaction_title(const Reaction &r);
void reaction_set_title(Reaction &r, rust::Str title);
rust::String reaction_comment(const Reaction &r);
void reaction_set_comment(Reaction &r, rust::Str comment);
bool reaction_is_reversible(const Reaction &r);
void reaction_set_reversible(Reaction &r, bool value);

// --- File I/O (format auto-detected from extension when `format` empty) ----
// Read the first molecule from a file; null on open/parse failure.
std::unique_ptr<Molecule> mol_read_file(rust::Str path, rust::Str format);
// Read every molecule from a (possibly multi-record) file; empty on failure.
std::unique_ptr<std::vector<Molecule>> mol_read_file_many(rust::Str path, rust::Str format);
// Write mol to a file; sets ok=false on unknown format / write failure.
void mol_write_file(const Molecule &mol, rust::Str path, rust::Str format, bool &ok);

// --- Subgraph isomorphism & automorphisms ---------------------------------
// Unique mappings of `query` as a substructure of `target`; sets `width` to the
// query atom count, flat result = `width` target atom indices (0-based) per
// mapping, ordered by query atom index. Empty if no match.
rust::Vec<uint32_t> mol_substructure_mappings(const Molecule &query, const Molecule &target,
                                              uint32_t &width);
// All graph automorphisms of `mol`; sets `width` to the atom count, flat result
// = one atom-index permutation (0-based) per automorphism.
rust::Vec<uint32_t> mol_automorphisms(const Molecule &mol, uint32_t &width);

// --- Geometry & topology (niche; atom indices 0-based) --------------------
// Set the a-b-c-d torsion to `radians`, rotating the b-c bond's far side.
void mol_set_torsion(Molecule &mol, uint32_t a, uint32_t b, uint32_t c, uint32_t d,
                     double radians);
// Atoms reachable from `to` without passing back through `from` (excludes both).
rust::Vec<uint32_t> mol_find_children(const Molecule &mol, uint32_t from, uint32_t to);
// Atom indices of the largest connected fragment.
rust::Vec<uint32_t> mol_largest_fragment(const Molecule &mol);
void mol_set_total_charge(Molecule &mol, int32_t charge);
void mol_set_total_spin(Molecule &mol, uint32_t spin);

// --- Per-atom / per-bond string data (OBPairData) -------------------------
// atom idx is 1-based, bond idx is 0-based. get sets ok=false when absent.
void atom_set_data(Molecule &mol, uint32_t idx, rust::Str key, rust::Str value);
rust::String atom_get_data(const Molecule &mol, uint32_t idx, rust::Str key, bool &ok);
void bond_set_data(Molecule &mol, uint32_t idx, rust::Str key, rust::Str value);
rust::String bond_get_data(const Molecule &mol, uint32_t idx, rust::Str key, bool &ok);

// --- Inter-atom distance & 2D wedge/hash bond stereo ----------------------
// Distance (Å) between atoms `i`/`j` (1-based); 0.0 for invalid indices.
double mol_distance(const Molecule &mol, uint32_t i, uint32_t j);
// Bond idx is 0-based. Wedge/hash mark 2D depiction stereo direction.
bool bond_is_wedge(const Molecule &mol, uint32_t idx);
bool bond_is_hash(const Molecule &mol, uint32_t idx);
void bond_set_wedge(Molecule &mol, uint32_t idx, bool value);
void bond_set_hash(Molecule &mol, uint32_t idx, bool value);

}  // namespace ob_shim
