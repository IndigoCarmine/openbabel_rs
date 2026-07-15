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

#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/phmodel.h>

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

}  // namespace ob_shim
