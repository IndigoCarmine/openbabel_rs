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
bool atom_is_aromatic(const Molecule &mol, uint32_t idx);
bool atom_is_in_ring(const Molecule &mol, uint32_t idx);

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

}  // namespace ob_shim
