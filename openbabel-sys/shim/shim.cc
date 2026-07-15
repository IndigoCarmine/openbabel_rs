#include "shim.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>

#include <cstdlib>
#include <sstream>
#include <string>

// `struct Molecule` is defined completely in shim.h (cxx needs it there).

namespace ob_shim {

namespace {

std::string to_std(rust::Str s) { return std::string(s.data(), s.size()); }

// Fetch atom `idx` (1-based). Returns nullptr if out of range.
const OpenBabel::OBAtom *atom_at(const Molecule &mol, uint32_t idx) {
  return const_cast<OpenBabel::OBMol &>(mol.mol).GetAtom(static_cast<int>(idx));
}

// Fetch bond `idx` (0-based). Returns nullptr if out of range.
const OpenBabel::OBBond *bond_at(const Molecule &mol, uint32_t idx) {
  return const_cast<OpenBabel::OBMol &>(mol.mol).GetBond(static_cast<int>(idx));
}

}  // namespace

rust::String release_version() {
  return rust::String(OpenBabel::OBReleaseVersion());
}

void set_env(rust::Str key, rust::Str value) {
  std::string k = to_std(key);
  std::string v = to_std(value);
#ifdef _WIN32
  _putenv_s(k.c_str(), v.c_str());
#else
  setenv(k.c_str(), v.c_str(), 1);
#endif
}

std::unique_ptr<Molecule> mol_new() {
  return std::unique_ptr<Molecule>(new Molecule());
}

std::unique_ptr<Molecule> mol_read(rust::Str format, rust::Str data) {
  try {
    OpenBabel::OBConversion conv;
    if (!conv.SetInFormat(to_std(format).c_str())) return nullptr;
    auto m = std::unique_ptr<Molecule>(new Molecule());
    std::istringstream iss(to_std(data));
    if (!conv.Read(&m->mol, &iss)) return nullptr;
    return m;
  } catch (...) {
    return nullptr;
  }
}

rust::String mol_write(const Molecule &mol, rust::Str format, bool &ok) {
  try {
    OpenBabel::OBConversion conv;
    if (!conv.SetOutFormat(to_std(format).c_str())) {
      ok = false;
      return rust::String();
    }
    ok = true;
    std::string out = conv.WriteString(&const_cast<Molecule &>(mol).mol);
    return rust::String(out);
  } catch (...) {
    ok = false;
    return rust::String();
  }
}

rust::String mol_formula(const Molecule &mol) {
  try {
    return rust::String(const_cast<Molecule &>(mol).mol.GetFormula());
  } catch (...) {
    return rust::String();
  }
}

double mol_mol_wt(const Molecule &mol) {
  try {
    return const_cast<Molecule &>(mol).mol.GetMolWt();
  } catch (...) {
    return 0.0;
  }
}

double mol_exact_mass(const Molecule &mol) {
  try {
    return const_cast<Molecule &>(mol).mol.GetExactMass();
  } catch (...) {
    return 0.0;
  }
}

int mol_total_charge(const Molecule &mol) {
  try {
    return const_cast<Molecule &>(mol).mol.GetTotalCharge();
  } catch (...) {
    return 0;
  }
}

uint32_t mol_num_atoms(const Molecule &mol) { return mol.mol.NumAtoms(); }
uint32_t mol_num_bonds(const Molecule &mol) { return mol.mol.NumBonds(); }

rust::String mol_title(const Molecule &mol) {
  return rust::String(const_cast<Molecule &>(mol).mol.GetTitle());
}

void mol_set_title(Molecule &mol, rust::Str title) {
  mol.mol.SetTitle(to_std(title).c_str());
}

void mol_add_hydrogens(Molecule &mol) {
  try {
    mol.mol.AddHydrogens();
  } catch (...) {
  }
}

void mol_delete_hydrogens(Molecule &mol) {
  try {
    mol.mol.DeleteHydrogens();
  } catch (...) {
  }
}

uint32_t atom_atomic_num(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetAtomicNum() : 0;
}
double atom_x(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetX() : 0.0;
}
double atom_y(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetY() : 0.0;
}
double atom_z(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetZ() : 0.0;
}
int atom_formal_charge(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetFormalCharge() : 0;
}
bool atom_is_aromatic(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsAromatic() : false;
}
bool atom_is_in_ring(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsInRing() : false;
}

uint32_t bond_begin_idx(const Molecule &mol, uint32_t idx) {
  const auto *b = bond_at(mol, idx);
  return b ? b->GetBeginAtomIdx() : 0;
}
uint32_t bond_end_idx(const Molecule &mol, uint32_t idx) {
  const auto *b = bond_at(mol, idx);
  return b ? b->GetEndAtomIdx() : 0;
}
uint32_t bond_order(const Molecule &mol, uint32_t idx) {
  const auto *b = bond_at(mol, idx);
  return b ? b->GetBondOrder() : 0;
}
bool bond_is_aromatic(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsAromatic() : false;
}
bool bond_is_in_ring(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsInRing() : false;
}

}  // namespace ob_shim
