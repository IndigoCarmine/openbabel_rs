#include "shim.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/chargemodel.h>
#include <openbabel/descriptor.h>
#include <openbabel/fingerprint.h>
#include <openbabel/conformersearch.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/align.h>
#include <openbabel/obconversion.h>
#include <openbabel/op.h>
#include <openbabel/plugin.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/tetrahedral.h>

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

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
double atom_partial_charge(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->GetPartialCharge() : 0.0;
}
bool atom_is_aromatic(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsAromatic() : false;
}
bool atom_is_in_ring(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsInRing() : false;
}
uint32_t atom_degree(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetExplicitDegree() : 0;
}
uint32_t atom_total_valence(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetTotalValence() : 0;
}
uint32_t atom_implicit_h_count(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetImplicitHCount() : 0;
}
uint32_t atom_hybridization(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetHyb() : 0;
}
bool atom_is_hbond_donor(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsHbondDonor() : false;
}
bool atom_is_hbond_acceptor(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsHbondAcceptor() : false;
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

// --- SMARTS --------------------------------------------------------------

std::unique_ptr<Smarts> smarts_new(rust::Str pattern) {
  try {
    auto s = std::unique_ptr<Smarts>(new Smarts());
    if (!s->pat.Init(to_std(pattern))) return nullptr;
    return s;
  } catch (...) {
    return nullptr;
  }
}

uint32_t smarts_atom_count(const Smarts &smarts) {
  return smarts.pat.NumAtoms();
}

bool smarts_matches(const Smarts &smarts, const Molecule &mol) {
  try {
    std::vector<std::vector<int>> mlist;
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    smarts.pat.Match(m, mlist, OpenBabel::OBSmartsPattern::Single);
    return !mlist.empty();
  } catch (...) {
    return false;
  }
}

rust::Vec<uint32_t> smarts_match_atoms(const Smarts &smarts,
                                       const Molecule &mol) {
  rust::Vec<uint32_t> out;
  try {
    std::vector<std::vector<int>> mlist;
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    smarts.pat.Match(m, mlist, OpenBabel::OBSmartsPattern::AllUnique);
    for (const auto &match : mlist)
      for (int idx : match) out.push_back(static_cast<uint32_t>(idx));
  } catch (...) {
  }
  return out;
}

// --- Fingerprints --------------------------------------------------------

rust::Vec<uint32_t> fingerprint(const Molecule &mol, rust::Str id) {
  rust::Vec<uint32_t> out;
  try {
    // Use OBPlugin::GetPlugin (exported from the OpenBabel DLL) rather than
    // OBFingerprint::FindFingerprint (an inline header function): the inline
    // path consults a plugin map local to *this* translation unit, which is
    // empty because the .obf plugins register into the DLL's own map.
    OpenBabel::OBFingerprint *fp = dynamic_cast<OpenBabel::OBFingerprint *>(
        OpenBabel::OBPlugin::GetPlugin("fingerprints", to_std(id).c_str()));
    if (!fp) return out;
    std::vector<unsigned int> v;
    if (!fp->GetFingerprint(&const_cast<Molecule &>(mol).mol, v)) return out;
    for (unsigned int x : v) out.push_back(static_cast<uint32_t>(x));
  } catch (...) {
  }
  return out;
}

double tanimoto(rust::Slice<const uint32_t> a, rust::Slice<const uint32_t> b) {
  try {
    std::vector<unsigned int> va(a.begin(), a.end());
    std::vector<unsigned int> vb(b.begin(), b.end());
    return OpenBabel::OBFingerprint::Tanimoto(va, vb);
  } catch (...) {
    return 0.0;
  }
}

// --- Descriptors ---------------------------------------------------------

double descriptor(const Molecule &mol, rust::Str id, bool &ok) {
  try {
    // GetPlugin (DLL-exported) instead of OBDescriptor::FindType (inline);
    // see the note in fingerprint() above.
    OpenBabel::OBDescriptor *d = dynamic_cast<OpenBabel::OBDescriptor *>(
        OpenBabel::OBPlugin::GetPlugin("descriptors", to_std(id).c_str()));
    if (!d) {
      ok = false;
      return 0.0;
    }
    ok = true;
    return d->Predict(&const_cast<Molecule &>(mol).mol);
  } catch (...) {
    ok = false;
    return 0.0;
  }
}

// --- Force fields --------------------------------------------------------

namespace {
// Look up a force field plugin by id via the DLL-exported plugin registry.
OpenBabel::OBForceField *find_forcefield(const std::string &id) {
  return dynamic_cast<OpenBabel::OBForceField *>(
      OpenBabel::OBPlugin::GetPlugin("forcefields", id.c_str()));
}
}  // namespace

double mol_energy(const Molecule &mol, rust::Str ff_id, bool &ok) {
  try {
    OpenBabel::OBForceField *ff = find_forcefield(to_std(ff_id));
    if (!ff) {
      ok = false;
      return std::nan("");
    }
    if (!ff->Setup(const_cast<Molecule &>(mol).mol)) {
      ok = false;
      return std::nan("");
    }
    ok = true;
    return ff->Energy(false);
  } catch (...) {
    ok = false;
    return std::nan("");
  }
}

rust::String forcefield_unit(rust::Str ff_id) {
  try {
    OpenBabel::OBForceField *ff = find_forcefield(to_std(ff_id));
    return ff ? rust::String(ff->GetUnit()) : rust::String();
  } catch (...) {
    return rust::String();
  }
}

double mol_optimize(Molecule &mol, rust::Str ff_id, uint32_t steps, bool &ok) {
  try {
    OpenBabel::OBForceField *ff = find_forcefield(to_std(ff_id));
    if (!ff || !ff->Setup(mol.mol)) {
      ok = false;
      return std::nan("");
    }
    ff->ConjugateGradients(static_cast<int>(steps));
    ff->GetCoordinates(mol.mol);
    ok = true;
    return ff->Energy(false);
  } catch (...) {
    ok = false;
    return std::nan("");
  }
}

// --- 3D structure generation ---------------------------------------------

bool mol_make_3d(Molecule &mol, rust::Str speed) {
  try {
    OpenBabel::OBOp *op = dynamic_cast<OpenBabel::OBOp *>(
        OpenBabel::OBPlugin::GetPlugin("ops", "gen3d"));
    if (!op) return false;
    return op->Do(&mol.mol, to_std(speed).c_str(), nullptr, nullptr);
  } catch (...) {
    return false;
  }
}

uint32_t mol_dimension(const Molecule &mol) { return mol.mol.GetDimension(); }

// --- Partial charges -----------------------------------------------------

bool mol_compute_charges(Molecule &mol, rust::Str model) {
  try {
    OpenBabel::OBChargeModel *cm = dynamic_cast<OpenBabel::OBChargeModel *>(
        OpenBabel::OBPlugin::GetPlugin("charges", to_std(model).c_str()));
    if (!cm) return false;
    return cm->ComputeCharges(mol.mol);
  } catch (...) {
    return false;
  }
}

// --- Structure alignment -------------------------------------------------

double mol_align(Molecule &mol, const Molecule &reference, bool include_h,
                 bool symmetry, bool &ok) {
  try {
    // OBAlign superposes atoms in order, so the molecules must have matching
    // atom counts; reject a mismatch rather than produce a meaningless fit.
    if (mol.mol.NumAtoms() != reference.mol.NumAtoms()) {
      ok = false;
      return 0.0;
    }
    OpenBabel::OBAlign aln(include_h, symmetry);
    aln.SetRefMol(reference.mol);
    aln.SetTargetMol(mol.mol);
    if (!aln.Align()) {
      ok = false;
      return 0.0;
    }
    aln.UpdateCoords(&mol.mol);
    ok = true;
    return aln.GetRMSD();
  } catch (...) {
    ok = false;
    return 0.0;
  }
}

// --- 2D depiction --------------------------------------------------------

bool mol_make_2d(Molecule &mol) {
  try {
    OpenBabel::OBOp *op = dynamic_cast<OpenBabel::OBOp *>(
        OpenBabel::OBPlugin::GetPlugin("ops", "gen2D"));
    if (!op) return false;
    return op->Do(&mol.mol, "", nullptr, nullptr);
  } catch (...) {
    return false;
  }
}

rust::String mol_to_svg(const Molecule &mol, bool all_carbons, bool atom_indices,
                        bool &ok) {
  try {
    OpenBabel::OBConversion conv;
    if (!conv.SetOutFormat("svg")) {
      ok = false;
      return rust::String();
    }
    // The SVG format generates 2D coordinates itself (via gen2D) when the
    // molecule has none, so a freshly parsed molecule renders directly.
    if (all_carbons) conv.AddOption("a", OpenBabel::OBConversion::OUTOPTIONS);
    if (atom_indices) conv.AddOption("i", OpenBabel::OBConversion::OUTOPTIONS);
    ok = true;
    return rust::String(conv.WriteString(&const_cast<Molecule &>(mol).mol));
  } catch (...) {
    ok = false;
    return rust::String();
  }
}

// --- Stereochemistry -----------------------------------------------------

void mol_perceive_stereo(Molecule &mol) {
  try {
    OpenBabel::PerceiveStereo(&mol.mol, true);
  } catch (...) {
  }
}

uint32_t mol_num_tetrahedral_stereo(const Molecule &mol) {
  try {
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    OpenBabel::OBStereoFacade facade(&m);
    uint32_t n = 0;
    for (uint32_t i = 1; i <= m.NumAtoms(); ++i) {
      OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(i));
      if (a && facade.HasTetrahedralStereo(a->GetId())) ++n;
    }
    return n;
  } catch (...) {
    return 0;
  }
}

uint32_t mol_num_cistrans_stereo(const Molecule &mol) {
  try {
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    OpenBabel::OBStereoFacade facade(&m);
    uint32_t n = 0;
    for (uint32_t i = 0; i < m.NumBonds(); ++i) {
      OpenBabel::OBBond *b = m.GetBond(static_cast<int>(i));
      if (b && facade.HasCisTransStereo(b->GetId())) ++n;
    }
    return n;
  } catch (...) {
    return 0;
  }
}

bool atom_is_tetrahedral_stereo(const Molecule &mol, uint32_t idx) {
  try {
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(idx));
    if (!a) return false;
    OpenBabel::OBStereoFacade facade(&m);
    return facade.HasTetrahedralStereo(a->GetId());
  } catch (...) {
    return false;
  }
}

int atom_tetrahedral_winding(const Molecule &mol, uint32_t idx) {
  try {
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(idx));
    if (!a) return 0;
    OpenBabel::OBStereoFacade facade(&m);
    if (!facade.HasTetrahedralStereo(a->GetId())) return 0;
    OpenBabel::OBTetrahedralStereo *ts = facade.GetTetrahedralStereo(a->GetId());
    if (!ts) return 0;
    OpenBabel::OBTetrahedralStereo::Config cfg = ts->GetConfig();
    if (!cfg.specified) return 0;
    return static_cast<int>(cfg.winding);  // 1 = clockwise, 2 = anticlockwise
  } catch (...) {
    return 0;
  }
}

bool bond_is_cistrans_stereo(const Molecule &mol, uint32_t idx) {
  try {
    OpenBabel::OBMol &m = const_cast<OpenBabel::OBMol &>(mol.mol);
    OpenBabel::OBBond *b = m.GetBond(static_cast<int>(idx));
    if (!b) return false;
    OpenBabel::OBStereoFacade facade(&m);
    return facade.HasCisTransStereo(b->GetId());
  } catch (...) {
    return false;
  }
}

// --- Reaction / SMIRKS-like transforms -----------------------------------

std::unique_ptr<Transform> transform_new(rust::Str reactant, rust::Str product) {
  try {
    auto t = std::unique_ptr<Transform>(new Transform());
    std::string start = to_std(reactant);
    std::string end = to_std(product);
    if (!t->tsfm.Init(start, end)) return nullptr;
    return t;
  } catch (...) {
    return nullptr;
  }
}

bool transform_apply(const Transform &t, Molecule &mol) {
  try {
    return const_cast<Transform &>(t).tsfm.Apply(mol.mol);
  } catch (...) {
    return false;
  }
}

// --- Conformer search ----------------------------------------------------

uint32_t mol_generate_conformers(Molecule &mol, uint32_t count) {
  try {
    OpenBabel::OBConformerSearch cs;
    if (!cs.Setup(mol.mol, static_cast<int>(count))) {
      // No rotatable bonds (or setup failed): the existing structure is the
      // only conformer.
      return mol.mol.NumConformers();
    }
    cs.Search();
    cs.GetConformers(mol.mol);
    return mol.mol.NumConformers();
  } catch (...) {
    return 0;
  }
}

uint32_t mol_num_conformers(const Molecule &mol) {
  return static_cast<uint32_t>(const_cast<Molecule &>(mol).mol.NumConformers());
}

void mol_set_conformer(Molecule &mol, uint32_t index) {
  try {
    if (index < static_cast<uint32_t>(mol.mol.NumConformers()))
      mol.mol.SetConformer(index);
  } catch (...) {
  }
}

}  // namespace ob_shim
