#include "shim.h"

#include <openbabel/atom.h>
#include <openbabel/bitvec.h>
#include <openbabel/bond.h>
#include <openbabel/canon.h>
#include <openbabel/chargemodel.h>
#include <openbabel/descriptor.h>
#include <openbabel/fingerprint.h>
#include <openbabel/conformersearch.h>
#include <openbabel/elements.h>
#include <openbabel/forcefield.h>
#include <openbabel/generic.h>
#include <openbabel/graphsym.h>
#include <openbabel/isomorphism.h>
#include <openbabel/math/align.h>
#include <openbabel/math/vector3.h>
#include <openbabel/obconversion.h>
#include <openbabel/op.h>
#include <openbabel/plugin.h>
#include <openbabel/query.h>
#include <openbabel/residue.h>
#include <openbabel/ring.h>
#include <openbabel/spectrophore.h>
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

// Build a rust::String from OpenBabel text WITHOUT ever throwing.
//
// rust::String's ordinary constructor validates UTF-8 and throws
// std::invalid_argument on any invalid byte. Text that originates in an input
// file — molecule/reaction titles, atom types, residue names, string
// properties, space-group names, … — may hold arbitrary bytes, so that
// exception is reachable with perfectly ordinary input. Because the
// cxx-generated FFI wrappers are `noexcept`, an escaping exception would call
// std::terminate() and abort the whole process. Lossy conversion substitutes
// U+FFFD for invalid sequences and is itself noexcept, so it both prevents the
// abort and preserves the (mostly valid) text instead of discarding it.
rust::String to_rust(const std::string &s) { return rust::String::lossy(s); }
rust::String to_rust(const char *s) {
  return s ? rust::String::lossy(s) : rust::String();
}

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
  return to_rust(OpenBabel::OBReleaseVersion());
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
    return to_rust(out);
  } catch (...) {
    ok = false;
    return rust::String();
  }
}

rust::String mol_formula(const Molecule &mol) {
  try {
    return to_rust(const_cast<Molecule &>(mol).mol.GetFormula());
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
  return to_rust(const_cast<Molecule &>(mol).mol.GetTitle());
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
    return ff ? to_rust(ff->GetUnit()) : rust::String();
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
    return to_rust(conv.WriteString(&const_cast<Molecule &>(mol).mol));
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

// Flat [x,y,z,...] coordinates of conformer `index` (3*NumAtoms doubles),
// read WITHOUT changing the active conformer. Empty if out of range.
rust::Vec<double> mol_conformer_coordinates(const Molecule &mol, uint32_t index) {
  rust::Vec<double> out;
  OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
  if (index >= static_cast<uint32_t>(m.NumConformers())) return out;
  const double *c = m.GetConformer(static_cast<int>(index));
  if (!c) return out;
  uint32_t n = static_cast<uint32_t>(m.NumAtoms());
  for (uint32_t i = 0; i < 3 * n; ++i) out.push_back(c[i]);
  return out;
}

// Energy of every stored conformer under force field `ff_id`, in conformer
// order. The active conformer is restored before returning, so this does not
// change the molecule. Empty if `ff_id` is unknown or there are no conformers.
rust::Vec<double> mol_conformer_energies(const Molecule &mol, rust::Str ff_id) {
  rust::Vec<double> out;
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    int n = m.NumConformers();
    if (n == 0) return out;
    OpenBabel::OBForceField *ff = find_forcefield(to_std(ff_id));
    if (!ff || !ff->Setup(m)) return out;
    // Remember the active conformer so we can restore it afterwards.
    const double *saved = m.GetCoordinates();
    int saved_idx = 0;
    for (int i = 0; i < n; ++i)
      if (m.GetConformer(i) == saved) {
        saved_idx = i;
        break;
      }
    for (int i = 0; i < n; ++i) {
      m.SetConformer(i);
      // Push the now-active conformer's coordinates into the field without a
      // full re-Setup (which would be skipped for this same-topology molecule
      // and keep stale coordinates).
      ff->SetCoordinates(m);
      out.push_back(ff->Energy(false));
    }
    m.SetConformer(saved_idx);
    return out;
  } catch (...) {
    return out;
  }
}

// --- Element data --------------------------------------------------------

rust::String element_symbol(uint32_t z) {
  return to_rust(OpenBabel::OBElements::GetSymbol(z));
}
rust::String element_name(uint32_t z) {
  return to_rust(OpenBabel::OBElements::GetName(z));
}
uint32_t element_atomic_number(rust::Str symbol) {
  return OpenBabel::OBElements::GetAtomicNum(to_std(symbol).c_str());
}
double element_mass(uint32_t z) { return OpenBabel::OBElements::GetMass(z); }
double element_exact_mass(uint32_t z) {
  return OpenBabel::OBElements::GetExactMass(z, 0);
}
double element_electronegativity(uint32_t z) {
  return OpenBabel::OBElements::GetElectroNeg(z);
}
double element_covalent_radius(uint32_t z) {
  return OpenBabel::OBElements::GetCovalentRad(z);
}
double element_vdw_radius(uint32_t z) {
  return OpenBabel::OBElements::GetVdwRad(z);
}
uint32_t element_max_bonds(uint32_t z) {
  return OpenBabel::OBElements::GetMaxBonds(z);
}

// --- More atom accessors -------------------------------------------------

rust::String atom_type(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? to_rust(a->GetType()) : rust::String();
}
uint32_t atom_isotope(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetIsotope() : 0;
}
double atom_atomic_mass(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetAtomicMass() : 0.0;
}
double atom_exact_mass(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetExactMass() : 0.0;
}
int atom_spin_multiplicity(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetSpinMultiplicity() : 0;
}
uint32_t atom_heavy_degree(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetHvyDegree() : 0;
}
uint32_t atom_hetero_degree(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->GetHeteroDegree() : 0;
}
bool atom_is_chiral(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsChiral() : false;
}
bool atom_is_heteroatom(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsHeteroatom() : false;
}
bool atom_is_metal(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsMetal() : false;
}
bool atom_is_polar_hydrogen(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsPolarHydrogen() : false;
}
uint32_t atom_member_of_ring_count(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->MemberOfRingCount() : 0;
}
uint32_t atom_member_of_ring_size(const Molecule &mol, uint32_t idx) {
  const auto *a = atom_at(mol, idx);
  return a ? a->MemberOfRingSize() : 0;
}
bool atom_is_in_ring_size(const Molecule &mol, uint32_t idx, uint32_t size) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->IsInRingSize(static_cast<int>(size)) : false;
}

// --- More bond accessors -------------------------------------------------

double bond_length(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->GetLength() : 0.0;
}
double bond_equilibrium_length(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->GetEquibLength() : 0.0;
}
bool bond_is_rotor(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsRotor() : false;
}
bool bond_is_amide(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsAmide() : false;
}
bool bond_is_ester(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsEster() : false;
}
bool bond_is_carbonyl(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsCarbonyl() : false;
}
bool bond_is_closure(const Molecule &mol, uint32_t idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsClosure() : false;
}

// --- More whole-molecule methods -----------------------------------------

uint32_t mol_num_heavy_atoms(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.NumHvyAtoms();
}
uint32_t mol_num_rotors(const Molecule &mol) {
  try {
    return const_cast<Molecule &>(mol).mol.NumRotors();
  } catch (...) {
    return 0;
  }
}
uint32_t mol_num_rings(const Molecule &mol) {
  try {
    return static_cast<uint32_t>(const_cast<Molecule &>(mol).mol.GetSSSR().size());
  } catch (...) {
    return 0;
  }
}
rust::String mol_spaced_formula(const Molecule &mol) {
  try {
    return to_rust(const_cast<Molecule &>(mol).mol.GetSpacedFormula());
  } catch (...) {
    return rust::String();
  }
}
uint32_t mol_spin_multiplicity(const Molecule &mol) {
  try {
    return const_cast<Molecule &>(mol).mol.GetTotalSpinMultiplicity();
  } catch (...) {
    return 0;
  }
}
void mol_center(Molecule &mol) {
  try {
    mol.mol.Center();
  } catch (...) {
  }
}
double mol_angle(const Molecule &mol, uint32_t i, uint32_t j, uint32_t k) {
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(i));
    OpenBabel::OBAtom *b = m.GetAtom(static_cast<int>(j));
    OpenBabel::OBAtom *c = m.GetAtom(static_cast<int>(k));
    if (!a || !b || !c) return 0.0;
    return m.GetAngle(a, b, c);
  } catch (...) {
    return 0.0;
  }
}
double mol_torsion(const Molecule &mol, uint32_t i, uint32_t j, uint32_t k,
                   uint32_t l) {
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(i));
    OpenBabel::OBAtom *b = m.GetAtom(static_cast<int>(j));
    OpenBabel::OBAtom *c = m.GetAtom(static_cast<int>(k));
    OpenBabel::OBAtom *d = m.GetAtom(static_cast<int>(l));
    if (!a || !b || !c || !d) return 0.0;
    return OpenBabel::CalcTorsionAngle(a->GetVector(), b->GetVector(),
                                       c->GetVector(), d->GetVector());
  } catch (...) {
    return 0.0;
  }
}
std::unique_ptr<Molecule> mol_clone(const Molecule &mol) {
  auto m = std::unique_ptr<Molecule>(new Molecule());
  try {
    m->mol = const_cast<Molecule &>(mol).mol;  // OBMol copy assignment
  } catch (...) {
    // Copying an OBMol only fails by exhausting memory. Return a valid (empty)
    // molecule rather than let the exception cross the noexcept FFI boundary:
    // the safe wrapper relies on this never being null.
    m->mol.Clear();
  }
  return m;
}
bool mol_strip_salts(Molecule &mol, uint32_t threshold) {
  try {
    return mol.mol.StripSalts(threshold);
  } catch (...) {
    return false;
  }
}
std::unique_ptr<std::vector<Molecule>> mol_separate(const Molecule &mol) {
  auto out = std::unique_ptr<std::vector<Molecule>>(new std::vector<Molecule>());
  try {
    std::vector<OpenBabel::OBMol> frags =
        const_cast<Molecule &>(mol).mol.Separate();
    out->reserve(frags.size());
    for (auto &f : frags) {
      Molecule wrapped;
      wrapped.mol = f;
      out->push_back(std::move(wrapped));
    }
  } catch (...) {
  }
  return out;
}
void mol_set_property(Molecule &mol, rust::Str key, rust::Str value) {
  try {
    std::string k = to_std(key);
    if (mol.mol.HasData(k)) mol.mol.DeleteData(k);
    OpenBabel::OBPairData *pd = new OpenBabel::OBPairData();
    pd->SetAttribute(k);
    pd->SetValue(to_std(value));
    mol.mol.SetData(pd);
  } catch (...) {
  }
}
rust::String mol_get_property(const Molecule &mol, rust::Str key, bool &ok) {
  try {
    OpenBabel::OBGenericData *d =
        const_cast<Molecule &>(mol).mol.GetData(to_std(key));
    OpenBabel::OBPairData *pd = dynamic_cast<OpenBabel::OBPairData *>(d);
    if (!pd) {
      ok = false;
      return rust::String();
    }
    ok = true;
    return to_rust(pd->GetValue());
  } catch (...) {
    ok = false;
    return rust::String();
  }
}

// --- Residues ------------------------------------------------------------

namespace {
// Fetch residue `ridx` (0-based). Returns nullptr if out of range.
OpenBabel::OBResidue *residue_at(const Molecule &mol, uint32_t ridx) {
  OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
  if (ridx >= m.NumResidues()) return nullptr;
  return m.GetResidue(static_cast<int>(ridx));
}
// Convert a single char to a rust::String, treating '\0' and ' ' as empty.
rust::String char_to_string(char c) {
  if (c == '\0' || c == ' ') return rust::String();
  return to_rust(std::string(1, c));
}
}  // namespace

uint32_t mol_num_residues(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.NumResidues();
}
rust::String residue_name(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? to_rust(r->GetName()) : rust::String();
}
int residue_number(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? r->GetNum() : 0;
}
rust::String residue_number_string(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? to_rust(r->GetNumString()) : rust::String();
}
rust::String residue_chain(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? char_to_string(r->GetChain()) : rust::String();
}
rust::String residue_insertion_code(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? char_to_string(r->GetInsertionCode()) : rust::String();
}
uint32_t residue_num_atoms(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? r->GetNumAtoms() : 0;
}
uint32_t residue_num_heavy_atoms(const Molecule &mol, uint32_t ridx) {
  auto *r = residue_at(mol, ridx);
  return r ? r->GetNumHvyAtoms() : 0;
}
rust::Vec<uint32_t> residue_atom_indices(const Molecule &mol, uint32_t ridx) {
  rust::Vec<uint32_t> out;
  auto *r = residue_at(mol, ridx);
  if (!r) return out;
  for (OpenBabel::OBAtom *a : r->GetAtoms())
    if (a) out.push_back(static_cast<uint32_t>(a->GetIdx()) - 1);  // 0-based
  return out;
}
// The atom_residue_* helpers deliberately use HasResidue() (a plain null check)
// rather than GetResidue() alone: OBAtom::GetResidue() force-perceives chains
// when none exist, which would synthesize a fake residue for a molecule parsed
// from SMILES. Guarding on HasResidue() reports residue data only when the
// input format actually carried it (PDB, mmCIF, …), keeping the result
// consistent with mol_num_residues().
int atom_residue_index(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (!a || !a->HasResidue()) return -1;
  OpenBabel::OBResidue *r = a->GetResidue();
  return r ? static_cast<int>(r->GetIdx()) : -1;
}
rust::String atom_residue_atom_id(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (!a || !a->HasResidue()) return rust::String();
  OpenBabel::OBResidue *r = a->GetResidue();
  return r ? to_rust(r->GetAtomID(a)) : rust::String();
}
bool atom_is_hetatm(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (!a || !a->HasResidue()) return false;
  OpenBabel::OBResidue *r = a->GetResidue();
  return r ? r->IsHetAtom(a) : false;
}
uint32_t atom_serial_number(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (!a || !a->HasResidue()) return 0;
  OpenBabel::OBResidue *r = a->GetResidue();
  return r ? r->GetSerialNum(a) : 0;
}

// --- Spectra -------------------------------------------------------------

rust::Vec<double> mol_spectrophore(const Molecule &mol) {
  rust::Vec<double> out;
  try {
    OpenBabel::OBSpectrophore engine;
    std::vector<double> v =
        engine.GetSpectrophore(&const_cast<Molecule &>(mol).mol);
    for (double x : v) out.push_back(x);
  } catch (...) {
  }
  return out;
}

namespace {
const OpenBabel::OBVibrationData *vibration_data(const Molecule &mol) {
  OpenBabel::OBGenericData *d =
      const_cast<Molecule &>(mol).mol.GetData(OpenBabel::OBGenericDataType::VibrationData);
  return dynamic_cast<OpenBabel::OBVibrationData *>(d);
}
}  // namespace

rust::Vec<double> mol_vibration_frequencies(const Molecule &mol) {
  rust::Vec<double> out;
  try {
    const auto *vd = vibration_data(mol);
    if (vd)
      for (double f : vd->GetFrequencies()) out.push_back(f);
  } catch (...) {
  }
  return out;
}
rust::Vec<double> mol_vibration_intensities(const Molecule &mol) {
  rust::Vec<double> out;
  try {
    const auto *vd = vibration_data(mol);
    if (vd)
      for (double x : vd->GetIntensities()) out.push_back(x);
  } catch (...) {
  }
  return out;
}

// --- Bulk coordinates ----------------------------------------------------

rust::Vec<double> mol_coordinates(const Molecule &mol) {
  rust::Vec<double> out;
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    unsigned int n = m.NumAtoms();
    out.reserve(n * 3);
    for (unsigned int i = 1; i <= n; ++i) {  // OBMol atoms are 1-based
      OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(i));
      if (!a) continue;
      out.push_back(a->GetX());
      out.push_back(a->GetY());
      out.push_back(a->GetZ());
    }
  } catch (...) {
  }
  return out;
}

// --- Force-field constraints ---------------------------------------------
// OBFFConstraints uses 1-based atom indices; the Rust API is 0-based, so add 1.

std::unique_ptr<Constraints> constraints_new() {
  return std::unique_ptr<Constraints>(new Constraints());
}
void constraints_add_ignore(Constraints &c, uint32_t atom) {
  c.c.AddIgnore(static_cast<int>(atom) + 1);
}
void constraints_add_atom(Constraints &c, uint32_t atom) {
  c.c.AddAtomConstraint(static_cast<int>(atom) + 1);
}
void constraints_add_atom_x(Constraints &c, uint32_t atom) {
  c.c.AddAtomXConstraint(static_cast<int>(atom) + 1);
}
void constraints_add_atom_y(Constraints &c, uint32_t atom) {
  c.c.AddAtomYConstraint(static_cast<int>(atom) + 1);
}
void constraints_add_atom_z(Constraints &c, uint32_t atom) {
  c.c.AddAtomZConstraint(static_cast<int>(atom) + 1);
}
void constraints_add_distance(Constraints &c, uint32_t a, uint32_t b, double length) {
  c.c.AddDistanceConstraint(static_cast<int>(a) + 1, static_cast<int>(b) + 1, length);
}
void constraints_add_angle(Constraints &c, uint32_t a, uint32_t b, uint32_t d, double angle) {
  c.c.AddAngleConstraint(static_cast<int>(a) + 1, static_cast<int>(b) + 1,
                         static_cast<int>(d) + 1, angle);
}
void constraints_add_torsion(Constraints &c, uint32_t a, uint32_t b, uint32_t d,
                             uint32_t e, double torsion) {
  c.c.AddTorsionConstraint(static_cast<int>(a) + 1, static_cast<int>(b) + 1,
                           static_cast<int>(d) + 1, static_cast<int>(e) + 1, torsion);
}
void constraints_set_factor(Constraints &c, double factor) {
  c.c.SetFactor(factor);
}

// --- Geometry optimization -----------------------------------------------

namespace {
// Set up the force-field plugin singleton over `mol` with `constraints`. Returns
// the singleton (NOT owned — never delete it), or nullptr on unknown force field
// / setup failure.
//
// We use the shared singleton, matching `mol_optimize` / `mol_energy`: the safe
// wrapper serializes every OpenBabel call and runs each whole optimization as one
// atomic call, so the singleton's per-run state is used entirely within one lock
// hold. (Private MakeNewInstance() instances are avoided here — creating and
// destroying them concurrently with molecule teardown has proven to corrupt the
// heap.)
//
// Constraints are installed with an explicit SetConstraints() AFTER Setup(),
// rather than via the two-argument Setup(mol, constraints), on purpose:
// OBForceField::Setup caches its previous setup and, when IsSetupNeeded() sees a
// topologically identical molecule, takes a fast path that silently keeps the
// PREVIOUS constraints. SetConstraints() unconditionally installs ours (and
// re-runs the term setup when the ignore set changes), so restraints are applied
// correctly no matter what was optimized just before.
OpenBabel::OBForceField *setup_forcefield(OpenBabel::OBMol &mol, const std::string &ff_id,
                                          OpenBabel::OBFFConstraints &constraints) {
  OpenBabel::OBForceField *ff = find_forcefield(ff_id);
  if (!ff) return nullptr;
  if (!ff->Setup(mol)) return nullptr;
  ff->SetConstraints(constraints);
  return ff;
}
void ff_initialize(OpenBabel::OBForceField *ff, uint32_t algorithm, int steps, double econv) {
  switch (algorithm) {
    case 1:
      ff->ConjugateGradientsInitialize(steps, econv);
      break;
    case 2:
      ff->LBFGSInitialize(steps, econv);
      break;
    case 0:
    default:
      ff->SteepestDescentInitialize(steps, econv);
      break;
  }
}
bool ff_take_steps(OpenBabel::OBForceField *ff, uint32_t algorithm, int n) {
  switch (algorithm) {
    case 1:
      return ff->ConjugateGradientsTakeNSteps(n);
    case 2:
      return ff->LBFGSTakeNSteps(n);
    case 0:
    default:
      return ff->SteepestDescentTakeNSteps(n);
  }
}
// Clear the shared static constraint set so later unconstrained force-field
// calls (mol_energy / mol_optimize) are not tripped up by our restraints. The
// force field is the shared singleton and must NOT be deleted.
void ff_finish(OpenBabel::OBForceField *ff) {
  if (!ff) return;
  try {
    OpenBabel::OBFFConstraints empty;
    ff->SetConstraints(empty);
  } catch (...) {
  }
}
}  // namespace

double optimizer_run_to_end(Molecule &mol, rust::Str ff_id, uint32_t algorithm,
                            uint32_t steps, double econv, const Constraints &constraints,
                            bool &ok) {
  OpenBabel::OBForceField *ff = nullptr;
  try {
    ff = setup_forcefield(mol.mol, to_std(ff_id), const_cast<Constraints &>(constraints).c);
    if (!ff) {
      ok = false;
      return std::nan("");
    }
    int n = static_cast<int>(steps);
    ff_initialize(ff, algorithm, n, econv);
    while (ff_take_steps(ff, algorithm, n)) {
    }
    ff->GetCoordinates(mol.mol);
    double e = ff->Energy(false);
    ff_finish(ff);
    ok = true;
    return e;
  } catch (...) {
    ff_finish(ff);
    ok = false;
    return std::nan("");
  }
}

rust::Vec<double> optimizer_run_trajectory(Molecule &mol, rust::Str ff_id, uint32_t algorithm,
                                           uint32_t steps, double econv,
                                           const Constraints &constraints,
                                           uint32_t frame_interval) {
  rust::Vec<double> out;
  OpenBabel::OBForceField *ff = nullptr;
  try {
    ff = setup_forcefield(mol.mol, to_std(ff_id), const_cast<Constraints &>(constraints).c);
    if (!ff) return out;
    int budget = static_cast<int>(steps);
    int chunk = frame_interval == 0 ? 1 : static_cast<int>(frame_interval);
    ff_initialize(ff, algorithm, budget, econv);
    bool more = true;
    while (more) {
      more = ff_take_steps(ff, algorithm, chunk);
      ff->GetCoordinates(mol.mol);
      out.push_back(ff->Energy(false));  // frame = [energy, x0,y0,z0, ...]
      unsigned int natoms = mol.mol.NumAtoms();
      for (unsigned int i = 1; i <= natoms; ++i) {
        OpenBabel::OBAtom *a = mol.mol.GetAtom(static_cast<int>(i));
        if (!a) continue;
        out.push_back(a->GetX());
        out.push_back(a->GetY());
        out.push_back(a->GetZ());
      }
    }
    ff_finish(ff);
  } catch (...) {
    ff_finish(ff);
  }
  return out;
}

// --- Molecule construction & editing --------------------------------------

uint32_t mol_add_atom(Molecule &mol, uint32_t atomic_num) {
  OpenBabel::OBAtom *a = mol.mol.NewAtom();
  a->SetAtomicNum(static_cast<int>(atomic_num));
  return static_cast<uint32_t>(a->GetIdx()) - 1;  // 0-based
}
bool mol_add_bond(Molecule &mol, uint32_t begin, uint32_t end, uint32_t order) {
  try {
    // OpenBabel bond endpoints are 1-based atom indices.
    return mol.mol.AddBond(static_cast<int>(begin) + 1, static_cast<int>(end) + 1,
                           static_cast<int>(order));
  } catch (...) {
    return false;
  }
}
bool mol_delete_atom(Molecule &mol, uint32_t idx) {
  OpenBabel::OBAtom *a = mol.mol.GetAtom(static_cast<int>(idx) + 1);
  if (!a) return false;
  return mol.mol.DeleteAtom(a);
}
bool mol_delete_bond(Molecule &mol, uint32_t idx) {
  OpenBabel::OBBond *b = mol.mol.GetBond(static_cast<int>(idx));
  if (!b) return false;
  return mol.mol.DeleteBond(b);
}
void mol_begin_modify(Molecule &mol) {
  try {
    mol.mol.BeginModify();
  } catch (...) {
  }
}
void mol_end_modify(Molecule &mol) {
  // EndModify re-runs perception (rings, aromaticity, …), which can throw.
  try {
    mol.mol.EndModify();
  } catch (...) {
  }
}
void mol_clear(Molecule &mol) { mol.mol.Clear(); }
void mol_translate(Molecule &mol, double x, double y, double z) {
  mol.mol.Translate(OpenBabel::vector3(x, y, z));
}
bool mol_set_coordinates(Molecule &mol, rust::Slice<const double> coords) {
  try {
    unsigned int n = mol.mol.NumAtoms();
    if (coords.size() != static_cast<size_t>(n) * 3) return false;
    if (n == 0) return true;
    // OBMol::SetCoordinates copies into the existing conformer buffer when one
    // is present, but ADOPTS the passed pointer (stores it in _vconf) when the
    // molecule has no coordinates yet. So allocate on the heap and only free
    // our buffer in the copy case; otherwise OpenBabel now owns it.
    bool had_coords = mol.mol.GetCoordinates() != nullptr;
    double *buf = new double[coords.size()];
    std::copy(coords.begin(), coords.end(), buf);
    mol.mol.SetCoordinates(buf);
    if (had_coords) delete[] buf;
    return true;
  } catch (...) {
    return false;
  }
}
void mol_set_dimension(Molecule &mol, uint32_t dim) {
  mol.mol.SetDimension(static_cast<unsigned short>(dim));
}
void mol_connect_the_dots(Molecule &mol) {
  try {
    mol.mol.ConnectTheDots();
  } catch (...) {
  }
}
void mol_perceive_bond_orders(Molecule &mol) {
  try {
    mol.mol.PerceiveBondOrders();
  } catch (...) {
  }
}
bool mol_add_polar_hydrogens(Molecule &mol) {
  try {
    return mol.mol.AddPolarHydrogens();
  } catch (...) {
    return false;
  }
}
bool mol_convert_dative_bonds(Molecule &mol) {
  try {
    return mol.mol.ConvertDativeBonds();
  } catch (...) {
    return false;
  }
}
bool mol_assign_spin_multiplicity(Molecule &mol) {
  try {
    return mol.mol.AssignSpinMultiplicity();
  } catch (...) {
    return false;
  }
}
bool mol_add_hydrogens_ph(Molecule &mol, double ph) {
  try {
    return mol.mol.AddHydrogens(false, true, ph);
  } catch (...) {
    return false;
  }
}

// --- Atom setters (idx is 1-based) ----------------------------------------

void atom_set_atomic_num(Molecule &mol, uint32_t idx, uint32_t atomic_num) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetAtomicNum(static_cast<int>(atomic_num));
}
void atom_set_formal_charge(Molecule &mol, uint32_t idx, int charge) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetFormalCharge(charge);
}
void atom_set_position(Molecule &mol, uint32_t idx, double x, double y, double z) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetVector(x, y, z);
}
void atom_set_isotope(Molecule &mol, uint32_t idx, uint32_t isotope) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetIsotope(isotope);
}
void atom_set_spin_multiplicity(Molecule &mol, uint32_t idx, int spin) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetSpinMultiplicity(static_cast<short>(spin));
}
void atom_set_partial_charge(Molecule &mol, uint32_t idx, double charge) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetPartialCharge(charge);
}
void atom_set_type(Molecule &mol, uint32_t idx, rust::Str type_name) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetType(to_std(type_name));
}
void atom_set_implicit_h(Molecule &mol, uint32_t idx, uint32_t count) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetImplicitHCount(count);
}

// --- Bond setters (idx is 0-based) ----------------------------------------

void bond_set_order(Molecule &mol, uint32_t idx, uint32_t order) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  if (b) b->SetBondOrder(static_cast<int>(order));
}
bool bond_set_length(Molecule &mol, uint32_t idx, double length) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  if (!b) return false;
  b->SetLength(length);
  return true;
}

// --- Multi-molecule input -------------------------------------------------

std::unique_ptr<std::vector<Molecule>> mol_read_many(rust::Str format, rust::Str data) {
  auto out = std::unique_ptr<std::vector<Molecule>>(new std::vector<Molecule>());
  try {
    OpenBabel::OBConversion conv;
    if (!conv.SetInFormat(to_std(format).c_str())) return out;
    std::istringstream iss(to_std(data));
    conv.SetInStream(&iss, false);
    while (true) {
      Molecule m;
      if (!conv.Read(&m.mol)) break;
      out->push_back(std::move(m));
    }
  } catch (...) {
  }
  return out;
}

// --- File I/O (format auto-detected from extension when `format` empty) ----

namespace {
// Resolve the input format from `format` (if given) or the file extension.
bool set_in_format(OpenBabel::OBConversion &conv, const std::string &format,
                   const std::string &path) {
  if (!format.empty()) return conv.SetInFormat(format.c_str());
  OpenBabel::OBFormat *f = OpenBabel::OBConversion::FormatFromExt(path);
  return f && conv.SetInFormat(f);
}
bool set_out_format(OpenBabel::OBConversion &conv, const std::string &format,
                    const std::string &path) {
  if (!format.empty()) return conv.SetOutFormat(format.c_str());
  OpenBabel::OBFormat *f = OpenBabel::OBConversion::FormatFromExt(path);
  return f && conv.SetOutFormat(f);
}
}  // namespace

std::unique_ptr<Molecule> mol_read_file(rust::Str path, rust::Str format) {
  try {
    OpenBabel::OBConversion conv;
    std::string p = to_std(path);
    if (!set_in_format(conv, to_std(format), p)) return nullptr;
    auto m = std::unique_ptr<Molecule>(new Molecule());
    if (!conv.ReadFile(&m->mol, p)) return nullptr;
    return m;
  } catch (...) {
    return nullptr;
  }
}

std::unique_ptr<std::vector<Molecule>> mol_read_file_many(rust::Str path, rust::Str format) {
  auto out = std::unique_ptr<std::vector<Molecule>>(new std::vector<Molecule>());
  try {
    OpenBabel::OBConversion conv;
    std::string p = to_std(path);
    if (!set_in_format(conv, to_std(format), p)) return out;
    Molecule first;
    // ReadFile opens the file, sets it as the input stream, and reads record #1.
    if (!conv.ReadFile(&first.mol, p)) return out;
    out->push_back(std::move(first));
    while (true) {
      Molecule m;
      if (!conv.Read(&m.mol)) break;  // continue on the same open stream
      out->push_back(std::move(m));
    }
  } catch (...) {
  }
  return out;
}

void mol_write_file(const Molecule &mol, rust::Str path, rust::Str format, bool &ok) {
  try {
    OpenBabel::OBConversion conv;
    std::string p = to_std(path);
    if (!set_out_format(conv, to_std(format), p)) {
      ok = false;
      return;
    }
    ok = conv.WriteFile(&const_cast<Molecule &>(mol).mol, p);
  } catch (...) {
    ok = false;
  }
}

// --- Ring access (SSSR; ring_idx is 0-based) ------------------------------

namespace {
OpenBabel::OBRing *ring_at(const Molecule &mol, uint32_t ring_idx) {
  std::vector<OpenBabel::OBRing *> &rings = const_cast<Molecule &>(mol).mol.GetSSSR();
  if (ring_idx >= rings.size()) return nullptr;
  return rings[ring_idx];
}
}  // namespace

// ring_at() (and hence every accessor below) calls OBMol::GetSSSR(), which runs
// ring perception on first use and can throw; guard each like mol_num_rings.
uint32_t ring_size(const Molecule &mol, uint32_t ring_idx) {
  try {
    OpenBabel::OBRing *r = ring_at(mol, ring_idx);
    return r ? static_cast<uint32_t>(r->Size()) : 0;
  } catch (...) {
    return 0;
  }
}
rust::Vec<uint32_t> ring_atom_indices(const Molecule &mol, uint32_t ring_idx) {
  rust::Vec<uint32_t> out;
  try {
    OpenBabel::OBRing *r = ring_at(mol, ring_idx);
    if (!r) return out;
    for (int idx : r->_path) out.push_back(static_cast<uint32_t>(idx) - 1);  // 0-based
  } catch (...) {
    out.clear();
  }
  return out;
}
bool ring_is_aromatic(const Molecule &mol, uint32_t ring_idx) {
  try {
    OpenBabel::OBRing *r = ring_at(mol, ring_idx);
    return r ? r->IsAromatic() : false;
  } catch (...) {
    return false;
  }
}

// --- Graph navigation -----------------------------------------------------

rust::Vec<uint32_t> atom_neighbor_indices(const Molecule &mol, uint32_t idx) {
  rust::Vec<uint32_t> out;
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (!a) return out;
  OpenBabel::OBBondIterator it;
  for (OpenBabel::OBAtom *n = a->BeginNbrAtom(it); n; n = a->NextNbrAtom(it))
    out.push_back(static_cast<uint32_t>(n->GetIdx()) - 1);
  return out;
}
rust::Vec<uint32_t> atom_bond_indices(const Molecule &mol, uint32_t idx) {
  rust::Vec<uint32_t> out;
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (!a) return out;
  OpenBabel::OBBondIterator it;
  for (OpenBabel::OBBond *b = a->BeginBond(it); b; b = a->NextBond(it))
    out.push_back(static_cast<uint32_t>(b->GetIdx()));  // bond idx already 0-based
  return out;
}
uint32_t atom_count_bonds_of_order(const Molecule &mol, uint32_t idx, uint32_t order) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->CountBondsOfOrder(order) : 0;
}
uint32_t atom_explicit_h_count(const Molecule &mol, uint32_t idx) {
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? a->ExplicitHydrogenCount() : 0;
}
int mol_bond_between(const Molecule &mol, uint32_t a, uint32_t b) {
  // OBMol::GetBond takes 1-based atom indices.
  OpenBabel::OBBond *bond =
      const_cast<Molecule &>(mol).mol.GetBond(static_cast<int>(a) + 1, static_cast<int>(b) + 1);
  return bond ? static_cast<int>(bond->GetIdx()) : -1;
}
int bond_other_atom(const Molecule &mol, uint32_t bond_idx, uint32_t atom_idx) {
  auto *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, bond_idx));
  auto *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, atom_idx + 1));  // atom_idx is 0-based here
  if (!b || !a) return -1;
  OpenBabel::OBAtom *other = b->GetNbrAtom(a);
  return other ? static_cast<int>(other->GetIdx()) - 1 : -1;
}

// --- Crystallography (unit cell) ------------------------------------------

namespace {
OpenBabel::OBUnitCell *unit_cell(const Molecule &mol) {
  OpenBabel::OBGenericData *d =
      const_cast<Molecule &>(mol).mol.GetData(OpenBabel::OBGenericDataType::UnitCell);
  return dynamic_cast<OpenBabel::OBUnitCell *>(d);
}
}  // namespace

bool mol_has_unit_cell(const Molecule &mol) { return unit_cell(mol) != nullptr; }
rust::Vec<double> mol_cell_parameters(const Molecule &mol) {
  rust::Vec<double> out;
  OpenBabel::OBUnitCell *c = unit_cell(mol);
  if (!c) return out;
  out.push_back(c->GetA());
  out.push_back(c->GetB());
  out.push_back(c->GetC());
  out.push_back(c->GetAlpha());
  out.push_back(c->GetBeta());
  out.push_back(c->GetGamma());
  return out;
}
double mol_cell_volume(const Molecule &mol) {
  OpenBabel::OBUnitCell *c = unit_cell(mol);
  return c ? c->GetCellVolume() : 0.0;
}
rust::String mol_cell_spacegroup(const Molecule &mol) {
  OpenBabel::OBUnitCell *c = unit_cell(mol);
  return c ? to_rust(c->GetSpaceGroupName()) : rust::String();
}
uint32_t mol_cell_lattice_type(const Molecule &mol) {
  OpenBabel::OBUnitCell *c = unit_cell(mol);
  return c ? static_cast<uint32_t>(c->GetLatticeType()) : 0;
}
rust::Vec<double> mol_cell_to_fractional(const Molecule &mol, double x, double y, double z) {
  rust::Vec<double> out;
  OpenBabel::OBUnitCell *c = unit_cell(mol);
  if (!c) return out;
  OpenBabel::vector3 f = c->CartesianToFractional(OpenBabel::vector3(x, y, z));
  out.push_back(f.x());
  out.push_back(f.y());
  out.push_back(f.z());
  return out;
}
rust::Vec<double> mol_cell_to_cartesian(const Molecule &mol, double x, double y, double z) {
  rust::Vec<double> out;
  OpenBabel::OBUnitCell *c = unit_cell(mol);
  if (!c) return out;
  OpenBabel::vector3 v = c->FractionalToCartesian(OpenBabel::vector3(x, y, z));
  out.push_back(v.x());
  out.push_back(v.y());
  out.push_back(v.z());
  return out;
}

// --- Symmetry & canonical ordering ----------------------------------------

// Topological symmetry class per atom, indexed in atom order (0-based); atoms
// sharing a value are graph-equivalent. Vector length == number of atoms.
rust::Vec<uint32_t> mol_symmetry_classes(const Molecule &mol) {
  rust::Vec<uint32_t> out;
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    OpenBabel::OBGraphSym gs(&m);
    std::vector<unsigned int> sym;
    gs.GetSymmetry(sym);
    for (unsigned int v : sym) out.push_back(static_cast<uint32_t>(v));
  } catch (...) {
    out.clear();  // never a partial (length must equal num_atoms, or be empty)
  }
  return out;
}

// Canonical rank (1-based) per atom, indexed in atom order (0-based) — a
// repeatable canonical labelling built from the symmetry classes.
rust::Vec<uint32_t> mol_canonical_ranks(const Molecule &mol) {
  rust::Vec<uint32_t> out;
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    OpenBabel::OBGraphSym gs(&m);
    std::vector<unsigned int> sym;
    gs.GetSymmetry(sym);
    std::vector<unsigned int> canon;
    // Empty mask == all atoms (see canon.cpp); default maxSeconds, not onlyOne.
    OpenBabel::CanonicalLabels(&m, sym, canon);
    for (unsigned int v : canon) out.push_back(static_cast<uint32_t>(v));
  } catch (...) {
    out.clear();  // never a partial (length must equal num_atoms, or be empty)
  }
  return out;
}

// --- Reactions ------------------------------------------------------------

std::unique_ptr<Reaction> reaction_new() {
  return std::unique_ptr<Reaction>(new Reaction());
}

std::unique_ptr<Reaction> reaction_read(rust::Str format, rust::Str data) {
  try {
    OpenBabel::OBConversion conv;
    if (!conv.SetInFormat(to_std(format).c_str())) return nullptr;
    auto r = std::unique_ptr<Reaction>(new Reaction());
    std::istringstream iss(to_std(data));
    if (!conv.Read(&r->rxn, &iss)) return nullptr;
    return r;
  } catch (...) {
    return nullptr;
  }
}

rust::String reaction_write(const Reaction &r, rust::Str format, bool &ok) {
  try {
    OpenBabel::OBConversion conv;
    if (!conv.SetOutFormat(to_std(format).c_str())) {
      ok = false;
      return rust::String();
    }
    std::ostringstream oss;
    ok = conv.Write(&const_cast<Reaction &>(r).rxn, &oss);
    return to_rust(oss.str());
  } catch (...) {
    ok = false;
    return rust::String();
  }
}

uint32_t reaction_num_reactants(const Reaction &r) {
  return static_cast<uint32_t>(const_cast<Reaction &>(r).rxn.NumReactants());
}
uint32_t reaction_num_products(const Reaction &r) {
  return static_cast<uint32_t>(const_cast<Reaction &>(r).rxn.NumProducts());
}
uint32_t reaction_num_agents(const Reaction &r) {
  return static_cast<uint32_t>(const_cast<Reaction &>(r).rxn.NumAgents());
}

// Copy a shared_ptr<OBMol> component out as a standalone Molecule.
static std::unique_ptr<Molecule> component_copy(std::shared_ptr<OpenBabel::OBMol> sp) {
  if (!sp) return nullptr;
  try {
    auto m = std::unique_ptr<Molecule>(new Molecule());
    m->mol = *sp;  // OBMol copy assignment
    return m;
  } catch (...) {
    return nullptr;  // out-of-memory: the safe wrapper maps null -> None
  }
}

std::unique_ptr<Molecule> reaction_reactant(const Reaction &r, uint32_t i) {
  return component_copy(const_cast<Reaction &>(r).rxn.GetReactant(i));
}
std::unique_ptr<Molecule> reaction_product(const Reaction &r, uint32_t i) {
  return component_copy(const_cast<Reaction &>(r).rxn.GetProduct(i));
}
std::unique_ptr<Molecule> reaction_agent(const Reaction &r, uint32_t i) {
  return component_copy(const_cast<Reaction &>(r).rxn.GetAgent(i));
}

void reaction_add_reactant(Reaction &r, const Molecule &mol) {
  try {
    r.rxn.AddReactant(std::make_shared<OpenBabel::OBMol>(const_cast<Molecule &>(mol).mol));
  } catch (...) {
  }
}
void reaction_add_product(Reaction &r, const Molecule &mol) {
  try {
    r.rxn.AddProduct(std::make_shared<OpenBabel::OBMol>(const_cast<Molecule &>(mol).mol));
  } catch (...) {
  }
}
void reaction_add_agent(Reaction &r, const Molecule &mol) {
  try {
    r.rxn.AddAgent(std::make_shared<OpenBabel::OBMol>(const_cast<Molecule &>(mol).mol));
  } catch (...) {
  }
}

rust::String reaction_title(const Reaction &r) {
  return to_rust(const_cast<Reaction &>(r).rxn.GetTitle());
}
void reaction_set_title(Reaction &r, rust::Str title) {
  r.rxn.SetTitle(to_std(title));
}
rust::String reaction_comment(const Reaction &r) {
  return to_rust(const_cast<Reaction &>(r).rxn.GetComment());
}
void reaction_set_comment(Reaction &r, rust::Str comment) {
  r.rxn.SetComment(to_std(comment));
}
bool reaction_is_reversible(const Reaction &r) {
  return const_cast<Reaction &>(r).rxn.IsReversible();
}
void reaction_set_reversible(Reaction &r, bool value) {
  r.rxn.SetReversible(value);
}

// --- Subgraph isomorphism & automorphisms ---------------------------------

// Flatten a set of mappings into `width` queried-atom indices per mapping,
// ordered by query-atom index (both 0-based). `pairs` are (query, queried).
static void flatten_mappings(const OpenBabel::OBIsomorphismMapper::Mappings &maps,
                             uint32_t width, rust::Vec<uint32_t> &out) {
  for (const auto &m : maps) {
    std::vector<uint32_t> row(width, 0);
    for (const auto &pr : m)
      if (pr.first < width) row[pr.first] = static_cast<uint32_t>(pr.second);
    for (uint32_t v : row) out.push_back(v);
  }
}

// All unique mappings of `query` as a substructure of `target`. Sets `width`
// to the number of query atoms; the flat result holds `width` target atom
// indices (0-based) per mapping. Empty if there is no match.
rust::Vec<uint32_t> mol_substructure_mappings(const Molecule &query, const Molecule &target,
                                              uint32_t &width) {
  rust::Vec<uint32_t> out;
  width = 0;
  try {
    OpenBabel::OBMol &q = const_cast<Molecule &>(query).mol;
    OpenBabel::OBMol &t = const_cast<Molecule &>(target).mol;
    // Own the query/mapper via unique_ptr so they are freed on every path,
    // including if MapUnique or flatten_mappings throws — the raw new/delete
    // version leaked both on the exception path.
    std::unique_ptr<OpenBabel::OBQuery> oq(OpenBabel::CompileMoleculeQuery(&q));
    if (!oq) return out;
    width = static_cast<uint32_t>(oq->GetAtoms().size());
    std::unique_ptr<OpenBabel::OBIsomorphismMapper> mapper(
        OpenBabel::OBIsomorphismMapper::GetInstance(oq.get()));
    if (!mapper) return out;
    OpenBabel::OBIsomorphismMapper::Mappings maps;
    mapper->MapUnique(&t, maps);
    flatten_mappings(maps, width, out);
  } catch (...) {
    out.clear();
  }
  return out;
}

// All graph automorphisms of `mol`; sets `width` to the atom count. Each
// automorphism is a permutation of atom indices (0-based), `width` per row.
rust::Vec<uint32_t> mol_automorphisms(const Molecule &mol, uint32_t &width) {
  rust::Vec<uint32_t> out;
  width = 0;
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    width = static_cast<uint32_t>(m.NumAtoms());
    OpenBabel::OBIsomorphismMapper::Mappings auts;
    OpenBabel::FindAutomorphisms(&m, auts);
    flatten_mappings(auts, width, out);
  } catch (...) {
    out.clear();
  }
  return out;
}

// --- Geometry & topology (niche) ------------------------------------------

// Set the a-b-c-d torsion to `radians`, rotating the b-c bond's far side.
// Atom indices are 0-based; a no-op if any is out of range.
void mol_set_torsion(Molecule &mol, uint32_t a, uint32_t b, uint32_t c, uint32_t d,
                     double radians) {
  try {
    auto *pa = const_cast<OpenBabel::OBAtom *>(atom_at(mol, a + 1));
    auto *pb = const_cast<OpenBabel::OBAtom *>(atom_at(mol, b + 1));
    auto *pc = const_cast<OpenBabel::OBAtom *>(atom_at(mol, c + 1));
    auto *pd = const_cast<OpenBabel::OBAtom *>(atom_at(mol, d + 1));
    if (pa && pb && pc && pd) mol.mol.SetTorsion(pa, pb, pc, pd, radians);
  } catch (...) {
  }
}

// 0-based indices of the atoms reachable from `to` without passing back
// through `from` (excludes both endpoints).
rust::Vec<uint32_t> mol_find_children(const Molecule &mol, uint32_t from, uint32_t to) {
  rust::Vec<uint32_t> out;
  std::vector<int> children;  // 1-based
  const_cast<Molecule &>(mol).mol.FindChildren(children, static_cast<int>(from) + 1,
                                               static_cast<int>(to) + 1);
  for (int i : children) out.push_back(static_cast<uint32_t>(i) - 1);
  return out;
}

// 0-based indices of the atoms in the largest connected fragment.
rust::Vec<uint32_t> mol_largest_fragment(const Molecule &mol) {
  rust::Vec<uint32_t> out;
  OpenBabel::OBBitVec frag;
  const_cast<Molecule &>(mol).mol.FindLargestFragment(frag);
  for (int i = frag.NextBit(-1); i != frag.EndBit(); i = frag.NextBit(i))
    out.push_back(static_cast<uint32_t>(i) - 1);  // bits are 1-based atom idx
  return out;
}

void mol_set_total_charge(Molecule &mol, int32_t charge) {
  mol.mol.SetTotalCharge(static_cast<int>(charge));
}
void mol_set_total_spin(Molecule &mol, uint32_t spin) {
  mol.mol.SetTotalSpinMultiplicity(static_cast<unsigned int>(spin));
}

// --- Per-atom / per-bond string data (OBPairData) -------------------------

namespace {
// Attach (replacing any existing) a key/value string pair to an OBBase.
void set_pair_data(OpenBabel::OBBase *b, const std::string &key, const std::string &value) {
  if (b->HasData(key)) b->DeleteData(key);
  OpenBabel::OBPairData *pd = new OpenBabel::OBPairData();
  pd->SetAttribute(key);
  pd->SetValue(value);
  b->SetData(pd);
}
// Read a key's string value from an OBBase; ok=false if absent / not a pair.
rust::String get_pair_data(OpenBabel::OBBase *b, const std::string &key, bool &ok) {
  OpenBabel::OBPairData *pd = dynamic_cast<OpenBabel::OBPairData *>(b->GetData(key));
  if (!pd) {
    ok = false;
    return rust::String();
  }
  ok = true;
  return rust::String(pd->GetValue());
}
}  // namespace

void atom_set_data(Molecule &mol, uint32_t idx, rust::Str key, rust::Str value) {
  try {
    OpenBabel::OBAtom *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
    if (a) set_pair_data(a, to_std(key), to_std(value));
  } catch (...) {
  }
}
rust::String atom_get_data(const Molecule &mol, uint32_t idx, rust::Str key, bool &ok) {
  try {
    OpenBabel::OBAtom *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
    if (!a) {
      ok = false;
      return rust::String();
    }
    return get_pair_data(a, to_std(key), ok);
  } catch (...) {
    ok = false;
    return rust::String();
  }
}
void bond_set_data(Molecule &mol, uint32_t idx, rust::Str key, rust::Str value) {
  try {
    OpenBabel::OBBond *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
    if (b) set_pair_data(b, to_std(key), to_std(value));
  } catch (...) {
  }
}
rust::String bond_get_data(const Molecule &mol, uint32_t idx, rust::Str key, bool &ok) {
  try {
    OpenBabel::OBBond *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
    if (!b) {
      ok = false;
      return rust::String();
    }
    return get_pair_data(b, to_std(key), ok);
  } catch (...) {
    ok = false;
    return rust::String();
  }
}

// --- Inter-atom distance & 2D wedge/hash bond stereo ----------------------

// Distance (Å) between atoms `i` and `j` (1-based); 0.0 for invalid indices.
double mol_distance(const Molecule &mol, uint32_t i, uint32_t j) {
  try {
    OpenBabel::OBMol &m = const_cast<Molecule &>(mol).mol;
    OpenBabel::OBAtom *a = m.GetAtom(static_cast<int>(i));
    OpenBabel::OBAtom *b = m.GetAtom(static_cast<int>(j));
    if (!a || !b) return 0.0;
    return a->GetDistance(b);
  } catch (...) {
    return 0.0;
  }
}

bool bond_is_wedge(const Molecule &mol, uint32_t idx) {
  OpenBabel::OBBond *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsWedge() : false;
}
bool bond_is_hash(const Molecule &mol, uint32_t idx) {
  OpenBabel::OBBond *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  return b ? b->IsHash() : false;
}
void bond_set_wedge(Molecule &mol, uint32_t idx, bool value) {
  OpenBabel::OBBond *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  if (b) b->SetWedge(value);
}
void bond_set_hash(Molecule &mol, uint32_t idx, bool value) {
  OpenBabel::OBBond *b = const_cast<OpenBabel::OBBond *>(bond_at(mol, idx));
  if (b) b->SetHash(value);
}

// --- Persistent atom ids, connectivity relations, LSSR --------------------

// Persistent atom id (survives atom deletion, unlike the 1-based index).
uint64_t atom_id(const Molecule &mol, uint32_t idx) {
  const OpenBabel::OBAtom *a = atom_at(mol, idx);
  return a ? static_cast<uint64_t>(a->GetId()) : 0;
}
void atom_set_id(Molecule &mol, uint32_t idx, uint64_t id) {
  OpenBabel::OBAtom *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  if (a) a->SetId(static_cast<unsigned long>(id));
}

// Connectivity relations between atoms `a` and `b` (both 1-based).
bool atom_is_connected(const Molecule &mol, uint32_t a, uint32_t b) {
  OpenBabel::OBAtom *pa = const_cast<OpenBabel::OBAtom *>(atom_at(mol, a));
  OpenBabel::OBAtom *pb = const_cast<OpenBabel::OBAtom *>(atom_at(mol, b));
  return (pa && pb) ? pa->IsConnected(pb) : false;
}
bool atom_is_one_three(const Molecule &mol, uint32_t a, uint32_t b) {
  OpenBabel::OBAtom *pa = const_cast<OpenBabel::OBAtom *>(atom_at(mol, a));
  OpenBabel::OBAtom *pb = const_cast<OpenBabel::OBAtom *>(atom_at(mol, b));
  return (pa && pb) ? pa->IsOneThree(pb) : false;
}
bool atom_is_one_four(const Molecule &mol, uint32_t a, uint32_t b) {
  OpenBabel::OBAtom *pa = const_cast<OpenBabel::OBAtom *>(atom_at(mol, a));
  OpenBabel::OBAtom *pb = const_cast<OpenBabel::OBAtom *>(atom_at(mol, b));
  return (pa && pb) ? pa->IsOneFour(pb) : false;
}

// Ring sizes from the Large Set of Smallest Rings (an alternative to SSSR).
rust::Vec<uint32_t> mol_lssr_sizes(const Molecule &mol) {
  rust::Vec<uint32_t> out;
  std::vector<OpenBabel::OBRing *> &rings = const_cast<Molecule &>(mol).mol.GetLSSR();
  for (OpenBabel::OBRing *r : rings)
    if (r) out.push_back(static_cast<uint32_t>(r->Size()));
  return out;
}

// --- Perception state flags & targeted hydrogen editing -------------------

bool mol_has_aromatic_perceived(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.HasAromaticPerceived();
}
bool mol_has_sssr_perceived(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.HasSSSRPerceived();
}
bool mol_has_ring_atoms_perceived(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.HasRingAtomsAndBondsPerceived();
}
bool mol_has_chains_perceived(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.HasChainsPerceived();
}
bool mol_has_hydrogens_added(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.HasHydrogensAdded();
}
bool mol_has_nonzero_coords(const Molecule &mol) {
  return const_cast<Molecule &>(mol).mol.HasNonZeroCoords();
}

// Add explicit hydrogens to just atom `idx` (1-based); false if out of range.
bool mol_add_hydrogens_to_atom(Molecule &mol, uint32_t idx) {
  OpenBabel::OBAtom *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? mol.mol.AddHydrogens(a) : false;
}
// Remove the explicit hydrogens attached to atom `idx` (1-based).
bool mol_delete_hydrogens_of_atom(Molecule &mol, uint32_t idx) {
  OpenBabel::OBAtom *a = const_cast<OpenBabel::OBAtom *>(atom_at(mol, idx));
  return a ? mol.mol.DeleteHydrogens(a) : false;
}

}  // namespace ob_shim
