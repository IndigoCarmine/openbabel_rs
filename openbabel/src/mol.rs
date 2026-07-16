//! The [`Molecule`] type — a safe, owning wrapper around OpenBabel's `OBMol`.

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::atom::{Atom, AtomMut};
use crate::bond::{Bond, BondMut};
use crate::error::Error;
use crate::residue::Residue;
use crate::ring::Ring;
use crate::with_ob;

/// Split a flat `Vec<u32>` of `width`-sized rows (as returned by the flattened
/// isomorphism shim calls) into a `Vec` of rows. A zero width yields no rows.
fn chunk(flat: Vec<u32>, width: usize) -> Vec<Vec<u32>> {
    if width == 0 {
        return Vec::new();
    }
    flat.chunks(width).map(<[u32]>::to_vec).collect()
}

/// Rendering options for [`Molecule::to_svg_with`].
///
/// The default (`SvgOptions::default()`, used by [`Molecule::to_svg`]) draws a
/// clean skeletal structure: only terminal carbons are labelled and atoms are
/// not indexed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SvgOptions {
    /// Draw a label on every carbon atom, not just terminal ones.
    pub all_carbons: bool,
    /// Annotate each atom with its index.
    pub atom_indices: bool,
}

/// A molecule.
///
/// Construct one by parsing (`Molecule::parse`) or building up from empty
/// (`Molecule::new`). Dropping it frees the underlying OpenBabel object.
pub struct Molecule {
    inner: UniquePtr<ffi::Molecule>,
}

impl Molecule {
    /// Create an empty molecule.
    pub fn new() -> Self {
        Molecule {
            inner: with_ob(ffi::mol_new),
        }
    }

    /// Parse `data` using the OpenBabel format id `format`
    /// (e.g. `"smi"`, `"mol"`, `"sdf"`, `"pdb"`, `"inchi"`).
    ///
    /// Returns [`Error::Parse`] if the data cannot be read in that format.
    pub fn parse(data: &str, format: &str) -> Result<Self, Error> {
        let inner = with_ob(|| ffi::mol_read(format, data));
        if inner.is_null() {
            Err(Error::Parse {
                format: format.to_string(),
            })
        } else {
            Ok(Molecule { inner })
        }
    }

    /// Parse **every** record from `data` in `format` — a multi-record SDF, one
    /// SMILES per line, a multi-model PDB, etc.
    ///
    /// Returns one [`Molecule`] per record, in file order. The vector is empty
    /// if `format` is unknown or `data` holds no records; reading stops at the
    /// first record that fails to parse.
    pub fn parse_many(data: &str, format: &str) -> Vec<Molecule> {
        with_ob(|| {
            ffi::mol_read_many(format, data)
                .iter()
                .map(|m| Molecule {
                    inner: ffi::mol_clone(m),
                })
                .collect()
        })
    }

    /// Read the first molecule from a file on disk.
    ///
    /// `format` is the OpenBabel format id (e.g. `"sdf"`); pass `None` to let
    /// OpenBabel infer it from the file's extension. Returns [`Error::Io`] if
    /// the file cannot be opened or its first record cannot be parsed.
    pub fn read_file(path: &str, format: Option<&str>) -> Result<Molecule, Error> {
        let inner = with_ob(|| ffi::mol_read_file(path, format.unwrap_or("")));
        if inner.is_null() {
            Err(Error::Io {
                path: path.to_string(),
            })
        } else {
            Ok(Molecule::from_inner(inner))
        }
    }

    /// Read **every** molecule from a (possibly multi-record) file on disk.
    ///
    /// `format` is the OpenBabel format id; pass `None` to infer it from the
    /// file's extension. Returns one [`Molecule`] per record, in file order.
    /// Returns [`Error::Io`] only if the file itself cannot be opened — a file
    /// that opens but yields no records gives an empty vector.
    pub fn read_file_many(path: &str, format: Option<&str>) -> Result<Vec<Molecule>, Error> {
        // The shim can't tell "missing file" from "no records"; check openability
        // here so a genuine I/O failure surfaces as an error.
        std::fs::File::open(path).map_err(|_| Error::Io {
            path: path.to_string(),
        })?;
        Ok(with_ob(|| {
            ffi::mol_read_file_many(path, format.unwrap_or(""))
                .iter()
                .map(|m| Molecule {
                    inner: ffi::mol_clone(m),
                })
                .collect()
        }))
    }

    /// Write this molecule to a file on disk.
    ///
    /// `format` is the OpenBabel format id; pass `None` to infer it from the
    /// file's extension. Returns [`Error::Io`] if the format cannot be resolved
    /// or the file cannot be written.
    pub fn write_file(&self, path: &str, format: Option<&str>) -> Result<(), Error> {
        let mut ok = false;
        with_ob(|| ffi::mol_write_file(self.as_inner(), path, format.unwrap_or(""), &mut ok));
        if ok {
            Ok(())
        } else {
            Err(Error::Io {
                path: path.to_string(),
            })
        }
    }

    /// Serialize this molecule to `format`, returning the text.
    ///
    /// Returns [`Error::UnknownFormat`] if OpenBabel doesn't know `format`.
    pub fn write(&self, format: &str) -> Result<String, Error> {
        let mut ok = true;
        let out = with_ob(|| ffi::mol_write(self.as_inner(), format, &mut ok));
        if ok {
            Ok(out)
        } else {
            Err(Error::UnknownFormat {
                format: format.to_string(),
            })
        }
    }

    /// Standard InChI identifier, e.g. `"InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"`.
    ///
    /// Returns `None` if InChI support is not compiled into the linked
    /// OpenBabel. Convenience wrapper over `write("inchi")`.
    pub fn inchi(&self) -> Option<String> {
        self.write("inchi").ok().map(|s| s.trim().to_string())
    }

    /// Standard InChIKey, e.g. `"LFQSCWFLJHTTHZ-UHFFFAOYSA-N"`.
    ///
    /// Returns `None` if InChI support is not compiled in. Convenience wrapper
    /// over `write("inchikey")`.
    pub fn inchikey(&self) -> Option<String> {
        self.write("inchikey").ok().map(|s| s.trim().to_string())
    }

    /// Molecular formula in Hill order, e.g. `"C2H6O"` (counts implicit H).
    pub fn formula(&self) -> String {
        with_ob(|| ffi::mol_formula(self.as_inner()))
    }

    /// Molecular formula with element counts spaced out, e.g. `"C 2 H 6 O 1"`.
    pub fn spaced_formula(&self) -> String {
        with_ob(|| ffi::mol_spaced_formula(self.as_inner()))
    }

    /// Standard molar mass in g/mol (counts implicit H).
    pub fn molar_mass(&self) -> f64 {
        with_ob(|| ffi::mol_mol_wt(self.as_inner()))
    }

    /// Monoisotopic exact mass (counts implicit H).
    pub fn exact_mass(&self) -> f64 {
        with_ob(|| ffi::mol_exact_mass(self.as_inner()))
    }

    /// Net formal charge of the molecule.
    pub fn total_charge(&self) -> i32 {
        with_ob(|| ffi::mol_total_charge(self.as_inner()))
    }

    /// Number of (explicit) atoms.
    pub fn num_atoms(&self) -> u32 {
        with_ob(|| ffi::mol_num_atoms(self.as_inner()))
    }

    /// Number of bonds.
    pub fn num_bonds(&self) -> u32 {
        with_ob(|| ffi::mol_num_bonds(self.as_inner()))
    }

    /// Number of heavy (non-hydrogen) atoms.
    pub fn num_heavy_atoms(&self) -> u32 {
        with_ob(|| ffi::mol_num_heavy_atoms(self.as_inner()))
    }

    /// Number of rings in the Smallest Set of Smallest Rings (SSSR).
    pub fn num_rings(&self) -> u32 {
        with_ob(|| ffi::mol_num_rings(self.as_inner()))
    }

    /// Number of rotatable bonds.
    pub fn num_rotatable_bonds(&self) -> u32 {
        with_ob(|| ffi::mol_num_rotors(self.as_inner()))
    }

    /// Total spin multiplicity of the molecule.
    pub fn spin_multiplicity(&self) -> u32 {
        with_ob(|| ffi::mol_spin_multiplicity(self.as_inner()))
    }

    /// The molecule's title (often a name or identifier; may be empty).
    pub fn title(&self) -> String {
        with_ob(|| ffi::mol_title(self.as_inner()))
    }

    /// Set the molecule's title.
    pub fn set_title(&mut self, title: &str) {
        with_ob(|| ffi::mol_set_title(self.inner.pin_mut(), title));
    }

    /// Make implicit hydrogens explicit (adds H atoms to the graph).
    pub fn add_hydrogens(&mut self) {
        with_ob(|| ffi::mol_add_hydrogens(self.inner.pin_mut()));
    }

    /// Remove explicit hydrogens (they become implicit again).
    pub fn remove_hydrogens(&mut self) {
        with_ob(|| ffi::mol_delete_hydrogens(self.inner.pin_mut()));
    }

    /// Evaluate a numeric descriptor plugin by id (e.g. `"logP"`, `"TPSA"`,
    /// `"MR"`, `"MW"`). Returns `None` if OpenBabel has no such descriptor.
    pub fn descriptor(&self, id: &str) -> Option<f64> {
        let mut ok = true;
        let value = with_ob(|| ffi::descriptor(self.as_inner(), id, &mut ok));
        if ok {
            Some(value)
        } else {
            None
        }
    }

    /// Predicted octanol/water partition coefficient (logP).
    pub fn logp(&self) -> Option<f64> {
        self.descriptor("logP")
    }

    /// Topological polar surface area (TPSA).
    pub fn tpsa(&self) -> Option<f64> {
        self.descriptor("TPSA")
    }

    /// Molar refractivity (MR).
    pub fn molar_refractivity(&self) -> Option<f64> {
        self.descriptor("MR")
    }

    /// Assign partial atomic charges using the named charge model
    /// (`"gasteiger"`, `"mmff94"`, `"eem"`, `"eqeq"`, `"qeq"`, `"qtpie"`).
    ///
    /// Returns `false` if the model is unknown or fails. Afterwards
    /// [`Atom::partial_charge`](crate::Atom::partial_charge) reflects the result.
    pub fn compute_charges(&mut self, model: &str) -> bool {
        with_ob(|| ffi::mol_compute_charges(self.inner.pin_mut(), model))
    }

    /// Least-squares superpose this molecule onto `reference` (Kabsch
    /// algorithm), updating this molecule's coordinates to the best-fit pose
    /// and returning the heavy-atom RMSD.
    ///
    /// Hydrogens are excluded from the fit and symmetry-equivalent atoms may be
    /// remapped for the best fit; use [`align_to_with`](Self::align_to_with) for
    /// control. Both molecules must have the same atoms in the same order (e.g.
    /// two conformers, or a structure and a transformed copy) and 3D
    /// coordinates. Returns `None` if the atom counts differ or alignment fails.
    pub fn align_to(&mut self, reference: &Molecule) -> Option<f64> {
        self.align_to_with(reference, false, true)
    }

    /// Like [`align_to`](Self::align_to), but with explicit control over whether
    /// hydrogens are included in the fit (`include_h`) and whether
    /// symmetry-equivalent atoms may be remapped (`symmetry`).
    pub fn align_to_with(
        &mut self,
        reference: &Molecule,
        include_h: bool,
        symmetry: bool,
    ) -> Option<f64> {
        let mut ok = true;
        let rmsd = with_ob(|| {
            ffi::mol_align(
                self.inner.pin_mut(),
                reference.as_inner(),
                include_h,
                symmetry,
                &mut ok,
            )
        });
        if ok {
            Some(rmsd)
        } else {
            None
        }
    }

    /// Coordinate dimension: `0` (no coordinates), `2`, or `3`.
    pub fn dimension(&self) -> u32 {
        with_ob(|| ffi::mol_dimension(self.as_inner()))
    }

    /// Whether the molecule has 3D coordinates.
    pub fn has_3d(&self) -> bool {
        self.dimension() == 3
    }

    /// Whether the molecule has 2D coordinates (e.g. for depiction).
    pub fn has_2d(&self) -> bool {
        self.dimension() == 2
    }

    /// Single-point energy of the molecule under the named force field
    /// (`"MMFF94"`, `"MMFF94s"`, `"UFF"`, `"GAFF"`, `"Ghemical"`).
    ///
    /// Returns `None` if the force field is unknown or cannot be set up.
    /// Only meaningful once the molecule has 3D coordinates (see
    /// [`generate_3d`](Self::generate_3d)). The unit is the force field's own
    /// (see [`forcefield_energy_unit`](crate::forcefield_energy_unit)).
    pub fn energy(&self, forcefield: &str) -> Option<f64> {
        let mut ok = true;
        let e = with_ob(|| ffi::mol_energy(self.as_inner(), forcefield, &mut ok));
        if ok { Some(e) } else { None }
    }

    /// Energy-minimize the geometry in place using `steps` conjugate-gradient
    /// steps of the named force field. Returns the final energy, or `None` if
    /// the force field is unknown or setup fails.
    ///
    /// For control over the algorithm, convergence, or constraints — or to
    /// stream the trajectory — use [`optimize_geometry_with`](Self::optimize_geometry_with)
    /// / [`minimize`](Self::minimize) with a [`Minimizer`](crate::Minimizer).
    pub fn optimize_geometry(&mut self, forcefield: &str, steps: u32) -> Option<f64> {
        let mut ok = true;
        let e = with_ob(|| ffi::mol_optimize(self.inner.pin_mut(), forcefield, steps, &mut ok));
        if ok { Some(e) } else { None }
    }

    /// Energy-minimize the geometry in place under a [`Minimizer`](crate::Minimizer)
    /// configuration — choosing the algorithm (steepest descent / conjugate
    /// gradients / L-BFGS), convergence threshold, and constraints — and return
    /// the final energy, or `None` if the force field is unknown or setup fails.
    ///
    /// Runs to completion without recording a trajectory; use
    /// [`minimize`](Self::minimize) to stream one instead. Only meaningful once
    /// the molecule has 3D coordinates (see [`generate_3d`](Self::generate_3d)).
    pub fn optimize_geometry_with(&mut self, config: &crate::Minimizer) -> Option<f64> {
        config.run(self)
    }

    /// Minimize the geometry and return its step-by-step trajectory.
    ///
    /// Runs the whole minimization under a [`Minimizer`](crate::Minimizer)
    /// configuration and returns an [`Optimization`](crate::Optimization) — an
    /// iterator over the recorded [`OptStep`](crate::OptStep) frames (step count,
    /// energy, coordinates), one per `steps_per_frame` steps until convergence or
    /// the step budget. The molecule is left holding the final geometry.
    ///
    /// The minimization runs eagerly when this is called (OpenBabel's shared
    /// force-field state means a run must be one atomic operation); the returned
    /// iterator replays the captured frames and does not borrow the molecule.
    ///
    /// ```no_run
    /// # use openbabel::{Minimizer, Molecule};
    /// # let mut mol = Molecule::parse("CCO", "smi").unwrap();
    /// # mol.generate_3d();
    /// let cfg = Minimizer::new("MMFF94");
    /// let trajectory: Vec<_> = mol.minimize(&cfg).collect();
    /// ```
    pub fn minimize(&mut self, config: &crate::Minimizer) -> crate::Optimization {
        crate::Optimization::new(self, config)
    }

    /// Generate 3D coordinates in place (like `obabel --gen3d`, "medium"
    /// quality: build from fragment templates, then force-field cleanup).
    /// Returns `false` if generation failed.
    pub fn generate_3d(&mut self) -> bool {
        with_ob(|| ffi::mol_make_3d(self.inner.pin_mut(), "med"))
    }

    /// Generate 3D coordinates with an explicit quality/speed setting, one of
    /// `"fastest"`, `"fast"`, `"med"`, `"slow"`, `"best"`.
    pub fn generate_3d_with(&mut self, speed: &str) -> bool {
        with_ob(|| ffi::mol_make_3d(self.inner.pin_mut(), speed))
    }

    /// Generate 2D coordinates in place (like `obabel --gen2d`), laying the
    /// molecule out for depiction. Returns `false` if generation failed.
    ///
    /// You don't need to call this before [`to_svg`](Self::to_svg): the SVG
    /// renderer generates 2D coordinates itself when they're absent.
    pub fn generate_2d(&mut self) -> bool {
        with_ob(|| ffi::mol_make_2d(self.inner.pin_mut()))
    }

    /// Render this molecule to an SVG document with default options.
    ///
    /// 2D coordinates are generated automatically if the molecule has none, so
    /// a freshly parsed molecule renders directly. Returns `None` on failure.
    pub fn to_svg(&self) -> Option<String> {
        self.to_svg_with(SvgOptions::default())
    }

    /// Render this molecule to an SVG document with explicit [`SvgOptions`].
    pub fn to_svg_with(&self, options: SvgOptions) -> Option<String> {
        let mut ok = true;
        let svg = with_ob(|| {
            ffi::mol_to_svg(
                self.as_inner(),
                options.all_carbons,
                options.atom_indices,
                &mut ok,
            )
        });
        if ok {
            Some(svg)
        } else {
            None
        }
    }

    /// Force (re)perception of stereochemistry from the molecule's structure
    /// (SMILES `@`/`@@` and `/\`, or 2D/3D coordinates).
    ///
    /// Stereo is perceived on demand by the query methods too; call this only
    /// to re-perceive after changing the structure.
    pub fn perceive_stereo(&mut self) {
        with_ob(|| ffi::mol_perceive_stereo(self.inner.pin_mut()));
    }

    /// Number of tetrahedral stereocenters perceived in the molecule.
    pub fn tetrahedral_stereo_count(&self) -> u32 {
        with_ob(|| ffi::mol_num_tetrahedral_stereo(self.as_inner()))
    }

    /// Number of cis/trans (double-bond) stereo units perceived.
    pub fn cistrans_stereo_count(&self) -> u32 {
        with_ob(|| ffi::mol_num_cistrans_stereo(self.as_inner()))
    }

    /// Run a genetic-algorithm conformer search targeting `count` diverse
    /// conformers, storing them in this molecule; returns the number now
    /// stored (see [`num_conformers`](Self::num_conformers)).
    ///
    /// The molecule must already have a 3D structure (call
    /// [`generate_3d`](Self::generate_3d) first). A rigid molecule with no
    /// rotatable bonds yields a single conformer.
    pub fn generate_conformers(&mut self, count: u32) -> u32 {
        with_ob(|| ffi::mol_generate_conformers(self.inner.pin_mut(), count))
    }

    /// Number of stored conformers (at least 1 once coordinates exist).
    pub fn num_conformers(&self) -> u32 {
        with_ob(|| ffi::mol_num_conformers(self.as_inner()))
    }

    /// Make conformer `index` the active coordinates. Out-of-range indices are
    /// ignored.
    ///
    /// To score every conformer, prefer [`conformer_energies`](Self::conformer_energies):
    /// it evaluates them all under one lock, whereas a `set_conformer` +
    /// [`energy`](Self::energy) loop leaves the force field's shared state open
    /// to interference from concurrent force-field calls between the two steps.
    pub fn set_conformer(&mut self, index: u32) {
        with_ob(|| ffi::mol_set_conformer(self.inner.pin_mut(), index));
    }

    /// The `(x, y, z)` coordinates of every atom in conformer `index`, in atom
    /// order — read without changing which conformer is active. `None` if
    /// `index` is out of range.
    pub fn conformer_coordinates(&self, index: u32) -> Option<Vec<[f64; 3]>> {
        if index >= self.num_conformers() {
            return None;
        }
        let flat = with_ob(|| ffi::mol_conformer_coordinates(self.as_inner(), index));
        Some(flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
    }

    /// The energy of every stored conformer under `forcefield`, in conformer
    /// order — the parallel of scoring each with [`energy`](Self::energy) but
    /// without disturbing the active conformer.
    ///
    /// Empty if `forcefield` is unknown (see
    /// [`energy`](Self::energy) for the valid ids) or no conformers exist.
    pub fn conformer_energies(&self, forcefield: &str) -> Vec<f64> {
        with_ob(|| ffi::mol_conformer_energies(self.as_inner(), forcefield))
    }

    /// Translate the molecule so its centroid sits at the origin.
    pub fn center(&mut self) {
        with_ob(|| ffi::mol_center(self.inner.pin_mut()));
    }

    /// Valence angle (in degrees) at atom `j` between atoms `i`, `j`, `k`
    /// (0-based indices). Requires coordinates; returns `0.0` for invalid
    /// indices.
    pub fn angle(&self, i: u32, j: u32, k: u32) -> f64 {
        with_ob(|| ffi::mol_angle(self.as_inner(), i + 1, j + 1, k + 1))
    }

    /// Torsion (dihedral) angle in degrees for atoms `i`, `j`, `k`, `l`
    /// (0-based indices). Requires coordinates; returns `0.0` for invalid
    /// indices.
    pub fn torsion(&self, i: u32, j: u32, k: u32, l: u32) -> f64 {
        with_ob(|| ffi::mol_torsion(self.as_inner(), i + 1, j + 1, k + 1, l + 1))
    }

    /// Remove disconnected fragments (e.g. counterions) smaller than
    /// `min_atoms` atoms; `0` keeps only the single largest fragment. Returns
    /// `true` if anything was removed.
    pub fn strip_salts(&mut self, min_atoms: u32) -> bool {
        with_ob(|| ffi::mol_strip_salts(self.inner.pin_mut(), min_atoms))
    }

    /// Split into its disconnected fragments, each returned as an owned
    /// molecule. A connected molecule yields a single element.
    pub fn separate(&self) -> Vec<Molecule> {
        with_ob(|| {
            ffi::mol_separate(self.as_inner())
                .iter()
                .map(|frag| Molecule {
                    inner: ffi::mol_clone(frag),
                })
                .collect()
        })
    }

    /// Read a string property previously attached under `key`, or `None`.
    pub fn property(&self, key: &str) -> Option<String> {
        let mut ok = true;
        let value = with_ob(|| ffi::mol_get_property(self.as_inner(), key, &mut ok));
        if ok {
            Some(value)
        } else {
            None
        }
    }

    /// Attach (or replace) a string property under `key`.
    pub fn set_property(&mut self, key: &str, value: &str) {
        with_ob(|| ffi::mol_set_property(self.inner.pin_mut(), key, value));
    }

    /// The atom at 0-based `index`, or `None` if out of range.
    pub fn atom(&self, index: u32) -> Option<Atom<'_>> {
        if index < self.num_atoms() {
            Some(Atom {
                mol: self.as_inner(),
                ob_idx: index + 1, // OpenBabel atoms are 1-based.
            })
        } else {
            None
        }
    }

    /// Iterate over the atoms in index order.
    pub fn atoms(&self) -> impl Iterator<Item = Atom<'_>> + '_ {
        let mol = self.as_inner();
        (0..self.num_atoms()).map(move |i| Atom { mol, ob_idx: i + 1 })
    }

    /// The bond at 0-based `index`, or `None` if out of range.
    pub fn bond(&self, index: u32) -> Option<Bond<'_>> {
        if index < self.num_bonds() {
            Some(Bond {
                mol: self.as_inner(),
                ob_idx: index, // OpenBabel bonds are 0-based.
            })
        } else {
            None
        }
    }

    /// Iterate over the bonds in index order.
    pub fn bonds(&self) -> impl Iterator<Item = Bond<'_>> + '_ {
        let mol = self.as_inner();
        (0..self.num_bonds()).map(move |i| Bond { mol, ob_idx: i })
    }

    /// The bond joining the atoms at 0-based indices `a` and `b`, or `None` if
    /// they are not directly bonded.
    pub fn bond_between(&self, a: u32, b: u32) -> Option<Bond<'_>> {
        let idx = with_ob(|| ffi::mol_bond_between(self.as_inner(), a, b));
        if idx < 0 {
            None
        } else {
            Some(Bond {
                mol: self.as_inner(),
                ob_idx: idx as u32,
            })
        }
    }

    /// Topological symmetry class of every atom, in atom order (parallel to
    /// [`atoms`](Self::atoms); entry `i` is atom `i`).
    ///
    /// Atoms that share a value are topologically equivalent (related by a graph
    /// automorphism) — e.g. the two methyl carbons of propane, or all six
    /// carbons of benzene. Computed with OpenBabel's `OBGraphSym`.
    pub fn symmetry_classes(&self) -> Vec<u32> {
        with_ob(|| ffi::mol_symmetry_classes(self.as_inner()))
    }

    /// A canonical rank (1-based) for every atom, in atom order (parallel to
    /// [`atoms`](Self::atoms); entry `i` is atom `i`).
    ///
    /// The ranks are a repeatable canonical labelling of the graph — the same
    /// molecule yields the same ranks regardless of input atom order — built
    /// from the [`symmetry_classes`](Self::symmetry_classes) via OpenBabel's
    /// canonical-labelling algorithm. Every atom gets a distinct rank.
    pub fn canonical_ranks(&self) -> Vec<u32> {
        with_ob(|| ffi::mol_canonical_ranks(self.as_inner()))
    }

    /// Find every unique way `query` occurs as a substructure of this molecule.
    ///
    /// Each returned mapping lists, for every atom of `query` (in `query`'s atom
    /// order), the 0-based index of the matching atom in `self`. Matching is by
    /// element and bond order (an exact subgraph isomorphism, `OBQuery` +
    /// VF2), not a flexible SMARTS query — use [`SmartsPattern`](crate::SmartsPattern)
    /// for that. The result is empty when `query` does not occur.
    ///
    /// "Unique" means no two mappings cover exactly the same set of `self`
    /// atoms (symmetry-equivalent hits are collapsed).
    pub fn substructure_search(&self, query: &Molecule) -> Vec<Vec<u32>> {
        let mut width = 0u32;
        let flat = with_ob(|| ffi::mol_substructure_mappings(query.as_inner(), self.as_inner(), &mut width));
        chunk(flat, width as usize)
    }

    /// Whether `query` occurs as a substructure of this molecule (see
    /// [`substructure_search`](Self::substructure_search)).
    pub fn has_substructure(&self, query: &Molecule) -> bool {
        !self.substructure_search(query).is_empty()
    }

    /// Every graph automorphism of this molecule.
    ///
    /// Each automorphism is a permutation of the atoms (a `Vec` indexed by
    /// 0-based atom index, whose entries are the 0-based indices the atoms map
    /// to) that preserves the molecular graph. The number of automorphisms is
    /// the order of the molecule's symmetry group — e.g. 12 for benzene. Built
    /// on the same symmetry perception as [`symmetry_classes`](Self::symmetry_classes).
    pub fn automorphisms(&self) -> Vec<Vec<u32>> {
        let mut width = 0u32;
        let flat = with_ob(|| ffi::mol_automorphisms(self.as_inner(), &mut width));
        chunk(flat, width as usize)
    }

    /// Set the torsion angle defined by the four atoms at 0-based indices `a`,
    /// `b`, `c`, `d` to `radians`, rotating the atoms on the far side of the
    /// `b`–`c` bond. Requires 3D coordinates; out-of-range indices are ignored.
    ///
    /// The angle is in **radians** (OpenBabel's `SetTorsion` convention), unlike
    /// [`torsion`](Self::torsion), which reports degrees.
    pub fn set_torsion(&mut self, a: u32, b: u32, c: u32, d: u32, radians: f64) {
        with_ob(|| ffi::mol_set_torsion(self.inner.pin_mut(), a, b, c, d, radians));
    }

    /// The 0-based indices of the atoms reachable from atom `to` without passing
    /// back through atom `from` — i.e. the fragment on `to`'s side of the
    /// `from`–`to` bond (excluding both endpoints). This is the set of atoms
    /// that move when [`set_torsion`](Self::set_torsion) rotates that bond.
    pub fn find_children(&self, from: u32, to: u32) -> Vec<u32> {
        with_ob(|| ffi::mol_find_children(self.as_inner(), from, to))
    }

    /// The 0-based indices of the atoms in the largest connected fragment.
    ///
    /// Useful for isolating the main molecule from counter-ions or solvent
    /// without modifying it (unlike [`strip_salts`](Self::strip_salts)).
    pub fn largest_fragment_atoms(&self) -> Vec<u32> {
        with_ob(|| ffi::mol_largest_fragment(self.as_inner()))
    }

    /// Set the molecule's total formal charge (used by some formats and
    /// perception steps). See [`total_charge`](Self::total_charge) to read it.
    pub fn set_total_charge(&mut self, charge: i32) {
        with_ob(|| ffi::mol_set_total_charge(self.inner.pin_mut(), charge));
    }

    /// Set the molecule's total spin multiplicity (2S+1). See
    /// [`spin_multiplicity`](Self::spin_multiplicity) to read it.
    pub fn set_total_spin_multiplicity(&mut self, spin: u32) {
        with_ob(|| ffi::mol_set_total_spin(self.inner.pin_mut(), spin));
    }

    /// Whether this molecule carries crystallographic unit-cell information
    /// (present after reading a crystal format such as CIF).
    pub fn has_unit_cell(&self) -> bool {
        with_ob(|| ffi::mol_has_unit_cell(self.as_inner()))
    }

    /// The crystallographic [`UnitCell`](crate::UnitCell), or `None` if this
    /// molecule has none.
    pub fn unit_cell(&self) -> Option<crate::UnitCell<'_>> {
        if self.has_unit_cell() {
            Some(crate::unitcell::UnitCell::new(self.as_inner()))
        } else {
            None
        }
    }

    /// Number of residues. `0` for molecules without biopolymer/PDB structure
    /// (e.g. anything parsed from SMILES).
    pub fn num_residues(&self) -> u32 {
        with_ob(|| ffi::mol_num_residues(self.as_inner()))
    }

    /// The residue at 0-based `index`, or `None` if out of range.
    pub fn residue(&self, index: u32) -> Option<Residue<'_>> {
        if index < self.num_residues() {
            Some(Residue::new(self.as_inner(), index))
        } else {
            None
        }
    }

    /// Iterate over the residues in index order.
    pub fn residues(&self) -> impl Iterator<Item = Residue<'_>> + '_ {
        let mol = self.as_inner();
        (0..self.num_residues()).map(move |i| Residue::new(mol, i))
    }

    /// The ring at 0-based `index` (into the SSSR — smallest set of smallest
    /// rings), or `None` if out of range. See [`num_rings`](Self::num_rings).
    pub fn ring(&self, index: u32) -> Option<Ring<'_>> {
        if index < self.num_rings() {
            Some(Ring::new(self.as_inner(), index))
        } else {
            None
        }
    }

    /// Iterate over the rings of the SSSR.
    pub fn rings(&self) -> impl Iterator<Item = Ring<'_>> + '_ {
        let mol = self.as_inner();
        (0..self.num_rings()).map(move |i| Ring::new(mol, i))
    }

    // --- Construction & editing ------------------------------------------

    /// Add an atom of atomic number `atomic_number` and return its 0-based
    /// index. Wrap a batch of edits in [`begin_modify`](Self::begin_modify) /
    /// [`end_modify`](Self::end_modify) to defer perception until the end.
    pub fn add_atom(&mut self, atomic_number: u32) -> u32 {
        with_ob(|| ffi::mol_add_atom(self.inner.pin_mut(), atomic_number))
    }

    /// Bond the atoms at 0-based indices `begin` and `end` with `order`
    /// (1 = single, 2 = double, 3 = triple). Returns `false` if either index is
    /// out of range.
    pub fn add_bond(&mut self, begin: u32, end: u32, order: u32) -> bool {
        with_ob(|| ffi::mol_add_bond(self.inner.pin_mut(), begin, end, order))
    }

    /// Delete the atom at 0-based `index` (and its bonds). Returns `false` if
    /// out of range. Note that deletion renumbers later atoms.
    pub fn delete_atom(&mut self, index: u32) -> bool {
        with_ob(|| ffi::mol_delete_atom(self.inner.pin_mut(), index))
    }

    /// Delete the bond at 0-based `index`. Returns `false` if out of range.
    pub fn delete_bond(&mut self, index: u32) -> bool {
        with_ob(|| ffi::mol_delete_bond(self.inner.pin_mut(), index))
    }

    /// Suspend perception (aromaticity, rings, …) while a batch of structural
    /// edits is applied; pair with [`end_modify`](Self::end_modify). Calls
    /// nest. Building a molecule inside a modify block is much faster than
    /// re-perceiving after every edit.
    pub fn begin_modify(&mut self) {
        with_ob(|| ffi::mol_begin_modify(self.inner.pin_mut()));
    }

    /// Resume perception after [`begin_modify`](Self::begin_modify).
    pub fn end_modify(&mut self) {
        with_ob(|| ffi::mol_end_modify(self.inner.pin_mut()));
    }

    /// Remove every atom and bond, leaving an empty molecule.
    pub fn clear(&mut self) {
        with_ob(|| ffi::mol_clear(self.inner.pin_mut()));
    }

    /// Translate every atom by `(dx, dy, dz)`.
    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        with_ob(|| ffi::mol_translate(self.inner.pin_mut(), dx, dy, dz));
    }

    /// Overwrite all atom coordinates from a flat `[x0, y0, z0, x1, …]` slice.
    /// Returns `false` unless `coords.len() == 3 * num_atoms`.
    pub fn set_coordinates(&mut self, coords: &[f64]) -> bool {
        with_ob(|| ffi::mol_set_coordinates(self.inner.pin_mut(), coords))
    }

    /// Set the coordinate dimension (0, 2, or 3). Mark a hand-built structure
    /// as `3` before calling [`connect_the_dots`](Self::connect_the_dots),
    /// which only runs on 3D structures.
    pub fn set_dimension(&mut self, dimension: u32) {
        with_ob(|| ffi::mol_set_dimension(self.inner.pin_mut(), dimension));
    }

    /// Infer connectivity (bonds) from 3D coordinates using covalent radii —
    /// the counterpart to reading a coordinates-only format. Requires the
    /// dimension to be 3 (see [`set_dimension`](Self::set_dimension)); follow
    /// with [`perceive_bond_orders`](Self::perceive_bond_orders) to assign
    /// orders.
    pub fn connect_the_dots(&mut self) {
        with_ob(|| ffi::mol_connect_the_dots(self.inner.pin_mut()));
    }

    /// Assign bond orders (single / double / triple, aromaticity) from the 3D
    /// geometry and connectivity.
    pub fn perceive_bond_orders(&mut self) {
        with_ob(|| ffi::mol_perceive_bond_orders(self.inner.pin_mut()));
    }

    /// Add only polar hydrogens (those on N, O, P, S). Returns `false` on
    /// failure.
    pub fn add_polar_hydrogens(&mut self) -> bool {
        with_ob(|| ffi::mol_add_polar_hydrogens(self.inner.pin_mut()))
    }

    /// Add hydrogens with pH-based (de)protonation: acidic/basic groups gain or
    /// lose H as appropriate for `ph`. Returns `false` on failure.
    pub fn add_hydrogens_for_ph(&mut self, ph: f64) -> bool {
        with_ob(|| ffi::mol_add_hydrogens_ph(self.inner.pin_mut(), ph))
    }

    /// Convert dative/coordinate bonds to their charge-separated form (e.g. a
    /// neutral nitro group to `-[N+](=O)[O-]`). Returns `false` if nothing
    /// changed.
    pub fn convert_dative_bonds(&mut self) -> bool {
        with_ob(|| ffi::mol_convert_dative_bonds(self.inner.pin_mut()))
    }

    /// (Re)assign radical spin multiplicities from the atomic valences. Returns
    /// `false` on failure.
    pub fn assign_spin_multiplicity(&mut self) -> bool {
        with_ob(|| ffi::mol_assign_spin_multiplicity(self.inner.pin_mut()))
    }

    /// A mutable handle to the atom at 0-based `index`, for setting its
    /// properties, or `None` if out of range.
    pub fn atom_mut(&mut self, index: u32) -> Option<AtomMut<'_>> {
        if index < self.num_atoms() {
            Some(AtomMut::new(self, index + 1)) // OpenBabel atoms are 1-based.
        } else {
            None
        }
    }

    /// A mutable handle to the bond at 0-based `index`, or `None` if out of
    /// range.
    pub fn bond_mut(&mut self, index: u32) -> Option<BondMut<'_>> {
        if index < self.num_bonds() {
            Some(BondMut::new(self, index)) // OpenBabel bonds are 0-based.
        } else {
            None
        }
    }

    /// Compute the Spectrophore™ descriptor.
    ///
    /// Returns the descriptor vector (48 values with the default settings). The
    /// molecule needs a 3D conformer — call [`generate_3d`](Self::generate_3d)
    /// first for structures without coordinates. Returns an empty vector on
    /// failure.
    pub fn spectrophore(&self) -> Vec<f64> {
        with_ob(|| ffi::mol_spectrophore(self.as_inner()))
    }

    /// Vibrational frequencies (cm⁻¹).
    ///
    /// Only populated when the molecule was read from a computational-chemistry
    /// output that carries vibration data (e.g. a Gaussian/ORCA log); otherwise
    /// empty.
    pub fn vibration_frequencies(&self) -> Vec<f64> {
        with_ob(|| ffi::mol_vibration_frequencies(self.as_inner()))
    }

    /// Vibrational IR intensities (km/mol), paired index-wise with
    /// [`vibration_frequencies`](Self::vibration_frequencies). Empty unless the
    /// molecule carries vibration data.
    pub fn vibration_intensities(&self) -> Vec<f64> {
        with_ob(|| ffi::mol_vibration_intensities(self.as_inner()))
    }

    /// Wrap a non-null `ffi::Molecule` owner. Callers must have checked the
    /// pointer is non-null (e.g. an out-of-range accessor returns null).
    pub(crate) fn from_inner(inner: UniquePtr<ffi::Molecule>) -> Self {
        debug_assert!(!inner.is_null(), "from_inner requires a non-null molecule");
        Molecule { inner }
    }

    pub(crate) fn as_inner(&self) -> &ffi::Molecule {
        self.inner.as_ref().expect("Molecule is never null")
    }

    pub(crate) fn as_inner_pin_mut(&mut self) -> std::pin::Pin<&mut ffi::Molecule> {
        self.inner.pin_mut()
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for Molecule {
    /// Deep-copy the molecule (a fully independent `OBMol`).
    fn clone(&self) -> Self {
        Molecule {
            inner: with_ob(|| ffi::mol_clone(self.as_inner())),
        }
    }
}

impl std::fmt::Debug for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Molecule")
            .field("formula", &self.formula())
            .field("num_atoms", &self.num_atoms())
            .field("num_bonds", &self.num_bonds())
            .finish()
    }
}
