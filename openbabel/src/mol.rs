//! The [`Molecule`] type — a safe, owning wrapper around OpenBabel's `OBMol`.

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::atom::Atom;
use crate::bond::Bond;
use crate::error::Error;
use crate::residue::Residue;
use crate::with_ob;

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
    pub fn optimize_geometry(&mut self, forcefield: &str, steps: u32) -> Option<f64> {
        let mut ok = true;
        let e = with_ob(|| ffi::mol_optimize(self.inner.pin_mut(), forcefield, steps, &mut ok));
        if ok { Some(e) } else { None }
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
    /// ignored. Combine with [`energy`](Self::energy) to score each conformer.
    pub fn set_conformer(&mut self, index: u32) {
        with_ob(|| ffi::mol_set_conformer(self.inner.pin_mut(), index));
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
