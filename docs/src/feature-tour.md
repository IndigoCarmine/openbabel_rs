# Feature Tour

A category-by-category map of what the `openbabel` crate can do. For the exact
signatures of every method mentioned here, see the
[API Reference](./api-reference.md).

## Core: molecules, atoms, bonds

- **I/O** — `Molecule::parse(data, format)` and `Molecule::write(format)` read
  and write any OpenBabel format (`smi`, `mol`, `sdf`, `pdb`, `inchi`, canonical
  SMILES `can`, …).
- **Properties** — `formula`, `molar_mass`, `exact_mass`, `total_charge`,
  `num_atoms`, `num_bonds`, `title`.
- **Atoms** (`Molecule::atoms` / `atom`) — element, coordinates, formal charge,
  aromaticity, ring membership, and much more (see *Extended queries* below).
- **Bonds** (`Molecule::bonds` / `bond`) — begin/end atoms, order, aromaticity,
  ring membership.
- **Hydrogens** — `add_hydrogens` / `remove_hydrogens`.

## Analysis

- **SMARTS substructure matching** — `SmartsPattern::new(pattern)`, then
  `matches`, `num_matches`, `match_indices`.
- **Fingerprints & similarity** — `Fingerprint::compute(&mol, id)` for `FP2`,
  `FP3`, `FP4`, `MACCS`, …, and `Fingerprint::tanimoto` for similarity.
- **Numeric descriptors** — `Molecule::descriptor(id)`, plus the convenience
  methods `logp`, `tpsa`, `molar_refractivity`.

## 3D structures & force fields

- **3D generation** — `generate_3d` / `generate_3d_with`; query with `dimension`
  and `has_3d`.
- **Force fields** — single-point `energy(forcefield)` and geometry
  optimization `optimize_geometry(forcefield, steps)` with `MMFF94`, `MMFF94s`,
  `UFF`, `GAFF`, `Ghemical`. The energy unit for a force field is available via
  the free function `forcefield_energy_unit`.
- **Configurable / constrained minimization** — `Minimizer` +
  `optimize_geometry_with`, and the pull-based `minimize` trajectory iterator.
  See [Key Concepts](./concepts.md#constrained-minimization).

## Charges, atom chemistry & identifiers

- **Partial charges** — `compute_charges(model)` with `gasteiger`, `mmff94`,
  `eem`, `eqeq`, `qeq`, `qtpie`; per-atom values via `Atom::partial_charge`.
- **Atom chemistry** — `degree`, `total_valence`, `implicit_hydrogens`,
  `hybridization`, `is_hbond_donor`, `is_hbond_acceptor`.
- **Identifiers** — `inchi` / `inchikey`, backed by the InChI library bundled
  with OpenBabel.

## Geometry & structure alignment

- **Alignment** — `align_to` / `align_to_with` least-squares-superpose one
  molecule onto another (Kabsch), returning the RMSD and updating coordinates in
  place. Backed by `OBAlign` (enabled by the vendored Eigen). See
  [Key Concepts](./concepts.md#structure-alignment).
- **Measurements** — `center`, `angle`, `torsion`.

## 2D depiction

- **SVG** — `to_svg` / `to_svg_with(SvgOptions)` render a molecule to an SVG
  image; 2D coordinates are laid out automatically. Also `generate_2d` /
  `has_2d`. The SVG painter is built in (no Cairo dependency).

## Stereochemistry

- Perceive and query stereocenters: `tetrahedral_stereo_count`,
  `cistrans_stereo_count`, `perceive_stereo`; `Atom::is_tetrahedral_stereo` /
  `stereo_winding` (a [`Winding`], clockwise/anticlockwise — OpenBabel's
  descriptor, not a CIP `R`/`S` label); `Bond::is_cistrans_stereo`.

## Reaction / SMIRKS-like transforms

- `Transform::new(reactant_smarts, product_smarts)` compiles a SMARTS→SMARTS
  edit; `Transform::apply(&mut mol)` rewrites every match in place. See
  [Key Concepts](./concepts.md#smarts-to-smarts-transforms).

## Conformer search

- `generate_conformers(count)` runs a genetic-algorithm conformer search (needs
  a 3D structure); walk results with `num_conformers` / `set_conformer` and
  score each with `energy`.

## Element data & extended queries

- **`elements` module** — periodic-table data: `symbol`, `name`,
  `atomic_number`, `mass`, `exact_mass`, `electronegativity`, `covalent_radius`,
  `vdw_radius`, `max_bonds`.
- **Richer atom queries** — `type_name`, `isotope`, `atomic_mass`, `exact_mass`,
  `spin_multiplicity`, `heavy_degree`, `hetero_degree`, `is_chiral`,
  `is_heteroatom`, `is_metal`, `is_polar_hydrogen`, and ring queries
  (`ring_count`, `smallest_ring_size`, `is_in_ring_size`).
- **Richer bond queries** — `length`, `equilibrium_length`, `is_rotor`,
  `is_amide`, `is_ester`, `is_carbonyl`, `is_closure`.
- **Molecule extras** — `num_heavy_atoms`, `num_rings`, `num_rotatable_bonds`,
  `spaced_formula`, `spin_multiplicity`, `clone`, `strip_salts`, `separate`
  (into fragments), and string `property` / `set_property` metadata.

## Residues (biopolymer / PDB substructure)

- `num_residues` / `residue` / `residues` expose the residue grouping that PDB
  and mmCIF carry. A [`Residue`] reports `name`, `number` / `number_string`,
  `chain`, `insertion_code`, atom counts, and its member `atoms`. From the atom
  side, `Atom::residue` links back, with `residue_atom_id`, `is_hetatm`, and
  `serial_number`. Molecules parsed from SMILES carry no residues — the API
  never synthesizes one.

## Spectra

- `spectrophore` computes the Spectrophore™ descriptor (48 values by default)
  from a 3D structure — a rotation-invariant shape/property fingerprint.
- `vibration_frequencies` / `vibration_intensities` read vibrational data when
  the molecule was parsed from a computational-chemistry output that carries it.
