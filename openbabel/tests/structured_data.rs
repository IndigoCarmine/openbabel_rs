//! Integration tests for structured torsion/angle enumeration, perception-flag
//! setters, and axial ring position (T23).

use openbabel::Molecule;

#[test]
fn find_angles_enumerates_valence_angles() {
    // Butane C0-C1-C2-C3 (no explicit H): only the two interior carbons have
    // two neighbours, so there are exactly two angles, centred on C1 and C2.
    let mol = Molecule::parse("CCCC", "smi").expect("butane");
    let mut angles = mol.find_angles();
    // Each triple is [vertex, a, b]; sort the endpoints for a stable compare.
    for a in &mut angles {
        let (v, mut ends) = (a[0], [a[1], a[2]]);
        ends.sort_unstable();
        *a = [v, ends[0], ends[1]];
    }
    angles.sort_unstable();
    assert_eq!(angles, vec![[1, 0, 2], [2, 1, 3]]);
}

#[test]
fn find_torsions_enumerates_dihedrals() {
    // Butane has a single heavy-atom torsion, around the central C1-C2 bond.
    let mol = Molecule::parse("CCCC", "smi").expect("butane");
    let torsions = mol.find_torsions();
    assert_eq!(torsions.len(), 1);

    let mut atoms = torsions[0];
    atoms.sort_unstable();
    assert_eq!(atoms, [0, 1, 2, 3]);
    // The central bond of the quad is C1-C2 (order may be forward or reverse).
    let t = torsions[0];
    assert!((t[1] == 1 && t[2] == 2) || (t[1] == 2 && t[2] == 1));

    // Ethanol has no heavy-atom torsion (no central bond with a neighbour on
    // both ends).
    let ethanol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert!(ethanol.find_torsions().is_empty());
}

#[test]
fn find_torsions_and_angles_stay_in_range_and_cache() {
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let n = mol.num_atoms();

    let angles = mol.find_angles();
    let torsions = mol.find_torsions();
    assert!(!angles.is_empty());
    assert!(!torsions.is_empty());
    for a in &angles {
        assert!(a.iter().all(|&i| i < n));
    }
    for t in &torsions {
        assert!(t.iter().all(|&i| i < n));
    }

    // OpenBabel caches the data on the molecule; a second call is identical.
    assert_eq!(mol.find_angles(), angles);
    assert_eq!(mol.find_torsions(), torsions);
}

#[test]
fn perception_flag_setters_round_trip() {
    let mut mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");

    // Force ring + aromatic perception, then flip the flags via the setters.
    let _ = mol.num_rings();
    let _ = mol.atom(0).unwrap().is_aromatic();
    assert!(mol.has_sssr_perceived());
    assert!(mol.has_aromatic_perceived());

    mol.set_sssr_perceived(false);
    mol.set_aromatic_perceived(false);
    assert!(!mol.has_sssr_perceived());
    assert!(!mol.has_aromatic_perceived());

    mol.set_sssr_perceived(true);
    mol.set_aromatic_perceived(true);
    assert!(mol.has_sssr_perceived());
    assert!(mol.has_aromatic_perceived());

    // The ring-atom and chain flags round-trip the same way.
    mol.set_ring_atoms_perceived(true);
    assert!(mol.has_ring_atoms_perceived());
    mol.set_ring_atoms_perceived(false);
    assert!(!mol.has_ring_atoms_perceived());

    mol.set_chains_perceived(true);
    assert!(mol.has_chains_perceived());
    mol.set_chains_perceived(false);
    assert!(!mol.has_chains_perceived());
}

#[test]
fn set_hydrogens_added_flag_only() {
    let mut mol = Molecule::parse("CCO", "smi").expect("ethanol");
    let n = mol.num_atoms();
    assert!(!mol.has_hydrogens_added());

    // Setting the flag is pure bookkeeping — it does not add atoms.
    mol.set_hydrogens_added(true);
    assert!(mol.has_hydrogens_added());
    assert_eq!(mol.num_atoms(), n);

    mol.set_hydrogens_added(false);
    assert!(!mol.has_hydrogens_added());
}

#[test]
fn is_axial_needs_geometry() {
    // Without coordinates, no atom can be axial.
    let mol = Molecule::parse("C1CCCCC1", "smi").expect("cyclohexane");
    assert!((0..mol.num_atoms()).all(|i| !mol.atom(i).unwrap().is_axial()));
}

#[test]
fn is_axial_finds_axial_substituents_in_a_chair() {
    // Build a full ideal-chair cyclohexane (C6H12) with explicit coordinates,
    // rather than relying on 3D generation — OpenBabel's fragment builder yields
    // NaN coordinates for bare cyclohexane. The ring is a true chair (C-C 1.54 Å,
    // C-C-C 111.4°, ring torsion ±55°), which puts each axial C-H bond at a ~65°
    // dihedral to the ring: inside OpenBabel's 55–75° axial window with margin,
    // while equatorial C-H bonds sit near 175°. Atoms 0..6 are the carbons;
    // atoms 6+2i and 7+2i are the two hydrogens on carbon i.
    #[rustfmt::skip]
    let coords: [[f64; 3]; 18] = [
        [ 1.4688,  0.0000,  0.2314], // C0
        [ 0.7344,  1.2720, -0.2314], // C1
        [-0.7344,  1.2720,  0.2314], // C2
        [-1.4688,  0.0000, -0.2314], // C3
        [-0.7344, -1.2720,  0.2314], // C4
        [ 0.7344, -1.2720, -0.2314], // C5
        [ 2.4757,  0.0000, -0.1860], // H on C0 (equatorial)
        [ 1.5267,  0.0000,  1.3199], // H on C0 (axial)
        [ 0.7633,  1.3222, -1.3199], // H on C1 (axial)
        [ 1.2379,  2.1440,  0.1860], // H on C1 (equatorial)
        [-1.2379,  2.1440, -0.1860], // H on C2 (equatorial)
        [-0.7633,  1.3222,  1.3199], // H on C2 (axial)
        [-1.5267,  0.0000, -1.3199], // H on C3 (axial)
        [-2.4757,  0.0000,  0.1860], // H on C3 (equatorial)
        [-1.2379, -2.1440, -0.1860], // H on C4 (equatorial)
        [-0.7633, -1.3222,  1.3199], // H on C4 (axial)
        [ 0.7633, -1.3222, -1.3199], // H on C5 (axial)
        [ 1.2379, -2.1440,  0.1860], // H on C5 (equatorial)
    ];

    let mut mol = Molecule::new();
    mol.begin_modify();
    for _ in 0..6 {
        mol.add_atom(6); // carbons
    }
    for _ in 0..12 {
        mol.add_atom(1); // hydrogens
    }
    for i in 0..6u32 {
        assert!(mol.add_bond(i, (i + 1) % 6, 1)); // ring bonds
        assert!(mol.add_bond(i, 6 + 2 * i, 1)); // two C-H bonds
        assert!(mol.add_bond(i, 7 + 2 * i, 1));
    }
    mol.end_modify();

    mol.set_dimension(3);
    let flat: Vec<f64> = coords.iter().flat_map(|c| *c).collect();
    assert!(mol.set_coordinates(&flat));

    // Ring carbons are never axial themselves (all their heavy bonds are ring
    // bonds); only the substituent hydrogens can be.
    assert!((0..6).all(|i| !mol.atom(i).unwrap().is_axial()));

    // A chair has exactly six axial and six equatorial C-H bonds.
    let axial = (0..mol.num_atoms())
        .filter(|&i| mol.atom(i).unwrap().is_axial())
        .count();
    assert_eq!(axial, 6, "a chair should have six axial hydrogens, got {axial}");
}
