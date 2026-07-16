//! Integration tests for the completed perception-flag, atom-predicate, and
//! bond-predicate surfaces (T24).

use openbabel::{Bond, Molecule};

#[test]
fn remaining_perception_flag_setters_round_trip() {
    let mut mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");

    macro_rules! round_trip {
        ($set:ident, $get:ident) => {
            mol.$set(true);
            assert!(mol.$get(), concat!(stringify!($get), " should be true"));
            mol.$set(false);
            assert!(!mol.$get(), concat!(stringify!($get), " should be false"));
        };
    }
    round_trip!(set_lssr_perceived, has_lssr_perceived);
    round_trip!(set_atom_types_perceived, has_atom_types_perceived);
    round_trip!(set_ring_types_perceived, has_ring_types_perceived);
    round_trip!(set_chirality_perceived, has_chirality_perceived);
    round_trip!(set_partial_charges_perceived, has_partial_charges_perceived);
    round_trip!(set_hybridization_perceived, has_hybridization_perceived);
    round_trip!(set_closure_bonds_perceived, has_closure_bonds_perceived);
    round_trip!(set_corrected_for_ph, is_corrected_for_ph);
    round_trip!(set_spin_multiplicity_assigned, has_spin_multiplicity_assigned);
}

#[test]
fn perception_flags_reflect_real_operations() {
    // Computing partial charges sets the cached flag.
    let mut mol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert!(!mol.has_partial_charges_perceived());
    assert!(mol.compute_charges("gasteiger"));
    assert!(mol.has_partial_charges_perceived());

    // Querying a hybridization perceives and caches it.
    let mol2 = Molecule::parse("CCO", "smi").expect("ethanol");
    assert!(!mol2.has_hybridization_perceived());
    let _ = mol2.atom(0).unwrap().hybridization();
    assert!(mol2.has_hybridization_perceived());
    // (The corrected-for-pH flag is exercised by the round-trip test above;
    // add_hydrogens_for_ph itself does not leave it set, because the internal
    // EndModify clears perception flags after CorrectForPH runs.)
}

#[test]
fn functional_group_atom_predicates() {
    // Acetic acid: both oxygens count as carboxyl oxygens; the carbonyl carbon
    // (index 1) has two "free" (single-heavy-valence) oxygens.
    let acid = Molecule::parse("CC(=O)O", "smi").expect("acetic acid");
    let carboxyl = (0..acid.num_atoms())
        .filter(|&i| acid.atom(i).unwrap().is_carboxyl_oxygen())
        .count();
    assert_eq!(carboxyl, 2);
    assert_eq!(acid.atom(1).unwrap().free_oxygen_count(), 2);

    // Nitromethane: two nitro oxygens.
    let nitro = Molecule::parse("C[N+](=O)[O-]", "smi").expect("nitromethane");
    let nitro_o = (0..nitro.num_atoms())
        .filter(|&i| nitro.atom(i).unwrap().is_nitro_oxygen())
        .count();
    assert_eq!(nitro_o, 2);

    // Acetamide: exactly one amide nitrogen.
    let amide = Molecule::parse("CC(=O)N", "smi").expect("acetamide");
    let amide_n = (0..amide.num_atoms())
        .filter(|&i| amide.atom(i).unwrap().is_amide_nitrogen())
        .count();
    assert_eq!(amide_n, 1);

    // Phosphate / sulfate oxygens are recognised.
    let phos = Molecule::parse("OP(=O)(O)O", "smi").expect("phosphate");
    assert!((0..phos.num_atoms()).any(|i| phos.atom(i).unwrap().is_phosphate_oxygen()));
    let sulf = Molecule::parse("OS(=O)(=O)O", "smi").expect("sulfate");
    assert!((0..sulf.num_atoms()).any(|i| sulf.atom(i).unwrap().is_sulfate_oxygen()));

    // Carbon disulfide: the central carbon has two free sulfurs.
    let cs2 = Molecule::parse("S=C=S", "smi").expect("carbon disulfide");
    let c = (0..cs2.num_atoms())
        .find(|&i| cs2.atom(i).unwrap().atomic_number() == 6)
        .unwrap();
    assert_eq!(cs2.atom(c).unwrap().free_sulfur_count(), 2);
}

#[test]
fn hydrogen_and_ring_atom_predicates() {
    // Ethanol with explicit H: the hydroxyl H is an H-bond donor H, and the
    // five C-H hydrogens are non-polar.
    let mut eth = Molecule::parse("CCO", "smi").expect("ethanol");
    eth.add_hydrogens();
    let donor_h = (0..eth.num_atoms())
        .filter(|&i| eth.atom(i).unwrap().is_hbond_donor_h())
        .count();
    assert_eq!(donor_h, 1);
    let nonpolar_h = (0..eth.num_atoms())
        .filter(|&i| eth.atom(i).unwrap().is_nonpolar_hydrogen())
        .count();
    assert_eq!(nonpolar_h, 5);

    // Benzene: every ring carbon has exactly two ring bonds.
    let benz = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    assert!((0..benz.num_atoms()).all(|i| benz.atom(i).unwrap().ring_bond_count() == 2));

    // Pyridine N-oxide: its ring nitrogen is an aromatic N-oxide. OpenBabel's
    // IsAromaticNOxide needs the exocyclic N=O as a real double bond, so convert
    // the charge-separated [n+]-[O-] SMILES form to the dative n=O first.
    let mut noxide = Molecule::parse("c1cc[n+]([O-])cc1", "smi").expect("pyridine N-oxide");
    assert!(noxide.convert_dative_bonds());
    assert!((0..noxide.num_atoms()).any(|i| noxide.atom(i).unwrap().is_aromatic_noxide()));
}

#[test]
fn bond_angles_from_explicit_geometry() {
    // Water built with a right-angle H-O-H (O at origin, H on +x and +y): the
    // one bond angle at O is exactly 90°.
    let mut w = Molecule::new();
    w.begin_modify();
    let o = w.add_atom(8);
    let h1 = w.add_atom(1);
    let h2 = w.add_atom(1);
    assert!(w.add_bond(o, h1, 1));
    assert!(w.add_bond(o, h2, 1));
    w.end_modify();
    w.set_dimension(3);
    #[rustfmt::skip]
    let coords = [0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0];
    assert!(w.set_coordinates(&coords));

    let o_atom = w.atom(0).unwrap();
    assert!((o_atom.smallest_bond_angle() - 90.0).abs() < 1e-6);
    assert!((o_atom.average_bond_angle() - 90.0).abs() < 1e-6);
}

#[test]
fn lewis_acid_base_counts_are_nonnegative() {
    let mut methane = Molecule::parse("C", "smi").expect("methane");
    methane.add_hydrogens();
    let (acid, base) = methane.atom(0).unwrap().lewis_acid_base_counts();
    // A saturated carbon has a full octet: no accepting or donating vacancies.
    assert!(acid >= 0 && base >= 0);
    assert_eq!((acid, base), (0, 0));
}

#[test]
fn amide_type_bond_predicates() {
    fn amide_bond(mol: &Molecule) -> Bond<'_> {
        (0..mol.num_bonds())
            .map(|i| mol.bond(i).unwrap())
            .find(Bond::is_amide)
            .expect("has an amide bond")
    }

    let primary = Molecule::parse("CC(=O)N", "smi").expect("acetamide");
    let b = amide_bond(&primary);
    assert!(b.is_primary_amide());
    assert!(!b.is_secondary_amide());
    assert!(!b.is_tertiary_amide());

    let secondary = Molecule::parse("CC(=O)NC", "smi").expect("N-methylacetamide");
    assert!(amide_bond(&secondary).is_secondary_amide());

    let tertiary = Molecule::parse("CC(=O)N(C)C", "smi").expect("N,N-dimethylacetamide");
    assert!(amide_bond(&tertiary).is_tertiary_amide());
}

#[test]
fn double_bond_geometry_predicates() {
    // is_double_bond_geometry is a 3D check: the C=C is planar (torsion near
    // 0/180°) → true; twisted out of plane (torsion near 90°) → false. Build a
    // C0-C1=C2-C3 skeleton and place C3 either in-plane or out of plane.
    fn butene(c3: [f64; 3]) -> Molecule {
        let mut m = Molecule::new();
        m.begin_modify();
        for _ in 0..4 {
            m.add_atom(6);
        }
        m.add_bond(0, 1, 1);
        m.add_bond(1, 2, 2); // the double bond, index 1
        m.add_bond(2, 3, 1);
        m.end_modify();
        m.set_dimension(3);
        #[rustfmt::skip]
        let coords = [
            -0.5, 0.87, 0.0, // C0
             0.0, 0.0,  0.0, // C1
             1.33, 0.0, 0.0, // C2
            c3[0], c3[1], c3[2], // C3
        ];
        assert!(m.set_coordinates(&coords));
        m
    }

    // Planar (C3 anti, in the z=0 plane): torsion ≈ 180°.
    let planar = butene([1.83, -0.87, 0.0]);
    assert!(planar.bond(1).unwrap().is_double_bond_geometry());

    // Twisted (C3 lifted out of plane): torsion ≈ 90°.
    let twisted = butene([1.83, 0.0, 0.87]);
    assert!(!twisted.bond(1).unwrap().is_double_bond_geometry());

    // A plain saturated bond carries neither 2D stereo flag.
    let ethane = Molecule::parse("CC", "smi").expect("ethane");
    let single = ethane.bond(0).unwrap();
    assert!(!single.is_wedge_or_hash());
    assert!(!single.is_cis_or_trans());
}
