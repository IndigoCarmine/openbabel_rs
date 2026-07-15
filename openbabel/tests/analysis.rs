//! Integration tests for the T2 "analysis" surface: SMARTS substructure
//! matching, fingerprints/similarity, and molecular descriptors.

use openbabel::{Fingerprint, Molecule, SmartsPattern};

#[test]
fn smarts_matches_hydroxyl_in_phenol() {
    let phenol = Molecule::parse("c1ccccc1O", "smi").expect("parse phenol");
    let hydroxyl = SmartsPattern::new("[OX2H]").expect("compile SMARTS");

    assert_eq!(hydroxyl.atom_count(), 1);
    assert!(hydroxyl.matches(&phenol), "phenol has a hydroxyl");
    assert_eq!(hydroxyl.num_matches(&phenol), 1);

    let matches = hydroxyl.match_indices(&phenol);
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].len(), 1);
    // The matched atom should be the oxygen (0-based index 6 here).
    let o_idx = matches[0][0];
    let o = phenol.atom(o_idx).expect("matched atom exists");
    assert_eq!(o.atomic_number(), 8, "hydroxyl match should be oxygen");
}

#[test]
fn smarts_counts_aromatic_carbons() {
    let toluene = Molecule::parse("Cc1ccccc1", "smi").expect("parse toluene");
    let aromatic_c = SmartsPattern::new("c").expect("compile SMARTS");
    // Six aromatic carbons in the ring; the methyl carbon is aliphatic.
    assert_eq!(aromatic_c.num_matches(&toluene), 6);
}

#[test]
fn invalid_smarts_is_an_error() {
    assert!(SmartsPattern::new("[[[bad").is_err());
}

#[test]
fn fingerprint_similarity_is_bounded_and_ordered() {
    let benzene = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let toluene = Molecule::parse("Cc1ccccc1", "smi").expect("toluene");
    let octane = Molecule::parse("CCCCCCCC", "smi").expect("octane");

    let fp_benzene = Fingerprint::compute(&benzene, "FP2").expect("fp benzene");
    let fp_toluene = Fingerprint::compute(&toluene, "FP2").expect("fp toluene");
    let fp_octane = Fingerprint::compute(&octane, "FP2").expect("fp octane");

    // Self-similarity is exactly 1.
    let self_sim = fp_benzene.tanimoto(&fp_benzene);
    assert!((self_sim - 1.0).abs() < 1e-9, "self similarity was {self_sim}");

    let sim_aromatic = fp_benzene.tanimoto(&fp_toluene);
    let sim_mixed = fp_benzene.tanimoto(&fp_octane);
    for s in [sim_aromatic, sim_mixed] {
        assert!((0.0..=1.0).contains(&s), "similarity {s} out of range");
    }
    // Two aromatics resemble each other more than an aromatic and an alkane.
    assert!(
        sim_aromatic > sim_mixed,
        "benzene~toluene ({sim_aromatic}) should exceed benzene~octane ({sim_mixed})"
    );
}

#[test]
fn descriptors_are_reasonable() {
    let ethanol = Molecule::parse("CCO", "smi").expect("ethanol");

    // TPSA of ethanol is one hydroxyl oxygen ~ 20.23 Å².
    let tpsa = ethanol.tpsa().expect("TPSA descriptor available");
    assert!((tpsa - 20.23).abs() < 1.0, "unexpected TPSA: {tpsa}");

    // logP of ethanol is small (slightly negative to near zero).
    let logp = ethanol.logp().expect("logP descriptor available");
    assert!(logp < 1.0, "ethanol logP should be low, got {logp}");

    // An unknown descriptor id yields None rather than panicking.
    assert!(ethanol.descriptor("definitely-not-a-descriptor").is_none());
}
