//! Port of OpenBabel's `test/ringtest.cpp` — ring perception (SSSR).
//!
//! Faithful reproduction: reads the same SMILES corpus (`attype.00.smi`) and
//! the same expected-results file (`ringresults.txt`) from OpenBabel's
//! `test/files/`, and for every molecule checks, exactly as the C++ test does:
//!   1. each bond's ring membership (`OBBond::IsInRing`),
//!   2. the SSSR ring count (`OBMol::GetSSSR().size()`),
//!   3. per-atom ring-membership counts (`OBAtom::MemberOfRingCount`).
//!
//! The C++ reference format (from `GenerateRingReference`) writes, per molecule,
//! three lines: the indices of in-ring bonds, the SSSR size, then one count per
//! atom.

mod common;

use common::{read_test_file, tokenize};
use openbabel::Molecule;

#[test]
fn ring_perception_matches_reference() {
    let smiles = read_test_file("attype.00.smi");
    let results = read_test_file("ringresults.txt");
    let mut ref_lines = results.lines();

    let mut molecules = 0usize;
    for line in smiles.lines() {
        // Mirror `conv.Read` + `if (mol.Empty()) continue;` — skip anything that
        // does not parse to a non-empty molecule, without consuming a reference
        // record.
        let mol = match Molecule::parse(line, "smi") {
            Ok(m) if m.num_atoms() > 0 => m,
            _ => continue,
        };
        molecules += 1;
        let title = mol.title();

        // --- (1) ring bonds -------------------------------------------------
        let bonds_line = ref_lines
            .next()
            .expect("ran out of reference data (ring bonds)");
        let ring_bonds: std::collections::HashSet<u32> = tokenize(bonds_line)
            .iter()
            .map(|t| t.parse::<u32>().expect("bond index"))
            .collect();
        for bond in mol.bonds() {
            let idx = bond.index();
            assert_eq!(
                ring_bonds.contains(&idx),
                bond.is_in_ring(),
                "ring bond data differs from reference for bond {idx} of {title:?}"
            );
        }

        // --- (2) SSSR size --------------------------------------------------
        let size_line = ref_lines
            .next()
            .expect("ran out of reference data (SSSR size)");
        let expected_sssr: u32 = size_line.trim().parse().expect("SSSR size");
        assert_eq!(
            mol.num_rings(),
            expected_sssr,
            "SSSR size differs from reference for {title:?}"
        );

        // --- (3) per-atom ring membership counts ----------------------------
        let counts_line = ref_lines
            .next()
            .expect("ran out of reference data (ring membership)");
        let counts = tokenize(counts_line);
        assert_eq!(
            counts.len() as u32,
            mol.num_atoms(),
            "reference has wrong number of atom ring-membership counts for {title:?}"
        );
        for (atom, expected) in mol.atoms().zip(counts.iter()) {
            let expected: u32 = expected.parse().expect("ring membership count");
            assert_eq!(
                atom.ring_count(),
                expected,
                "ring membership differs from reference for atom {} of {title:?}",
                atom.index()
            );
        }
    }

    assert!(molecules > 0, "no molecules were read from attype.00.smi");
}
