//! Integration tests for spectra descriptors (T9): Spectrophore and vibration
//! data.

use openbabel::Molecule;

/// Same `00T` ligand as the residue tests: a 22-atom structure carrying 3D
/// coordinates, enough for a Spectrophore calculation.
const LIGAND_PDB: &str = r#"ATOM      1  N19 00T A   1      -2.126  -2.241  -0.764  1.00 10.00           N
ATOM      2  C10 00T A   1      -1.652  -1.628   0.484  1.00 10.00           C
ATOM      3  C6  00T A   1      -0.631  -0.566   0.166  1.00 10.00           C
ATOM      4  C5  00T A   1      -1.043   0.686  -0.251  1.00 10.00           C
ATOM      5  C32 00T A   1      -2.513   0.988  -0.389  1.00 10.00           C
ATOM      6  N33 00T A   1      -3.078   1.283   0.934  1.00 10.00           N
ATOM      7  C4  00T A   1      -0.107   1.661  -0.542  1.00 10.00           C
ATOM      8  C3  00T A   1       1.241   1.384  -0.416  1.00 10.00           C
ATOM      9  C2  00T A   1       1.654   0.131   0.002  1.00 10.00           C
ATOM     10 CL1  00T A   1       3.348  -0.217   0.160  1.00 10.00          CL
ATOM     11  C7  00T A   1       0.717  -0.842   0.298  1.00 10.00           C
ATOM     12 HN19 00T A   1      -1.355  -2.609  -1.301  1.00 10.00           H
ATOM     13  H10 00T A   1      -2.494  -1.177   1.010  1.00 10.00           H
ATOM     14 H10A 00T A   1      -1.197  -2.392   1.115  1.00 10.00           H
ATOM     15  H32 00T A   1      -2.647   1.850  -1.043  1.00 10.00           H
ATOM     16 H32A 00T A   1      -3.023   0.124  -0.818  1.00 10.00           H
ATOM     17 HN33 00T A   1      -2.580   2.039   1.379  1.00 10.00           H
ATOM     18 HN3A 00T A   1      -4.064   1.489   0.868  1.00 10.00           H
ATOM     19  H4  00T A   1      -0.429   2.638  -0.868  1.00 10.00           H
ATOM     20  H3  00T A   1       1.972   2.146  -0.644  1.00 10.00           H
ATOM     21  H7  00T A   1       1.039  -1.819   0.629  1.00 10.00           H
ATOM     22  H2  00T A   1      -2.812  -2.957  -0.577  1.00 10.00           H
END
"#;

#[test]
fn spectrophore_from_3d_structure() {
    let mol = Molecule::parse(LIGAND_PDB, "pdb").expect("parse PDB");
    let sp = mol.spectrophore();
    // The default Spectrophore is 48 real values (4 properties × 12 probes).
    assert_eq!(sp.len(), 48, "default Spectrophore length");
    assert!(sp.iter().all(|v| v.is_finite()), "all values finite");
    assert!(sp.iter().any(|&v| v != 0.0), "descriptor is not all-zero");
}

#[test]
fn spectrophore_needs_enough_atoms() {
    // OpenBabel refuses to compute a Spectrophore for fewer than 3 atoms and
    // returns an empty result.
    let mol = Molecule::parse("O", "smi").expect("parse SMILES");
    assert!(mol.spectrophore().is_empty());
}

#[test]
fn vibration_data_absent_without_compchem_input() {
    // A structure not read from a computational-chemistry output carries no
    // vibrational data.
    let mol = Molecule::parse(LIGAND_PDB, "pdb").expect("parse PDB");
    assert!(mol.vibration_frequencies().is_empty());
    assert!(mol.vibration_intensities().is_empty());
}
