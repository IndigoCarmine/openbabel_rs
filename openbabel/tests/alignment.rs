//! Integration tests for least-squares structure alignment (T5, `OBAlign`).
//!
//! These require Eigen to be compiled into OpenBabel (`EIGEN3_INCLUDE_DIR` in
//! `openbabel-sys/build.rs`); without it `OBAlign` is not built.

use openbabel::Molecule;

/// Build an XYZ-format string from `(element, x, y, z)` rows.
fn xyz(atoms: &[(&str, f64, f64, f64)]) -> String {
    let mut s = format!("{}\ngenerated\n", atoms.len());
    for (el, x, y, z) in atoms {
        s.push_str(&format!("{el} {x:.6} {y:.6} {z:.6}\n"));
    }
    s
}

/// Reference geometry: four distinct, non-coplanar atoms (a rank-3 point set,
/// so the 3D alignment is well determined). Distinct elements mean there are no
/// symmetry-equivalent atoms, so the atom mapping is unique.
fn reference_atoms() -> Vec<(&'static str, f64, f64, f64)> {
    vec![
        ("C", 0.0, 0.0, 0.0),
        ("O", 1.3, 0.0, 0.0),
        ("N", 0.0, 1.4, 0.0),
        ("F", 0.0, 0.0, 1.5),
    ]
}

#[test]
fn rigid_transform_aligns_to_zero_rmsd() {
    let reference = Molecule::parse(&xyz(&reference_atoms()), "xyz").expect("ref");

    // Rotate each atom 90° about z: (x, y, z) -> (-y, x, z); then translate.
    // A pure rotation + translation can be undone exactly, so RMSD -> ~0.
    let (tx, ty, tz) = (5.0, -3.0, 2.0);
    let moved: Vec<(&str, f64, f64, f64)> = reference_atoms()
        .iter()
        .map(|&(el, x, y, z)| (el, -y + tx, x + ty, z + tz))
        .collect();
    let mut target = Molecule::parse(&xyz(&moved), "xyz").expect("target");

    let rmsd = target.align_to(&reference).expect("alignment should succeed");
    assert!(rmsd < 1e-3, "rigid transform should align to ~0 RMSD, got {rmsd}");

    // After alignment the target overlays the reference atom-for-atom.
    for (t, r) in target.atoms().zip(reference.atoms()) {
        let (ax, ay, az) = t.coords();
        let (bx, by, bz) = r.coords();
        let d = ((ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2)).sqrt();
        assert!(d < 1e-2, "atom {} not overlaid onto reference: d={d}", t.index());
    }
}

#[test]
fn distinct_structures_have_positive_rmsd() {
    let reference = Molecule::parse(&xyz(&reference_atoms()), "xyz").expect("ref");

    // Move one atom substantially — no rigid transform can overlay this, so the
    // best-fit RMSD must stay clearly above zero.
    let mut distorted = reference_atoms();
    distorted[3].3 = 3.0; // F: z 1.5 -> 3.0
    let mut target = Molecule::parse(&xyz(&distorted), "xyz").expect("target");

    let rmsd = target.align_to(&reference).expect("alignment should succeed");
    assert!(rmsd > 0.1, "distinct structures should not align to 0, got {rmsd}");
}

#[test]
fn atom_count_mismatch_returns_none() {
    let reference = Molecule::parse(&xyz(&reference_atoms()), "xyz").expect("ref");
    let refs = reference_atoms();
    let mut target = Molecule::parse(&xyz(&refs[..3]), "xyz").expect("target");
    assert!(
        target.align_to(&reference).is_none(),
        "aligning molecules with different atom counts should fail"
    );
}
