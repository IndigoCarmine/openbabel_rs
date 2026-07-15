//! Integration tests for 2D depiction / SVG rendering (T6).

use openbabel::{Molecule, SvgOptions};

#[test]
fn to_svg_produces_svg_document() {
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let svg = mol.to_svg().expect("SVG rendering should succeed");
    assert!(svg.contains("<svg"), "should contain an <svg> element");
    assert!(svg.contains("</svg>"), "should be a complete SVG document");
    assert!(
        svg.len() > 200,
        "expected a non-trivial SVG, got {} bytes",
        svg.len()
    );
}

#[test]
fn generate_2d_sets_dimension_to_two() {
    let mut mol = Molecule::parse("CCO", "smi").expect("ethanol");
    assert_eq!(mol.dimension(), 0, "SMILES starts with no coordinates");
    assert!(!mol.has_2d());

    assert!(mol.generate_2d(), "2D generation should succeed");
    assert_eq!(mol.dimension(), 2);
    assert!(mol.has_2d());
}

#[test]
fn atom_index_option_annotates_output() {
    let mol = Molecule::parse("CCO", "smi").expect("ethanol");
    let plain = mol.to_svg().expect("svg");
    let indexed = mol
        .to_svg_with(SvgOptions {
            atom_indices: true,
            ..Default::default()
        })
        .expect("svg");
    assert_ne!(plain, indexed, "atom index labels should change the SVG");
    assert!(
        indexed.len() > plain.len(),
        "index labels should add markup"
    );
}

#[test]
fn all_carbons_option_changes_output() {
    // Benzene's carbons are all non-terminal, so the default labels none of
    // them while `all_carbons` labels all six — a clear, deterministic diff.
    let mol = Molecule::parse("c1ccccc1", "smi").expect("benzene");
    let plain = mol.to_svg().expect("svg");
    let all_c = mol
        .to_svg_with(SvgOptions {
            all_carbons: true,
            ..Default::default()
        })
        .expect("svg");
    assert_ne!(plain, all_c, "drawing all carbons should change the SVG");
}
