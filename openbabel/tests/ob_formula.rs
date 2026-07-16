//! Port of OpenBabel's `test/formula.cpp` — molecular formula, weight and exact
//! mass.
//!
//! Reads the same corpus (`attype.00.smi`) and reference (`formularesults.txt`,
//! one `"formula molwt exactmass"` line per molecule). For each molecule it
//! checks the three values, then makes hydrogens explicit and checks that all
//! three are unchanged (regression PR#1485580). Tolerances match the C++: 1e-3
//! for the two masses; the formula compares exactly.

mod common;

use common::{read_test_file, tokenize};
use openbabel::Molecule;

#[test]
fn formula_weight_and_exact_mass_match_reference() {
    let smiles = read_test_file("attype.00.smi");
    let results = read_test_file("formularesults.txt");
    let mut ref_lines = results.lines();

    let mut molecules = 0usize;
    for line in smiles.lines() {
        let mut mol = match Molecule::parse(line, "smi") {
            Ok(m) if m.num_atoms() > 0 => m,
            _ => continue,
        };
        molecules += 1;
        let title = mol.title();

        let ref_line = ref_lines
            .next()
            .expect("ran out of reference data (formula)");
        let vs = tokenize(ref_line);
        assert_eq!(vs.len(), 3, "reference data has incorrect format");
        let ref_formula = vs[0];
        let ref_molwt: f64 = vs[1].parse().expect("molwt");
        let ref_exact: f64 = vs[2].parse().expect("exact mass");

        // Before adding explicit hydrogens.
        assert_eq!(mol.formula(), ref_formula, "formula for {title:?}");
        assert!(
            (mol.molar_mass() - ref_molwt).abs() <= 1.0e-3,
            "molar mass for {title:?}: expected {ref_molwt}, got {}",
            mol.molar_mass()
        );
        assert!(
            (mol.exact_mass() - ref_exact).abs() <= 1.0e-3,
            "exact mass for {title:?}: expected {ref_exact}, got {}",
            mol.exact_mass()
        );

        // After making implicit hydrogens explicit — must be identical.
        mol.add_hydrogens();
        assert_eq!(
            mol.formula(),
            ref_formula,
            "formula after add_hydrogens for {title:?}"
        );
        assert!(
            (mol.molar_mass() - ref_molwt).abs() <= 1.0e-3,
            "molar mass after add_hydrogens for {title:?}"
        );
        assert!(
            (mol.exact_mass() - ref_exact).abs() <= 1.0e-3,
            "exact mass after add_hydrogens for {title:?}"
        );
    }

    assert!(molecules > 0, "no molecules were read from attype.00.smi");
}
