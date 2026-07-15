//! Integration tests for reaction I/O (T14): the `Reaction` type.

use openbabel::{Molecule, Reaction};
use std::collections::HashSet;

#[test]
fn parse_reaction_smiles_counts() {
    // Ethene + water -> ethanol.
    let rxn = Reaction::parse("C=C.O>>CCO", "rsmi").expect("parse rsmi");
    assert_eq!(rxn.num_reactants(), 2);
    assert_eq!(rxn.num_products(), 1);
    assert_eq!(rxn.num_agents(), 0);
}

#[test]
fn extracted_components_are_real_molecules() {
    let rxn = Reaction::parse("C=C.O>>CCO", "rsmi").expect("parse rsmi");

    let reactant_formulas: HashSet<String> = (0..rxn.num_reactants())
        .map(|i| rxn.reactant(i).unwrap().formula())
        .collect();
    assert!(reactant_formulas.contains("C2H4"), "{reactant_formulas:?}");
    assert!(reactant_formulas.contains("H2O"), "{reactant_formulas:?}");

    assert_eq!(rxn.product(0).unwrap().formula(), "C2H6O");
    // Out-of-range access yields None.
    assert!(rxn.product(5).is_none());
    assert!(rxn.reactant(2).is_none());
}

#[test]
fn reaction_smiles_with_agents() {
    // Ethene + H2, palladium as an agent, -> ethane.
    let rxn = Reaction::parse("C=C.[H][H]>[Pd]>CC", "rsmi").expect("parse rsmi");
    assert_eq!(rxn.num_reactants(), 2);
    assert_eq!(rxn.num_agents(), 1);
    assert_eq!(rxn.num_products(), 1);
}

#[test]
fn round_trip_reaction_smiles() {
    let rxn = Reaction::parse("C=C.O>>CCO", "rsmi").expect("parse");
    let text = rxn.write("rsmi").expect("write rsmi");
    assert!(text.contains(">>"), "not a reaction SMILES: {text:?}");

    // Re-parsing the written form preserves the component counts.
    let rxn2 = Reaction::parse(text.trim(), "rsmi").expect("re-parse");
    assert_eq!(rxn2.num_reactants(), 2);
    assert_eq!(rxn2.num_products(), 1);
}

#[test]
fn build_reaction_programmatically() {
    let ethene = Molecule::parse("C=C", "smi").expect("ethene");
    let ethanol = Molecule::parse("CCO", "smi").expect("ethanol");

    let mut rxn = Reaction::new();
    rxn.add_reactant(&ethene).add_product(&ethanol);
    assert_eq!(rxn.num_reactants(), 1);
    assert_eq!(rxn.num_products(), 1);

    let text = rxn.write("rsmi").expect("write");
    assert!(text.contains(">>"), "{text:?}");
}

#[test]
fn title_comment_and_reversible() {
    let mut rxn = Reaction::new();
    assert!(rxn.title().is_empty());
    assert!(!rxn.is_reversible());

    rxn.set_title("esterification")
        .set_comment("acid-catalysed")
        .set_reversible(true);
    assert_eq!(rxn.title(), "esterification");
    assert_eq!(rxn.comment(), "acid-catalysed");
    assert!(rxn.is_reversible());
}

#[test]
fn unknown_format_is_rejected() {
    assert!(Reaction::parse("C=C>>CC", "not-a-format").is_err());
    let rxn = Reaction::new();
    assert!(rxn.write("not-a-format").is_err());
}
