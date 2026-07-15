//! Demo CLI for the `openbabel` Rust bindings.
//!
//! Usage:
//!   openbabel-demo [SMILES]
//!
//! Parses the given SMILES (default: ethanol "CCO"), then prints core
//! properties and the canonical SMILES — exercising the MVP binding surface.

use openbabel::{Molecule, SmartsPattern};

fn main() {
    let smiles = std::env::args().nth(1).unwrap_or_else(|| "CCO".to_string());

    println!("OpenBabel version: {}", openbabel::version());
    println!("Input SMILES:      {smiles}");

    let mut mol = match Molecule::parse(&smiles, "smi") {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };

    println!("Formula:           {}", mol.formula());
    println!("Molar mass:        {:.3} g/mol", mol.molar_mass());
    println!("Exact mass:        {:.4}", mol.exact_mass());
    println!("Total charge:      {}", mol.total_charge());
    println!("Heavy atoms:       {}", mol.num_atoms());
    println!("Bonds:             {}", mol.num_bonds());

    // Descriptors (T2).
    if let Some(v) = mol.logp() {
        println!("logP:              {v:.3}");
    }
    if let Some(v) = mol.tpsa() {
        println!("TPSA:              {v:.2}");
    }

    // A SMARTS query (T2): count hydroxyl groups.
    if let Ok(oh) = SmartsPattern::new("[OX2H]") {
        println!("Hydroxyl groups:   {}", oh.num_matches(&mol));
    }

    mol.add_hydrogens();
    println!("Atoms (with H):    {}", mol.num_atoms());

    match mol.write("can") {
        Ok(can) => println!("Canonical SMILES:  {}", can.trim()),
        Err(e) => eprintln!("canonical SMILES failed: {e}"),
    }

    println!("\nAtoms:");
    for atom in mol.atoms() {
        let (x, y, z) = atom.coords();
        println!(
            "  #{:<3} Z={:<3} aromatic={:<5} ring={:<5} coords=({x:.3}, {y:.3}, {z:.3})",
            atom.index(),
            atom.atomic_number(),
            atom.is_aromatic(),
            atom.is_in_ring(),
        );
    }
}
