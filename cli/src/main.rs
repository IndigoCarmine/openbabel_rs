//! Demo CLI for the `openbabel` Rust bindings.
//!
//! Usage:
//!   openbabel-demo [SMILES] [SVG_OUTPUT_PATH]
//!
//! Parses the given SMILES (default: ethanol "CCO"), then prints core
//! properties and the canonical SMILES — exercising the MVP binding surface.
//! If a second argument is given, the 2D depiction is written there as SVG.

use openbabel::{Algorithm, Minimizer, Molecule, SmartsPattern};

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
    println!("Rings (SSSR):      {}", mol.num_rings());
    println!("Rotatable bonds:   {}", mol.num_rotatable_bonds());

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

    // Stereochemistry (T7): perceived stereocenters.
    println!(
        "Stereocenters:     {} tetrahedral, {} cis/trans",
        mol.tetrahedral_stereo_count(),
        mol.cistrans_stereo_count(),
    );

    // Rings (T11): SSSR sizes and aromaticity ("a" = aromatic).
    let rings: Vec<String> = mol
        .rings()
        .map(|r| format!("{}{}", r.size(), if r.is_aromatic() { "a" } else { "" }))
        .collect();
    if !rings.is_empty() {
        println!("Ring sizes:        [{}]", rings.join(", "));
    }

    // Graph navigation (T12): describe the first heavy atom's neighbourhood.
    if let Some(a0) = mol.atom(0) {
        let neighbors: Vec<String> = a0
            .neighbors()
            .iter()
            .map(|n| format!("#{}(Z{})", n.index(), n.atomic_number()))
            .collect();
        println!(
            "Atom #0 neighbours: {} bond(s) -> [{}]",
            a0.bonds().len(),
            neighbors.join(", "),
        );
    }

    // Topological symmetry (T13): how many symmetry-distinct atom environments.
    {
        use std::collections::HashSet;
        let classes = mol.symmetry_classes();
        let distinct = classes.iter().collect::<HashSet<_>>().len();
        println!(
            "Symmetry classes:  {distinct} distinct among {} atoms",
            classes.len(),
        );
    }

    // Graph automorphisms (T15): the order of the molecule's symmetry group.
    println!("Automorphisms:     {}", mol.automorphisms().len());

    // 2D depiction (T6): render the skeletal structure to SVG (before adding
    // explicit H, for a cleaner drawing). A second CLI argument saves it.
    if let Some(svg) = mol.to_svg() {
        match std::env::args().nth(2) {
            Some(path) => match std::fs::write(&path, &svg) {
                Ok(()) => println!("Wrote SVG:         {path} ({} bytes)", svg.len()),
                Err(e) => eprintln!("failed to write SVG to {path}: {e}"),
            },
            None => println!(
                "2D depiction:      {} bytes of SVG (pass an output path to save)",
                svg.len()
            ),
        }
    }

    mol.add_hydrogens();
    println!("Atoms (with H):    {}", mol.num_atoms());

    match mol.write("can") {
        Ok(can) => println!("Canonical SMILES:  {}", can.trim()),
        Err(e) => eprintln!("canonical SMILES failed: {e}"),
    }

    // Structure identifiers (T4).
    if let Some(inchi) = mol.inchi() {
        println!("InChI:             {inchi}");
    }
    if let Some(key) = mol.inchikey() {
        println!("InChIKey:          {key}");
    }

    // Generate a 3D structure and report a force-field energy (T3).
    if mol.generate_3d() {
        let unit = openbabel::forcefield_energy_unit("MMFF94").unwrap_or_default();
        if let Some(e) = mol.energy("MMFF94") {
            println!("MMFF94 energy:     {e:.3} {unit}");
        }
        if let Some(e) = mol.optimize_geometry("MMFF94", 500) {
            println!("  after optimize:  {e:.3} {unit}");
        }

        // Superpose a fresh copy back onto the original (T5): identical
        // structures align with ~0 RMSD.
        if let Ok(mut copy) = Molecule::parse(&mol.write("mol").unwrap_or_default(), "mol") {
            if let Some(rmsd) = copy.align_to(&mol) {
                println!("Self-align RMSD:   {rmsd:.4}");
            }
        }

        // Conformer search (T7): count diverse low-energy conformers.
        println!("Conformers:        {}", mol.generate_conformers(10));

        // Score every conformer at once (T17): report the lowest energy found.
        if let Some(lowest) = mol
            .conformer_energies("MMFF94")
            .into_iter()
            .filter(|e| e.is_finite())
            .reduce(f64::min)
        {
            println!("  lowest E:        {lowest:.3} {unit}");
        }

        // Spectrophore descriptor (T9): a 48-value 3D shape/property fingerprint.
        let sp = mol.spectrophore();
        if !sp.is_empty() {
            println!("Spectrophore:      {} values (first {:+.2})", sp.len(), sp[0]);
        }

        // Configurable minimization with a trajectory (T10): steepest descent,
        // capturing the energy every few steps.
        let mut cfg = Minimizer::new("MMFF94");
        cfg.algorithm(Algorithm::SteepestDescent)
            .max_steps(100)
            .steps_per_frame(20);
        let traj: Vec<_> = mol.minimize(&cfg).collect();
        if let (Some(first), Some(last)) = (traj.first(), traj.last()) {
            println!(
                "Minimize (SD):     {} frames, E {:.3} -> {:.3} {unit}",
                traj.len(),
                first.energy,
                last.energy,
            );
        }
    }

    // Assign partial atomic charges so the per-atom listing can show them (T4).
    let charged = mol.compute_charges("gasteiger");
    if charged {
        println!("Charge model:      gasteiger");
    }

    println!("\nAtoms:");
    for atom in mol.atoms() {
        let (x, y, z) = atom.coords();
        println!(
            "  #{:<3} Z={:<3} q={:+.3} hyb={} arom={:<5} ring={:<5} coords=({x:.3}, {y:.3}, {z:.3})",
            atom.index(),
            atom.atomic_number(),
            atom.partial_charge(),
            atom.hybridization(),
            atom.is_aromatic(),
            atom.is_in_ring(),
        );
    }
}
