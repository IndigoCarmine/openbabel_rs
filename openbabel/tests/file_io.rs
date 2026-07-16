//! Integration tests for file-path I/O (T16).

use openbabel::{write_many, Molecule};
use std::path::PathBuf;

/// A unique temp path for this test process (avoids cross-run collisions).
fn tmp(name: &str) -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("obrs_{}_{}", std::process::id(), name));
    p
}

#[test]
fn write_then_read_roundtrip_by_extension() {
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    let path = tmp("ethanol.mol");
    let ps = path.to_str().unwrap();

    // Format inferred from the ".mol" extension on both sides.
    mol.write_file(ps, None).expect("write .mol");
    let back = Molecule::read_file(ps, None).expect("read .mol");
    assert_eq!(back.formula(), "C2H6O");

    std::fs::remove_file(&path).ok();
}

#[test]
fn explicit_format_overrides_extension() {
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    let path = tmp("data.txt"); // extension the reader/writer would not guess as SMILES
    let ps = path.to_str().unwrap();

    mol.write_file(ps, Some("smi")).expect("write smi");
    let back = Molecule::read_file(ps, Some("smi")).expect("read smi");
    assert_eq!(back.formula(), "C2H6O");

    std::fs::remove_file(&path).ok();
}

#[test]
fn read_many_from_multi_record_file() {
    let mols = vec![
        Molecule::parse("CCO", "smi").expect("parse"),
        Molecule::parse("c1ccccc1", "smi").expect("parse"),
    ];
    let sdf = write_many(&mols, "sdf").expect("write_many");
    let path = tmp("multi.sdf");
    std::fs::write(&path, sdf).expect("write file");
    let ps = path.to_str().unwrap();

    let read = Molecule::read_file_many(ps, None).expect("read many");
    assert_eq!(read.len(), 2);
    assert_eq!(read[0].formula(), "C2H6O");
    assert_eq!(read[1].formula(), "C6H6");

    std::fs::remove_file(&path).ok();
}

#[test]
fn missing_file_is_an_error() {
    let ps = tmp("does_not_exist.mol");
    let s = ps.to_str().unwrap();
    assert!(Molecule::read_file(s, None).is_err());
    assert!(Molecule::read_file_many(s, None).is_err());
}

#[test]
fn unresolvable_format_is_an_error() {
    let mol = Molecule::parse("CCO", "smi").expect("parse");
    // No extension and no explicit format → the format cannot be resolved.
    let path = tmp("no_extension_here");
    assert!(mol.write_file(path.to_str().unwrap(), None).is_err());
    std::fs::remove_file(&path).ok();
}
