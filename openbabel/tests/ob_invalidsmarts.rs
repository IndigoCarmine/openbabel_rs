//! Port of OpenBabel's `test/invalidsmarts.cpp` — invalid SMARTS must be
//! rejected.
//!
//! Every line of `invalid-smarts.txt` must fail to compile, and each of the
//! three "random garbage" files (`random`, `random2`, `random3`), read whole
//! and newline-stripped into one pattern, must also fail to compile.

mod common;

use common::{ob_test_file, read_test_file};
use openbabel::SmartsPattern;

fn is_rejected(pattern: &str) -> bool {
    SmartsPattern::new(pattern).is_err()
}

#[test]
fn invalid_smarts_patterns_rejected() {
    let content = read_test_file("invalid-smarts.txt");

    let mut lines = 0usize;
    for line in content.lines() {
        lines += 1;
        assert!(
            is_rejected(line),
            "line {lines} compiled but should be an invalid SMARTS: {line:?}"
        );
    }

    assert!(lines > 0, "invalid-smarts.txt was empty");
}

#[test]
fn random_data_rejected() {
    for name in ["random", "random2", "random3"] {
        // The files may hold non-UTF-8 bytes; read raw and decode lossily (the
        // content is invalid SMARTS regardless of the exact bytes). The C++
        // getline-concatenates every line, dropping the '\n' separators.
        let bytes = std::fs::read(ob_test_file(name)).expect("read random file");
        let joined: String = String::from_utf8_lossy(&bytes).split('\n').collect();
        assert!(
            is_rejected(&joined),
            "random file {name} compiled but should be invalid SMARTS"
        );
    }
}
