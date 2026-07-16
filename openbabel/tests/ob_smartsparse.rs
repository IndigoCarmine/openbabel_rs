//! Port of OpenBabel's `test/smartsparse.cpp` — valid SMARTS must compile.
//!
//! Every non-comment line of `validsmarts.txt` must compile. As in the C++,
//! only lines whose first character is `#` are skipped (so, per that test,
//! patterns like `#[C~C]` are treated as comments and not parsed).

mod common;

use common::read_test_file;
use openbabel::SmartsPattern;

#[test]
fn valid_smarts_patterns_compile() {
    let content = read_test_file("validsmarts.txt");

    let mut patterns = 0usize;
    for line in content.lines() {
        if line.starts_with('#') {
            continue; // comment line
        }
        patterns += 1;
        assert!(
            SmartsPattern::new(line).is_ok(),
            "valid SMARTS failed to compile: {line:?}"
        );
    }

    assert!(patterns > 0, "no patterns were read from validsmarts.txt");
}
