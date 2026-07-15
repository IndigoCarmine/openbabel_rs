#!/usr/bin/env bash
#
# Build the Japanese API reference for the `openbabel` crate.
#
# rustdoc has no localization, so we translate the crate's doc comments on a
# *throwaway copy* of the sources and run `cargo doc` on that. The English
# sources are backed up and restored on exit (even on failure), independent of
# git state. Because only `openbabel/src` changes, `openbabel-sys` (the slow,
# builds-OpenBabel-from-source crate) stays cached and the doc build is fast.
#
# Output: target/doc/ contains the Japanese `openbabel` docs; this script copies
# them to the directory given as $1 (default: target/ja-api).
#
# NOTE: while this runs it briefly rewrites openbabel/src/*.rs and restores them.
# In CI (fresh checkout) that is perfectly safe; locally, avoid editing the crate
# while it runs.
set -euo pipefail

cd "$(dirname "$0")/../.."          # repo root
OUT="${1:-target/ja-api}"
JA="docs/api-i18n/ja.json"

PY="$(command -v python3 || command -v python)"

HDR="$PWD/docs/api-i18n/rustdoc-lang.html"
if command -v cygpath >/dev/null 2>&1; then HDR="$(cygpath -m "$HDR")"; fi

BK="$(mktemp -d)"
cp -r openbabel/src/. "$BK/"
restore() { cp -r "$BK/." openbabel/src/; rm -rf "$BK"; }
trap restore EXIT

"$PY" docs/api-i18n/doci18n.py inject openbabel/src "$JA"

RUSTDOCFLAGS="--html-in-header $HDR" cargo doc -p openbabel --no-deps

rm -rf "$OUT"
mkdir -p "$OUT"
cp -r target/doc/. "$OUT/"
echo "Japanese API reference written to $OUT"
