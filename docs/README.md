# openbabel-rs guide (mdBook, EN + JA)

This directory holds the prose guide for `openbabel-rs`, built with
[mdBook](https://rust-lang.github.io/mdBook/). It is authored in **English**
under `src/`; **Japanese** is produced by the
[mdbook-i18n-helpers](https://github.com/google/mdbook-i18n-helpers) gettext
workflow.

## Prerequisites

```sh
cargo install mdbook mdbook-i18n-helpers
```

`mdbook-i18n-helpers` provides the pure-Rust `mdbook-gettext` (applies a
translation at build time) and `mdbook-xgettext` (extracts messages). Updating a
`.po` against a regenerated `.pot` uses GNU gettext's `msgmerge` — available on
Linux/macOS and in CI; on Windows it is easiest to run the merge step in CI.

## Build

```sh
# English (source language) -> book/
mdbook build

# Japanese -> book/ja/
MDBOOK_BOOK__LANGUAGE=ja mdbook build -d book/ja

# Live preview (English)
mdbook serve
```

The site is assembled for GitHub Pages as:

```
/            English book
/ja/         Japanese book
/api/        rustdoc API reference (cargo doc --workspace --no-deps)
```

## Translation workflow (gettext)

1. **Write/edit English** in `src/*.md`.
2. **Extract messages** to `po/messages.pot`:

   ```sh
   MDBOOK_OUTPUT__XGETTEXT__ENABLE=true mdbook build -d po
   ```

3. **Create or update the Japanese catalog** `po/ja.po`:

   ```sh
   # first time
   msginit -i po/messages.pot -l ja -o po/ja.po --no-translator
   # subsequently, merge new/changed messages
   msgmerge --update po/ja.po po/messages.pot
   ```

4. **Translate**: fill each `msgstr ""` in `po/ja.po`. When an English message
   changes, gettext marks the entry `#, fuzzy` — review and refresh those.
5. **Build Japanese** (step above) and check `book/ja/`.

Because translations are keyed by the exact English source text, an edit to the
English automatically flags the corresponding Japanese as stale, so the two
never silently drift apart.

## Japanese API reference (`api-i18n/`)

rustdoc has no localization, so the Japanese API reference is built by
translating the crate's **doc comments** and running `cargo doc` on the result.
The English doc comments in `openbabel/src/*.rs` are authoritative; Japanese
lives in `api-i18n/ja.json` (keyed by the whitespace-normalized English
paragraph). Fenced code blocks — including doctests — are never translated.

```sh
# Build the Japanese API reference into target/ja-api/ (English sources are
# rewritten on a copy and restored on exit; openbabel-sys stays cached).
bash docs/api-i18n/build-ja-api.sh target/ja-api

# See which paragraphs are not yet translated (English will show through):
python docs/api-i18n/doci18n.py stats openbabel/src docs/api-i18n/ja.json
```

**Keeping it in sync.** When you add or change a doc comment, `stats` lists the
new/changed English paragraphs; add an entry for each to `api-i18n/ja.json`
(key = the exact English paragraph, value = the Japanese). Anything left
untranslated simply falls back to English, so the reference is never broken by a
missing translation.

The English (`/api/`) and Japanese (`/ja/api/`) references are deployed together
by `docs.yml`, with a 🌐 switch (injected via `api-i18n/rustdoc-lang.html`) in
the top-right corner of every page.
