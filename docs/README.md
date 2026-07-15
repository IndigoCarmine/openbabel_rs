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
