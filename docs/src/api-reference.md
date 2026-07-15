# API Reference

This guide explains the concepts and architecture. For the **per-item
reference** — every type, method, argument, and return value — use the
`rustdoc`-generated API documentation.

👉 **[Open the API reference](https://indigocarmine.github.io/openbabel_rs/api/openbabel/index.html)**

The reference is generated directly from the doc comments in the source, so it
is always in sync with the code. It is deployed alongside this book at
[`/api/`](https://indigocarmine.github.io/openbabel_rs/api/openbabel/index.html).

## Generating it locally

```sh
cargo doc --workspace --no-deps --open
```

`--no-deps` documents only the workspace crates (not their dependencies). The
output lands in `target/doc/`; `--open` launches it in your browser. The first
run builds OpenBabel from source and is slow — see
[Building & Platform Notes](./building.md).

## Where to start in the reference

- **[`Molecule`](https://indigocarmine.github.io/openbabel_rs/api/openbabel/struct.Molecule.html)** —
  the central type; almost everything hangs off it.
- **[`Atom`](https://indigocarmine.github.io/openbabel_rs/api/openbabel/struct.Atom.html)**
  and **[`Bond`](https://indigocarmine.github.io/openbabel_rs/api/openbabel/struct.Bond.html)** —
  borrowed views into a molecule.
- **[`Minimizer`](https://indigocarmine.github.io/openbabel_rs/api/openbabel/struct.Minimizer.html)**
  / **[`Constraints`](https://indigocarmine.github.io/openbabel_rs/api/openbabel/struct.Constraints.html)** —
  force-field minimization (see [Key Concepts](./concepts.md#constrained-minimization)).

> **A note on language.** Both the prose guide (this book) and the API reference
> are available in English and Japanese. rustdoc has no built-in language toggle,
> so the Japanese API reference is a separately built site (its doc comments are
> translated from the English source before `cargo doc` runs); a 🌐 button in the
> top-right corner of every page switches between them. You can also go straight
> to the Japanese reference at
> <https://indigocarmine.github.io/openbabel_rs/ja/api/openbabel/index.html>.
