#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Doc-comment i18n for the `openbabel` crate.

rustdoc has no native localization, so we translate the crate's doc comments
themselves: the English source is authoritative, Japanese lives in a catalog
(ja.json, keyed by the whitespace-normalized English paragraph), and this tool
rewrites a *copy* of the source with the Japanese before `cargo doc` runs.

Fenced code blocks inside doc comments (```...```) are passed through verbatim,
so doctests and example code stay exactly as written.

Usage:
  doci18n.py extract <src_dir> [<src_dir> ...]      # dump unique paragraphs -> stdout JSON
  doci18n.py inject  <src_dir> <ja.json>            # translate .rs files in place (on a copy!)
  doci18n.py stats   <src_dir> <ja.json>            # coverage report
"""
import json
import re
import sys
from pathlib import Path

MARKER = re.compile(r"^(\s*)(///|//!)(.*)$")
BLOCKISH = re.compile(r"^(-|\*|>|#{1,6}\s|\d+\.)")


def split_marker(line):
    m = MARKER.match(line)
    if not m:
        return None
    indent, marker, rest = m.group(1), m.group(2), m.group(3)
    if rest.startswith(" "):
        rest = rest[1:]
    return indent, marker, rest


def norm(s):
    return " ".join(s.split())


def parse_file(path):
    """Return a list of segments for the file, each:
    ('code'|'blank'|'para'|'raw', payload). 'para' payload is (indent, marker, text).
    'raw' is a verbatim source line (non-doc or code-fence line)."""
    lines = Path(path).read_text(encoding="utf-8").split("\n")
    segs = []
    i, n = 0, len(lines)
    in_fence = False
    cur = None  # (indent, marker, [texts])

    def flush():
        nonlocal cur
        if cur:
            indent, marker, texts = cur
            segs.append(("para", (indent, marker, " ".join(texts).strip())))
            cur = None

    while i < n:
        sm = split_marker(lines[i])
        if sm is None:
            flush()
            segs.append(("raw", lines[i]))
            i += 1
            continue
        indent, marker, text = sm
        stripped = text.strip()
        if in_fence:
            segs.append(("raw", lines[i]))
            if stripped.startswith("```"):
                in_fence = False
            i += 1
            continue
        if stripped.startswith("```"):
            flush()
            segs.append(("raw", lines[i]))
            in_fence = True
            i += 1
            continue
        if stripped == "":
            flush()
            segs.append(("blank", (indent, marker)))
            i += 1
            continue
        if BLOCKISH.match(stripped):
            # list item / heading / blockquote: translate as its own single-line paragraph
            flush()
            segs.append(("para", (indent, marker, stripped)))
            i += 1
            continue
        if cur is None:
            cur = (indent, marker, [stripped])
        else:
            cur[2].append(stripped)
        i += 1
    flush()
    return segs


def iter_rs(dirs):
    for d in dirs:
        for p in sorted(Path(d).rglob("*.rs")):
            yield p


def cmd_extract(dirs):
    seen = {}
    order = []
    for p in iter_rs(dirs):
        for kind, payload in parse_file(p):
            if kind == "para":
                _, _, text = payload
                k = norm(text)
                if k and k not in seen:
                    seen[k] = str(p)
                    order.append(text)
    print(json.dumps(order, ensure_ascii=False, indent=1))
    print(f"# {len(order)} unique paragraphs", file=sys.stderr)


def wrap_doc(marker, indent, text, width=76):
    """Emit a translated paragraph as marker lines, wrapping on spaces but never
    inside `code`/[links]. Simple greedy word wrap; Japanese has few spaces so it
    mostly ends up on one line, which rustdoc renders fine."""
    words = text.split(" ")
    out, line = [], ""
    for w in words:
        cand = w if not line else line + " " + w
        if len(cand) > width and line:
            out.append(f"{indent}{marker} {line}".rstrip())
            line = w
        else:
            line = cand
    out.append(f"{indent}{marker} {line}".rstrip())
    return out


def cmd_inject(dirs, ja_path):
    ja = json.loads(Path(ja_path).read_text(encoding="utf-8"))
    ja = {norm(k): v for k, v in ja.items()}
    total = hit = 0
    for p in iter_rs(dirs):
        segs = parse_file(p)
        out = []
        for kind, payload in segs:
            if kind == "raw":
                out.append(payload)
            elif kind == "blank":
                indent, marker = payload
                out.append(f"{indent}{marker}".rstrip())
            elif kind == "para":
                indent, marker, text = payload
                total += 1
                tr = ja.get(norm(text))
                if tr:
                    hit += 1
                    out.extend(wrap_doc(marker, indent, tr))
                else:
                    out.append(f"{indent}{marker} {text}".rstrip())
        Path(p).write_text("\n".join(out), encoding="utf-8")
    print(f"injected {hit}/{total} paragraphs", file=sys.stderr)


def cmd_stats(dirs, ja_path):
    ja = {norm(k): v for k, v in json.loads(Path(ja_path).read_text(encoding="utf-8")).items()}
    total = hit = 0
    missing = []
    for p in iter_rs(dirs):
        for kind, payload in parse_file(p):
            if kind == "para":
                _, _, text = payload
                total += 1
                if norm(text) in ja:
                    hit += 1
                else:
                    missing.append((str(p), text))
    print(f"coverage: {hit}/{total}")
    for f, t in missing:
        print("  MISSING", f, ascii(t[:70]))


def main():
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8")
        except Exception:
            pass
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    cmd = sys.argv[1]
    if cmd == "extract":
        cmd_extract(sys.argv[2:])
    elif cmd == "inject":
        cmd_inject(sys.argv[2:-1], sys.argv[-1])
    elif cmd == "stats":
        cmd_stats(sys.argv[2:-1], sys.argv[-1])
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
