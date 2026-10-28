#!/usr/bin/env python3
"""Print a status table of all experiments from ``experiments/NNN_*/result.md``.

Each ``result.md`` starts with a heading, a ``Stav:`` line and a ``Souhrn:`` paragraph (see
``experiments/_template/result.md``).  Experiments without ``result.md`` are listed as ``zadáno``.
Nothing is written; run from the repository root::

    python scripts/status.py            # markdown table
    python scripts/status.py --short    # one line per experiment
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXPERIMENTS = ROOT / "experiments"


def read_result(folder: Path) -> dict[str, str]:
    """Return ``{'status', 'date', 'summary'}`` parsed from ``folder/result.md`` (or 'zadáno')."""
    path = folder / "result.md"
    if not path.exists():
        return {"status": "zadáno", "date": "", "summary": ""}
    text = path.read_text(encoding="utf-8")
    status, date = "?", ""
    m = re.search(r"^\**Stav:\**\s*(.+)$", text, flags=re.MULTILINE)
    if m:
        parts = [p.strip() for p in re.split(r"\s·\s|\s\|\s", m.group(1))]
        status = parts[0].strip("*` ")
        for p in parts[1:]:
            if re.fullmatch(r"\d{4}-\d{2}-\d{2}", p):
                date = p
    m = re.search(r"^\**Souhrn:\**\s*(.+?)(?:\n\s*\n|\Z)", text, flags=re.MULTILINE | re.DOTALL)
    summary = " ".join(m.group(1).split()) if m else ""
    return {"status": status, "date": date, "summary": summary}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--short", action="store_true", help="one line per experiment, no table")
    ap.add_argument("--width", type=int, default=160, help="truncate summaries to this many characters (0 = no limit)")
    args = ap.parse_args(argv)

    folders = sorted(p for p in EXPERIMENTS.iterdir() if p.is_dir() and re.match(r"\d{3}_", p.name))
    rows = []
    for folder in folders:
        num, _, name = folder.name.partition("_")
        rows.append((num, name, read_result(folder)))

    if args.short:
        for num, name, r in rows:
            print(f"{num} {name}: {r['status']} {r['date']}".rstrip())
        return 0

    def cut(s: str) -> str:
        return s if not args.width or len(s) <= args.width else s[: args.width - 1] + "…"

    print("| # | Název | Stav | Datum | Souhrn |")
    print("|---|---|---|---|---|")
    for num, name, r in rows:
        summary = cut(r["summary"]).replace("|", "\\|")
        print(f"| {num} | {name} | {r['status']} | {r['date']} | {summary} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
