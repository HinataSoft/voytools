"""Example 2: compare sections of the manuscript (Currier A vs B, labels vs paragraph text).

Usage::

    python examples/02_compare_sections.py [path/to/ZL3b-n.txt]

Shows vocabulary and glyph-entropy summaries per Currier language, the vocabulary overlap
of A and B, Dunning keyness of tokens and stems, prefix/suffix inventories per language and
how label words relate to the paragraph vocabulary (cf. Currier 1976; Zattera 2022).
"""

from __future__ import annotations

import argparse

import pandas as pd

from _paths import resolve_data, show

from voytools import filter_tokens, load_tables
from voytools.morphology import affix_inventory, decompose_table
from voytools.stats import (
    glyph_entropy_report,
    keyness,
    labels_vs_text,
    positional_profile,
    top_keyness,
    vocabulary_overlap,
    vocabulary_summary,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path", nargs="?", help="IVTFF file (default: search for ZL3b-n.txt)")
    args = ap.parse_args()

    path = resolve_data(args.path)
    tables = load_tables(path)
    print(f"File: {path}")

    # paragraph text of pages whose language is known
    text = filter_tokens(tables.tokens, generic_type="P")
    ab = text[text["language"].isin(["A", "B"])].copy()
    print(f"Paragraph text with a known language: {len(ab)} of {len(text)} clean tokens")

    show("Vocabulary per language", vocabulary_summary(ab, by="language").set_index("language"))

    rows = []
    for lang, sub in ab.groupby("language", sort=True):
        r = glyph_entropy_report(sub, "extended", scope="line").iloc[0]
        rows.append({"language": lang, "n_glyphs": int(r["n_glyphs"]), "H1": r["H1"], "H2": r["H2"], "H3": r["H3"]})
    show("Glyph entropy per language (extended preset, bits)", pd.DataFrame(rows).set_index("language"))

    ov = vocabulary_overlap(ab, "language", "A", "B")
    show("Vocabulary overlap A vs B", pd.Series(ov))

    # which whole tokens, and which stems, characterise each language?
    kt = keyness(ab, "language", unit="token", min_count=20)
    show("Tokens most characteristic of each language (G2)",
         top_keyness(kt, k=10).loc[:, ["group", "item", "count_in", "count_out", "rel_in", "rel_out", "log2_ratio", "G2"]]
         .set_index(["group", "item"]))

    dec = decompose_table(ab, prefixes="reference", suffixes="reference")
    ks = keyness(dec, "language", unit="stem", min_count=20)
    show("Stems most characteristic of each language (G2)",
         top_keyness(ks, k=8).loc[:, ["group", "item", "count_in", "count_out", "log2_ratio", "G2"]].set_index(["group", "item"]))

    for side in ("prefix", "suffix"):
        inv = affix_inventory(ab, side, max_len=3, min_count=100, by="language")
        wide = inv.pivot(index="affix", columns="language", values="share").fillna(0.0)
        wide["B/A"] = wide["B"] / wide["A"].where(wide["A"] > 0)
        wide = wide.sort_values("A", ascending=False).head(12)
        show(f"Share of tokens by {side} (top 12 for A)", wide)

    prof = {lang: positional_profile(sub, which="first") for lang, sub in ab.groupby("language", sort=True)}
    show("Token-initial glyph at the start of a line: shares in A and B",
         pd.DataFrame({f"{lang}_line_first": p.shares["line_first"] for lang, p in prof.items()})
         .fillna(0.0).sort_values("A_line_first", ascending=False).head(8))

    # labels (generic type L) against paragraph text
    lab = filter_tokens(tables.tokens, generic_type="P,L")
    lt = labels_vs_text(lab)
    show("Labels vs paragraph text", pd.Series(lt))
    print(f"\n{lt['a_token_in_b'] * 100:.1f}% of the label tokens and {lt['a_type_in_b'] * 100:.1f}% of the "
          f"label types also occur in the paragraph text; only {lt['b_type_in_a'] * 100:.1f}% of the "
          "paragraph types occur as labels")


if __name__ == "__main__":
    main()
