"""Example 1: overview of a transliteration (vocabulary, Zipf, lengths, glyph entropy).

Usage::

    python examples/01_summary.py [path/to/ZL3b-n.txt] [--plot zipf.png]

Without a path the file is searched as described in ``_paths.py``.
"""

from __future__ import annotations

import argparse

import pandas as pd

from _paths import resolve_data, show

from voytools import filter_tokens, load_tables
from voytools.plots import have_matplotlib, plot_zipf
from voytools.stats import (
    glyph_entropy_report,
    length_distribution,
    token_frequency,
    vocabulary_summary,
    zipf_fit,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path", nargs="?", help="IVTFF file (default: search for ZL3b-n.txt)")
    ap.add_argument("--plot", metavar="PNG", help="save a rank-frequency plot (needs matplotlib)")
    args = ap.parse_args()

    path = resolve_data(args.path)
    tables = load_tables(path)
    print(f"File: {path}")
    print(f"Alphabet {tables.header.alphabet}, IVTFF {tables.header.version}: "
          f"{len(tables.pages)} pages, {len(tables.loci)} loci, {len(tables.tokens)} tokens, "
          f"{len(tables.warnings)} parser warnings")

    show("Loci by generic type", tables.loci.groupby("generic_type").size().rename("loci"))

    # paragraph text only, without tokens that contain unreadable / alternative / high-ASCII characters
    text = filter_tokens(tables.tokens, generic_type="P")
    print(f"\nParagraph text: {len(text)} clean tokens")

    vs = pd.concat([
        vocabulary_summary(filter_tokens(tables.tokens, generic_type="all", clean=True)).assign(selection="all loci"),
        vocabulary_summary(text).assign(selection="paragraph text"),
    ]).set_index("selection")
    show("Vocabulary", vs)
    show("Vocabulary by language (paragraph text)", vocabulary_summary(text, by="language").set_index("language"))

    freq = token_frequency(text)
    show("Most frequent tokens", freq.head(15).set_index("rank"))
    fit = zipf_fit(freq["count"], min_rank=1, max_rank=500)
    print(f"\nZipf fit on ranks 1-500: slope {fit['slope']:.3f}, R^2 {fit['r2']:.3f}")

    show("Token length distribution (characters)", length_distribution(text, unit="chars").head(14).set_index("length"))
    print(f"Mean token length: {text['len_chars'].mean():.3f} characters, {text['len_glyphs'].mean():.3f} glyphs (basic preset)")

    # glyph entropies of the text without spaces; 'line' keeps bigrams inside lines, 'all' is one long stream
    rows = []
    for preset in ("basic", "extended"):
        for scope in ("line", "all"):
            r = glyph_entropy_report(text, preset, scope=scope).iloc[0]
            rows.append({"preset": preset, "scope": scope, **r[["n_glyphs", "n_symbols", "H1", "H2", "H3"]]})
    ent = pd.DataFrame(rows).set_index(["preset", "scope"])
    ent[["n_glyphs", "n_symbols"]] = ent[["n_glyphs", "n_symbols"]].astype(int)
    show("Glyph entropy of the space-less text (bits)", ent)

    if args.plot:
        if have_matplotlib():
            plot_zipf({"all paragraph tokens": freq}, path=args.plot)
            print(f"\nwrote {args.plot}")
        else:
            print("\nmatplotlib is not installed: pip install voytools[plot]")


if __name__ == "__main__":
    main()
