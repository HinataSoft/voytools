"""Example 3: prefix-stem-suffix morphology of the real text against two null models.

Usage::

    python examples/03_morphology_and_nulls.py [path/to/ZL3b-n.txt] [--seed N] [--shuffles N]

Steps
-----
1. learn prefix/suffix lists from the paragraph text and split every token;
2. paradigm metrics (stems with several suffixes, I(prefix; suffix | stem));
3. the same metrics on a within-line shuffle and on an auto-copy text of identical size,
   decomposed with the SAME affix lists;
4. order-dependent statistics (adjacent-token MI, edge-glyph MI, neighbour similarity,
   glyph entropy of the space-less text, BPE dependence gap).

Type-level morphology is a property of the vocabulary and its frequencies, so the within-line
shuffle reproduces it exactly: only the order-dependent block in step 4 can tell the shuffle
from the real text.  The auto-copy generator (Timm & Schinner 2020 style, deliberately simple)
reproduces word order effects but only partly the morphology.
"""

from __future__ import annotations

import argparse

import pandas as pd

from _paths import resolve_data, show

from voytools import filter_tokens, load_tables
from voytools.morphology import bpe_dependence_curve, decompose_table, paradigm_metrics, paradigm_table
from voytools.nulls import autocopy_like, shuffle_within_lines
from voytools.stats import (
    adjacent_token_mi,
    edge_glyph_mi,
    glyph_entropy_report,
    neighbor_similarity,
    vocabulary_summary,
)


def lines_of(df: pd.DataFrame) -> list[list[str]]:
    """Tokens of a table as a list of lines (one list of tokens per locus)."""
    return df.groupby("locus_id", sort=False)["token"].apply(list).tolist()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path", nargs="?", help="IVTFF file (default: search for ZL3b-n.txt)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--shuffles", type=int, default=10, help="shuffles for the bias baselines (default 10)")
    args = ap.parse_args()

    path = resolve_data(args.path)
    tables = load_tables(path)
    text = filter_tokens(tables.tokens, generic_type="P")
    print(f"File: {path}\nParagraph text: {len(text)} clean tokens in {text['locus_id'].nunique()} lines")

    # 1. affix lists are learned on the real text only and reused for the null texts
    real = decompose_table(text)
    pre, suf = real.attrs["prefixes"], real.attrs["suffixes"]
    print(f"\nLearned prefixes: {' '.join(pre)}\nLearned suffixes: {' '.join(suf)}")
    print("Example splits: " + ", ".join(
        f"{t}={p}+{s}+{x}" for t, p, s, x in real.drop_duplicates("token").head(8)[["token", "prefix", "stem", "suffix"]].itertuples(index=False)))

    # 2. paradigm of the most frequent stems
    pt = paradigm_table(real, min_stem_count=100, top_stems=10, top_suffixes=8)
    show("Stem x suffix token counts (most frequent stems)", pt)

    # 3. null texts of the same size and line layout
    shuffled = decompose_table(shuffle_within_lines(text, seed=args.seed), prefixes=pre, suffixes=suf)
    # the generator learns which glyph pairs are legal from all real words (allowed_from)
    gen = autocopy_like(text, allowed_from=text["token"].tolist(), seed=args.seed)
    copied = decompose_table(gen, prefixes=pre, suffixes=suf)
    corpora = {"real": real, "shuffled": shuffled, "autocopy": copied}

    metrics: dict[str, dict[str, float]] = {}
    for name, df in corpora.items():
        vs = vocabulary_summary(df).iloc[0]
        pm = paradigm_metrics(df, min_stem_count=5, shuffles=args.shuffles, seed=args.seed).summary
        metrics[name] = {
            "types": vs["n_types"], "hapax_share": vs["hapax_share"], "mean_len_glyphs": vs["mean_len_glyphs"],
            "stems_with_5+": pm["n_stems_k"], "suffixes_per_stem": pm["mean_suffixes_per_stem"],
            "H(suffix|stem)": pm["H_suffix_given_stem"], "I(P;S|stem)": pm["cmi_prefix_suffix_given_stem"],
            "I(P;S|stem) excess": pm["cmi_excess"],
        }
    show("Morphology (paradigm) metrics - real and shuffled agree by construction (up to the random bias baseline)",
         pd.DataFrame(metrics), digits=3)

    # 4. order-dependent statistics
    order: dict[str, dict[str, float]] = {}
    for name, df in corpora.items():
        mi = adjacent_token_mi(df, shuffles=args.shuffles, seed=args.seed).iloc[0]
        em = edge_glyph_mi(df, shuffles=args.shuffles, seed=args.seed).iloc[0]
        ns = neighbor_similarity(df, window=5, shuffles=args.shuffles, seed=args.seed).iloc[0]
        h = glyph_entropy_report(df, "extended", scope="all").iloc[0]
        order[name] = {
            "token MI excess share": mi["excess_share"], "edge-glyph MI excess (bits)": em["excess_bits"],
            "near-neighbour share (window 5)": ns["share_near"], "near / shuffled": ns["ratio"],
            "H(glyph_t+1 | glyph_t), extended": h["H2"],
        }
    show("Order-dependent statistics", pd.DataFrame(order), digits=4)

    # BPE dependence gap H(U) - H(U_t+1 | U_t) of the space-less text (Rozanova & Temerev 2026)
    curves = {name: bpe_dependence_curve(lines_of(df), (0, 16, 32, 64, 128), preset="basic").set_index("merges")["gap"]
              for name, df in corpora.items()}
    show("BPE dependence gap (bits) after k merges", pd.DataFrame(curves), digits=4)

    print("\nReading the tables: the shuffle keeps the vocabulary, so every type-level number agrees with the real "
          "text, while the two mutual-information excesses fall to zero and the glyph entropy, near-neighbour "
          "share and BPE gap move towards their order-free values.  The auto-copy text shows how much of the "
          "sequence structure simple copy-and-modify dynamics can produce (too many near neighbours, too little "
          "glyph-level regularity) and that its vocabulary and paradigms differ from the real ones.")

if __name__ == "__main__":
    main()
