"""Regenerate VERIFICATION.md: numbers measured with voytools on the real transliteration files.

Usage::

    python tools/make_verification.py [DATA_DIR] [--out VERIFICATION.md] [--tests "155 passed in 20s"]

``DATA_DIR`` is the folder with ZL3b-n.txt, GC2a-n.txt, RF1b-er.txt, IT2a-n.txt and
LSI_ivtff_0d.txt (default: the same search as the examples).  Everything in the report is
computed by this script; only the commentary text is fixed.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import platform
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "examples"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _paths import candidates  # noqa: E402
from voytools import __version__, filter_tokens, load_tables  # noqa: E402
from voytools import cli  # noqa: E402
from voytools.ivtff import parse_ivtff  # noqa: E402
from voytools.morphology import affix_inventory, decompose_table, paradigm_metrics  # noqa: E402
from voytools.stats import (  # noqa: E402
    adjacent_token_mi,
    edge_glyph_mi,
    glyph_entropy_report,
    neighbor_similarity,
    positional_profile,
    vocabulary_summary,
)

SHUFFLES = 100


def md_table(df: pd.DataFrame, digits: int = 4, index: bool = True) -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table."""
    d = df.reset_index() if index else df
    cols = [str(c) for c in d.columns]

    def fmt(v) -> str:
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return "NA"
            return f"{v:.{digits}f}"
        if isinstance(v, (int, np.integer)):
            return f"{int(v):,}".replace(",", " ")
        return str(v)

    rows = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for rec in d.itertuples(index=False):
        rows.append("| " + " | ".join(fmt(v) for v in rec) + " |")
    return "\n".join(rows)


def find_data(arg: str | None) -> Path:
    if arg:
        p = Path(arg)
        if (p / "ZL3b-n.txt").exists():
            return p
        sys.exit(f"ZL3b-n.txt not found in {p}")
    for c in candidates():
        if c.exists():
            return c.parent
    sys.exit("data folder not found; pass it as the first argument")


def measure(path: Path, presets: list, transcriber: str | None = None, full: bool = True) -> dict:
    """All headline numbers for one transliteration file."""
    t0 = time.perf_counter()
    translit = parse_ivtff(path)
    tables = load_tables(path, transcriber=transcriber)
    t_load = time.perf_counter() - t0
    tok = tables.tokens
    text = filter_tokens(tok, generic_type="P")
    res: dict = {
        "name": path.name,
        "alphabet": tables.header.alphabet,
        "version": tables.header.version,
        "pages": len(tables.pages),
        "loci": len(tables.loci),
        "tokens_all": len(tok),
        "tokens_all_clean": len(filter_tokens(tok, generic_type="all")),
        "tokens_P": int((tok["generic_type"] == "P").sum()),
        "tokens_P_clean": len(text),
        "tokens_P0_clean": len(filter_tokens(tok, generic_type="P", locus_type="P0")),
        "warnings": len(tables.warnings),
        "t_load": t_load,
        "interlinear": translit.is_interlinear,
    }
    vs = vocabulary_summary(text).iloc[0]
    res.update(types=int(vs["n_types"]), hapax=int(vs["n_hapax"]), hapax_share=float(vs["hapax_share"]),
               mean_len=float(vs["mean_len"]), mean_len_glyphs=float(vs["mean_len_glyphs"]))
    ent = []
    for pr in presets:
        for scope in ("line", "all"):
            r = glyph_entropy_report(text, pr, scope=scope).iloc[0]
            ent.append({"preset": pr or "none", "scope": scope, "n_glyphs": int(r["n_glyphs"]),
                        "n_symbols": int(r["n_symbols"]), "H1": r["H1"], "H2": r["H2"], "H3": r["H3"]})
    res["entropy"] = pd.DataFrame(ent)
    t0 = time.perf_counter()
    res["mi_tok"] = adjacent_token_mi(text, shuffles=SHUFFLES, seed=0).iloc[0]
    res["mi_edge"] = edge_glyph_mi(text, shuffles=SHUFFLES, seed=0).iloc[0]
    res["t_mi"] = time.perf_counter() - t0
    if full:
        t0 = time.perf_counter()
        res["neighbor"] = neighbor_similarity(text, window=5, shuffles=SHUFFLES, seed=0).iloc[0]
        res["t_nb"] = time.perf_counter() - t0
        res["pre"] = affix_inventory(text, "prefix", max_len=3, min_count=20).head(15).reset_index(drop=True)
        res["suf"] = affix_inventory(text, "suffix", max_len=3, min_count=20).head(15).reset_index(drop=True)
        dec = decompose_table(text)
        res["learned"] = (dec.attrs["prefixes"], dec.attrs["suffixes"])
        res["paradigm"] = paradigm_metrics(dec, min_stem_count=5, shuffles=20, seed=0).summary
        pos = positional_profile(text, which="first", permutations=100, seed=0)
        res["pos_first"] = pos
        res["pos_last"] = positional_profile(text, which="last", permutations=100, seed=0)
    return res


def extra_checks(d: Path) -> dict:
    """Checks that explain differences to the task brief and test the estimators on a null."""
    import re

    from voytools.nulls import shuffle_within_lines

    out: dict = {}
    zl = d / "ZL3b-n.txt"
    rows = []
    for pol in ("split", "join"):
        tok = load_tables(zl, space_policy=pol).tokens
        for clean in (False, True):
            rows.append({
                "space policy": pol, "tokens": "clean" if clean else "unfiltered",
                "all loci": len(filter_tokens(tok, generic_type="all", clean=clean)),
                "P": len(filter_tokens(tok, generic_type="P", clean=clean)),
                "P0": len(filter_tokens(tok, generic_type="P", locus_type="P0", clean=clean)),
            })
    out["counts"] = pd.DataFrame(rows)
    # reference-style count: the tokenizer of voynich-units/analysis/reproduce_headlines.py
    locus = re.compile(r"^<(f[0-9rv]+[0-9]*)\.(\d+),([@+*=~])([A-Za-z])[^>]*>")
    bad = set("?*<>{}[]()|")
    n_all = n_ok = 0
    ref_ids = set()
    for raw in zl.open(encoding="utf-8", errors="replace"):
        m = locus.match(raw)
        if not m or m.group(4) != "P":
            continue
        ref_ids.add(f"{m.group(1)}.{m.group(2)}")
        text = raw[raw.index(">") + 1:].strip()
        text = re.sub(r"<[^>]*>", "", text)
        text = re.sub(r"\[([^:\]]*):[^\]]*\]", r"\1", text)
        text = re.sub(r"\{[^}]*\}", "", text)
        words = [w for w in re.split(r"[.,]", text) if w]
        n_all += len(words)
        n_ok += sum(1 for w in words if not set(w) & bad)
    out["reference_count"] = (n_all, n_ok)
    full = load_tables(zl).tokens
    ptok = full[full["generic_type"] == "P"]
    out["diff_breaks"] = int(ptok["after_drawing_break"].sum())
    out["diff_unmatched"] = int((~ptok["locus_id"].isin(ref_ids)).sum())
    # null check of the estimators: a within-line shuffle must remove both MI excesses
    text = filter_tokens(load_tables(zl).tokens, generic_type="P")
    sh = shuffle_within_lines(text, seed=0)
    out["shuffle_mi"] = adjacent_token_mi(sh, shuffles=20, seed=1).iloc[0]
    out["shuffle_edge"] = edge_glyph_mi(sh, shuffles=20, seed=1).iloc[0]
    out["real_h2"] = glyph_entropy_report(text, "extended", scope="all").iloc[0]["H2"]
    out["shuffle_h2"] = glyph_entropy_report(sh, "extended", scope="all").iloc[0]["H2"]
    out["vocab_same"] = vocabulary_summary(text).iloc[0].tolist() == vocabulary_summary(sh).iloc[0].tolist()
    # Currier languages
    rows = []
    for lang, sub in text[text["language"].isin(["A", "B"])].groupby("language", sort=True):
        g = glyph_entropy_report(sub, "extended", scope="line").iloc[0]
        v = vocabulary_summary(sub).iloc[0]
        rows.append({"language": lang, "tokens": int(v["n_tokens"]), "types": int(v["n_types"]),
                     "hapax share": v["hapax_share"], "H1": g["H1"], "H2": g["H2"]})
    out["languages"] = pd.DataFrame(rows)
    return out


def run_cli(name: str, path: Path, argv_extra: list[str]) -> tuple[int, float, int]:
    """Run one CLI command in-process; returns (exit code, seconds, bytes of output)."""
    buf, err = io.StringIO(), io.StringIO()
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
        code = cli.main([name, str(path), *argv_extra])
    return code, time.perf_counter() - t0, len(buf.getvalue())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("data_dir", nargs="?")
    ap.add_argument("--out", default=str(ROOT / "VERIFICATION.md"))
    ap.add_argument("--tests", default="(not recorded; run `python -m pytest`)", help="pytest summary line to include")
    args = ap.parse_args()

    d = find_data(args.data_dir)
    zl = measure(d / "ZL3b-n.txt", ["basic", "extended"])
    gc = measure(d / "GC2a-n.txt", [None])
    others = []
    for fname, tr in (("RF1b-er.txt", None), ("IT2a-n.txt", None), ("LSI_ivtff_0d.txt", "H")):
        if (d / fname).exists():
            others.append(measure(d / fname, ["basic", "extended"], transcriber=tr, full=False))

    ex = extra_checks(d)
    cli_rows = []
    for fname in ("ZL3b-n.txt", "GC2a-n.txt"):
        for cmd, extra in (("summary", []), ("entropy", []), ("mi", ["--shuffles", str(SHUFFLES)]),
                           ("positions", ["--permutations", "100"]), ("affixes", ["--top", "15"])):
            code, secs, nbytes = run_cli(cmd, d / fname, extra)
            cli_rows.append({"file": fname, "command": f"{cmd} {' '.join(extra)}".strip(),
                             "exit code": code, "seconds": secs, "output bytes": nbytes})

    L: list[str] = []
    w = L.append
    w(f"# Verification of voytools {__version__}\n")
    w("Numbers measured with `tools/make_verification.py` on the real transliteration files "
      f"(Python {platform.python_version()}, numpy {np.__version__}, pandas {pd.__version__}).  "
      "Unless stated otherwise the statistics use paragraph text (generic type `P`), clean tokens "
      "(no unreadable `?`, no alternative readings, no high-ASCII characters), the default "
      "`split` space policy, seed 0 and "
      f"{SHUFFLES} shuffles for the shuffle baselines.\n")

    w("## 1. Test suite\n")
    w(f"`{args.tests}`\n")

    w("## 2. Corpus size\n")
    rows = []
    for r in [zl, gc, *others]:
        rows.append({"file": r["name"], "alphabet": r["alphabet"], "IVTFF": r["version"], "pages": r["pages"],
                     "loci": r["loci"], "tokens (all)": r["tokens_all"], "all, clean": r["tokens_all_clean"],
                     "P": r["tokens_P"], "P, clean": r["tokens_P_clean"], "P0, clean": r["tokens_P0_clean"],
                     "parser warnings": r["warnings"]})
    w(md_table(pd.DataFrame(rows).set_index("file")))
    w("\nThe LSI file is interlinear (18 transcribers); the row shows transcriber H, the most complete one.\n")

    w("## 3. Vocabulary (paragraph text, clean tokens)\n")
    rows = [{"file": r["name"], "tokens": r["tokens_P_clean"], "types": r["types"], "hapax": r["hapax"],
             "hapax share": r["hapax_share"], "mean length (chars)": r["mean_len"],
             "mean length (glyphs)": r["mean_len_glyphs"]} for r in [zl, gc, *others]]
    w(md_table(pd.DataFrame(rows).set_index("file")))
    w("\nHapax share = hapax legomena / types.  Glyph lengths use the `basic` Eva preset (identity for v101).\n")

    w("## 4. Glyph entropy of the space-less text (bits)\n")
    w("`H(X)` = unigram entropy `H1`, `H(X_{t+1} given X_t)` = conditional entropy `H2`, `H3` = the same given two preceding glyphs.  Scope `line`: glyph pairs "
      "do not cross line breaks; scope `all`: the whole text is one stream (this is what Rozanova & Temerev "
      "do).\n")
    for r in [zl, gc]:
        w(f"### {r['name']}\n")
        w(md_table(r["entropy"], index=False))
        w("")
    w("### Same measure on the other transcriptions (extended preset, scope `all`)\n")
    rows = []
    for r in [zl, *others]:
        e = r["entropy"]
        e = e[(e["preset"] == "extended") & (e["scope"] == "all")].iloc[0]
        rows.append({"file": r["name"], "H1": e["H1"], "H2": e["H2"], "H3": e["H3"]})
    w(md_table(pd.DataFrame(rows).set_index("file")))
    w("")

    w("## 5. Mutual information\n")
    rows = []
    for r in [zl, gc, *others]:
        t = r["mi_tok"]
        rows.append({"file": r["name"], "token pairs": int(t["n_pairs"]), "H(X) capped": t["H_cap"],
                     "MI obs": t["mi_obs"], "MI shuffled": t["mi_shuffled_mean"],
                     "excess (bits)": t["excess_bits"], "excess / H(X)": t["excess_share"],
                     "excess / MI obs": t["excess_bits"] / t["mi_obs"], "z": t["z"]})
    w("### Adjacent tokens (within lines; types capped at 2000 + `<other>`)\n")
    w(f"Excess = observed MI minus the mean of {SHUFFLES} within-line shuffles.  The *excess share* used in the "
      "comparison below is excess / H(X), where H(X) is the entropy of the (capped) token distribution, i.e. the "
      "fraction of the information of a token that is explained by its left neighbour; the column "
      "`excess / MI obs` is the share of the observed plug-in MI that survives the shuffle correction.\n")
    w(md_table(pd.DataFrame(rows).set_index("file")))
    rows = []
    for r in [zl, gc, *others]:
        e = r["mi_edge"]
        rows.append({"file": r["name"], "pairs": int(e["n_pairs"]), "MI obs": e["mi_obs"],
                     "MI shuffled": e["mi_shuffled_mean"], "excess (bits)": e["excess_bits"], "z": e["z"]})
    w("\n### Edge glyphs (last glyph of token t against first glyph of token t+1)\n")
    w(md_table(pd.DataFrame(rows).set_index("file")))
    w("\n### Neighbour similarity (window 5 tokens, edit distance <= 1, shuffle = within-page)\n")
    rows = []
    for r in [zl, gc]:
        n = r["neighbor"]
        rows.append({"file": r["name"], "share near": n["share_near"], "shuffled": n["share_near_shuffled"],
                     "ratio": n["ratio"], "z": n["z"], "share identical": n["share_ident"],
                     "identical shuffled": n["share_ident_shuffled"]})
    w(md_table(pd.DataFrame(rows).set_index("file")))
    w("")

    w("## 6. Positional profile of the first glyph (ZL3b)\n")
    pf = zl["pos_first"]
    w("Share of tokens in each position class that begin with the glyph (the eight glyphs with the highest share "
      "in any class; `parfirst_*` = lines that begin a paragraph).\n")
    top = pf.shares.max(axis=1).sort_values(ascending=False).head(8).index
    vals = pf.shares.loc[top]
    tab = pd.DataFrame([[f"{v:.3f}" for v in row] for row in vals.to_numpy(dtype=float)],
                       index=vals.index, columns=vals.columns)
    tab.loc["tokens in class"] = [f"{int(pf.n[c]):,}".replace(",", " ") for c in tab.columns]
    w(md_table(tab.rename_axis("glyph")))
    w("\nJensen-Shannon divergence between classes of the first glyph (bits):\n")
    w(md_table(pf.jsd.round(3).rename_axis("class")))
    w("")

    w("## 7. Prefixes and suffixes (top 15 by token count, up to 3 glyphs, at least 20 tokens)\n")
    for r in [zl, gc]:
        w(f"### {r['name']}\n")
        cols = ["affix", "n_glyphs", "count", "share", "n_types", "n_free", "ext_share"]
        w("Prefixes\n")
        w(md_table(r["pre"][cols], index=False))
        w("\nSuffixes\n")
        w(md_table(r["suf"][cols], index=False))
        w("")
    w("Affix lists learned by `decompose_table` from ZL3b paragraph text "
      "(ranked by the number of attested stems, dominated affixes dropped):\n")
    w(f"- prefixes: `{' '.join(zl['learned'][0])}`")
    w(f"- suffixes: `{' '.join(zl['learned'][1])}`\n")
    pm = zl["paradigm"]
    w("Paradigm metrics with these lists (stems with >= 5 tokens):\n")
    ser = pd.Series({k: (int(v) if k in ("n_tokens", "n_stems", "n_stems_k") else v) for k, v in pm.items()},
                    name="ZL3b", dtype=object)
    w(md_table(ser.to_frame().rename_axis("metric")))
    w("")

    w("## 8. Comparison with published values\n")
    ext = zl["entropy"].query("preset == 'extended' and scope == 'all'").iloc[0]
    bas = zl["entropy"].query("preset == 'basic' and scope == 'all'").iloc[0]
    rows = [
        {"quantity": "glyph conditional entropy H(X_{t+1} given X_t), composites collapsed, ZL3b",
         "published": "about 2.69 (+-0.1 accepted)", "voytools": f"{ext['H2']:.4f} (extended, scope all); "
         f"{bas['H2']:.4f} (basic, scope all)", "verdict": "agrees" if abs(ext["H2"] - 2.69) <= 0.1 else "DIFFERS"},
        {"quantity": "adjacent-token MI excess share = excess / H(X) (Rozanova & Temerev 2026)", "published": "< 1 %",
         "voytools": f"{zl['mi_tok']['excess_share'] * 100:.2f} % (excess {zl['mi_tok']['excess_bits']:.4f} bits)",
         "verdict": "agrees" if zl["mi_tok"]["excess_share"] < 0.01 else "DIFFERS"},
        {"quantity": "edge-glyph MI", "published": "about 0.2 bits",
         "voytools": f"{zl['mi_edge']['mi_obs']:.4f} observed; {zl['mi_edge']['excess_bits']:.4f} excess over shuffles",
         "verdict": "agrees" if 0.1 <= zl["mi_edge"]["mi_obs"] <= 0.3 else "DIFFERS"},
    ]
    w(md_table(pd.DataFrame(rows), index=False))
    w("")

    w("## 9. Command line interface\n")
    w("All commands run in-process through `voytools.cli.main` (the same code as `python -m voytools`); "
      "`mi` includes the neighbour-similarity test.\n")
    w(md_table(pd.DataFrame(cli_rows), index=False, digits=2))
    w("")
    w("## 10. Timing\n")
    w(f"Parsing and building the tables of ZL3b takes {zl['t_load']:.1f} s; adjacent-token plus edge-glyph MI with "
      f"{SHUFFLES} shuffles takes {zl['t_mi']:.1f} s (GC2a: {gc['t_mi']:.1f} s) and the neighbour-similarity test "
      f"with {SHUFFLES} within-page shuffles {zl['t_nb']:.1f} s (GC2a: {gc['t_nb']:.1f} s), which together make up "
      "the time of the `mi` command in section 9.\n")

    w("## 11. Additional checks\n")
    w("### Token counts of ZL3b under different conventions\n")
    w("The task brief expected roughly 36-38 thousand tokens in total and about 32.7 thousand in paragraph text.  "
      "The counts depend on how the uncertain space `,` and unreadable tokens are treated:\n")
    w(md_table(ex["counts"], index=False))
    n_all, n_ok = ex["reference_count"]
    cj = ex["counts"].set_index(["space policy", "tokens"])
    lig = "{...}"
    diff = zl["tokens_P"] - n_all
    rest = diff - ex["diff_breaks"] - ex["diff_unmatched"]
    w(f"\nThe tokenizer of `voynich-units/analysis/reproduce_headlines.py` (regular expressions; commas and "
      f"periods split; `<...>` tags and `{lig}` ligatures deleted) finds {n_all:,} paragraph tokens in the same "
      f"file, {n_ok:,} of them without any of the characters `?*<>{{}}[]()|`; voytools finds "
      f"{zl['tokens_P']:,} and {zl['tokens_P_clean']:,}.  The difference of {diff} unfiltered tokens consists of "
      f"{ex['diff_breaks']} `<->` / `<~>` drawing breaks inside lines (voytools starts a new token there, as IVTFF 2.0 "
      f"prescribes; the reference merges the two halves), {ex['diff_unmatched']} tokens in loci that the "
      f"reference's pattern does not match (page `fRos`, `Pb` loci) and {rest} further tokens (the remainder; "
      f"ligature-only tokens, which the reference deletes, are the likely cause).  Numbers inside the brief's ranges are obtained with the `join` policy "
      f"(unfiltered {cj.loc[('join', 'unfiltered'), 'all loci']:,} tokens in all loci, "
      f"{cj.loc[('join', 'unfiltered'), 'P']:,} in paragraph text) and, with the default `split` policy, for clean "
      f"tokens ({cj.loc[('split', 'clean'), 'all loci']:,} in all loci, {cj.loc[('split', 'clean'), 'P0']:,} in "
      f"`P0` loci); the unfiltered `split` counts ({cj.loc[('split', 'unfiltered'), 'all loci']:,} and "
      f"{cj.loc[('split', 'unfiltered'), 'P']:,}) lie above them.  The reference script itself splits on both "
      "`.` and `,`, i.e. it uses the same policy as the voytools default.\n")
    w("### Estimator check on a null text\n")
    sm, se = ex["shuffle_mi"], ex["shuffle_edge"]
    w(f"After a within-line shuffle of the ZL3b paragraph text the vocabulary statistics are identical "
      f"({ex['vocab_same']}), the adjacent-token MI excess is {sm['excess_bits']:.4f} bits "
      f"(share {sm['excess_share']:.4f}, z {sm['z']:.2f}) and the edge-glyph MI excess is {se['excess_bits']:.4f} bits "
      f"(z {se['z']:.2f}), as it must be.  The glyph entropy H(X_t+1|X_t) of the space-less text (extended preset, "
      f"scope all) is {ex['real_h2']:.4f} for the real text and {ex['shuffle_h2']:.4f} after the shuffle: it is not a "
      "bag-of-words statistic because glyph pairs at word boundaries depend on word order.\n")
    w("### Currier languages (paragraph text, extended preset, scope line)\n")
    w(md_table(ex["languages"], index=False))
    w("")

    Path(args.out).write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
