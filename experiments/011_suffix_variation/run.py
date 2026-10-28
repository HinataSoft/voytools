"""Experiment 011: decomposition of within-word variation - is the "paradigm" a positional effect?

How much of H(suffix | stem) (and H(prefix | stem)) is explained by position in the line, Currier
language, scribe, section and page, and does the rest have structure (agreement with the
neighbour, consistency inside the page/paragraph, section after position + language + scribe)?

Measured with the same code (``voytools.variation``) on

* the Voynich paragraph text (``ZL_split``, ``ZL_join``, ``GC_split``, ``GC_join``), decomposition
  (a) learned affixes and (b) learned affixes + positional glyphs (experiment 008, variant b);
* real-language controls of experiment 009 (``voytools.controls`` via the loaders of
  ``experiments/009_register_profile/run.py``), laid out in the ZL split layout; decomposition (a);
* nulls: ts2020 calB, structured template 005, auto-copy 86 (the generators of 009, seeds 0..29).

Usage (from the repository root)::

    python experiments/011_suffix_variation/run.py measure --source ZL_split          # point + 100 subsamples
    python experiments/011_suffix_variation/run.py measure --source ctl:celsus_prose  # point + 30 subsamples
    python experiments/011_suffix_variation/run.py measure --source ts2020:ZL_split   # 30 seeds
    python experiments/011_suffix_variation/run.py top                                # top_stems.csv
    python experiments/011_suffix_variation/run.py summarize

Seeds: real text and controls ``-1`` = full text (point estimate, 100 permutations per correction),
``0..99`` (controls ``0..29``) = subsamples of 80 % of the pages (20 permutations); nulls ``0..29``
(20 permutations).  Every command is resumable; runs go to ``results/011_suffix_variation/runs/<source>/``.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import run_seeds  # noqa: E402
from voytools import variation as V  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M9 = _load("exp009", "009_register_profile")

OUT = ROOT / "results" / "011_suffix_variation"
REAL = ("ZL_split", "ZL_join", "GC_split", "GC_join")
CONTROLS = ("caesar_prose", "celsus_prose", "pliny_prose", "culpeper_prose", "linnaeus_prose", "linnaeus_items")
NULLS = {"ts2020": ("ZL_split", "ZL_join"), "tmpl": ("ZL_split", "GC_split"), "ac86": ("ZL_split", "GC_split")}
N_SUB = 100
N_SUB_CTL = 30
N_NULL = 30
N_PERM = {"point": 100, "sub": 20, "null": 20}
MIN_STEM = 20           # main: stems with >= 20 tokens
MIN_STEM_CHECK = 5      # check (point estimate only)
# (b) positional glyphs added to the learned affix lists (008, POS_AFFIXES)
POS_AFFIXES = {
    "ZL": {"prefixes": ["p", "f", "t", "k", "cph", "cfh", "cth", "ckh"], "suffixes": ["m", "g", "am", "om"]},
    "GC": {"prefixes": ["g", "j", "f", "u", "k", "h", "G", "J", "F", "U", "K", "H"],
           "suffixes": ["p", "P", "*", "ap", "op"]},
}
ALL_NO_PAGE = ["cv_pos", "cv_relpos", "cv_place", "cv_lip", "cv_lang", "cv_hand", "cv_sec"]
BY = {
    "pos": "cv_pos", "relpos": "cv_relpos", "place": "cv_place", "lip": "cv_lip",
    "lang": "cv_lang", "hand": "cv_hand", "sec": "cv_sec", "page": "cv_page",
    "pos+lang": ["cv_pos", "cv_lang"], "pos+hand": ["cv_pos", "cv_hand"],
    "pos+lang+sec": ["cv_pos", "cv_lang", "cv_sec"], "all-page": ALL_NO_PAGE, "all": ALL_NO_PAGE + ["cv_page"],
}
TARGETS = ("suffix", "prefix")
FIELDS = ("H_given", "H_cond", "bias", "reduction", "explained", "explained_raw", "z")
COND = ["cv_pos", "cv_lang", "cv_hand"]   # residual structure: after position + language + scribe
RESID_BY = {"sec": "cv_sec", "page": "cv_page"}


# --------------------------------------------------------------------------- #
# sources
# --------------------------------------------------------------------------- #


def kind_of(source: str) -> str:
    if source.startswith("ctl:"):
        return "control"
    return source.split(":")[0] if ":" in source else "real"


def variants(ctx) -> dict:
    """Affix lists of decomposition (a) and, for transliterations, (b) with positional glyphs."""
    out = {"a": (list(ctx.prefixes), list(ctx.suffixes))}
    if not ctx.source.startswith("ctl:"):
        add = POS_AFFIXES[ctx.base_tag.split("_")[0]]
        out["b"] = (list(ctx.prefixes) + [a for a in add["prefixes"] if a not in ctx.prefixes],
                    list(ctx.suffixes) + [a for a in add["suffixes"] if a not in ctx.suffixes])
    return out


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #


def measure(df: pd.DataFrame, ctx, seed: int, level: str) -> dict:
    """All 011 metrics of one corpus, for each decomposition variant (column prefix ``a__`` / ``b__``)."""
    k = N_PERM[level]
    s = max(seed, 0)
    df = V.add_covariates(df.reset_index(drop=True))
    out = {"n_tokens": len(df), "n_pages": int(df["cv_page"].nunique())}
    for var, (pre, suf) in variants(ctx).items():
        dfd = decompose_table(df, prefixes=pre, suffixes=suf)
        out[f"{var}__share_empty_stem"] = float((dfd["stem"] == "").mean())
        for target in TARGETS:
            dec = V.conditional_entropy_decomposition(dfd, target, "stem", BY, n_perm=k, min_count=MIN_STEM, seed=s)
            cv = V.crossval_explained(dfd, target, "stem", BY, min_count=MIN_STEM, seed=s).set_index("X")
            out[f"{var}__{target}__n_stems"] = int(dec["n_stems"].iloc[0])
            out[f"{var}__{target}__n_tok"] = int(dec["n_tokens"].iloc[0])
            for _, r in dec.iterrows():
                for f in FIELDS:
                    out[f"{var}__{target}__{r['X']}__{f}"] = float(r[f])
                out[f"{var}__{target}__{r['X']}__explained_cv"] = float(cv.loc[r["X"], "explained_cv"])
            if level == "point":
                d5 = V.conditional_entropy_decomposition(dfd, target, "stem", BY, n_perm=k,
                                                         min_count=MIN_STEM_CHECK, seed=s)
                out[f"{var}__{target}__n_stems_min5"] = int(d5["n_stems"].iloc[0])
                for _, r in d5.iterrows():
                    out[f"{var}__{target}__{r['X']}__explained_min5"] = float(r["explained"])
                    out[f"{var}__{target}__{r['X']}__H_given_min5"] = float(r["H_given"])
        # residual structure of the suffix choice, after position class + language + scribe
        given = ["stem"] + COND
        for name, nb, off in (("nb_prev_suffix", "suffix", -1), ("nb_next_prefix", "prefix", 1)):
            r = V.neighbour_conditional_mi(dfd, "suffix", nb, off, given, unit="cv_line", n_perm=k,
                                           min_count=MIN_STEM, seed=s)
            for f in ("mi_raw", "mi_perm", "mi", "perm_sd", "p", "n_pairs"):
                out[f"{var}__{name}__{f}"] = r[f]
        for unit in ("page", "par"):
            r = V.mode_consistency(dfd, "suffix", "stem", f"cv_{unit}", groups=COND, min_count=MIN_STEM,
                                   n_perm=k, seed=s)
            for f in ("within", "expected", "excess", "perm_sd", "z", "p", "global", "n_cells", "n_tokens"):
                out[f"{var}__cons_{unit}__{f}"] = r[f]
        dec = V.conditional_entropy_decomposition(dfd, "suffix", given, RESID_BY, n_perm=k, min_count=MIN_STEM,
                                                  seed=s)
        for _, r in dec.iterrows():
            for f in ("H_given", "reduction_raw", "bias", "reduction", "explained", "z"):
                out[f"{var}__resid_{r['X']}__{f}"] = float(r[f])
    return out


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def safe(source: str) -> str:
    return source.replace(":", "__")


def progress(i: int, n: int, s: int) -> None:
    if i % 5 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_measure(args) -> None:
    source = args.source
    kind = kind_of(source)
    ctx = M9.Ctx(source if kind in ("real", "control") else f"real:{source.split(':')[1]}")
    if kind in ("real", "control"):
        def one(seed: int) -> dict:
            df = ctx.real if seed < 0 else M9.page_subsample(ctx.real, seed)
            level = "point" if seed < 0 else "sub"
            return {"source": source, "kind": level, **measure(df, ctx, seed, level)}
        default = [-1] + list(range(N_SUB if kind == "real" else N_SUB_CTL))
    else:
        def one(seed: int) -> dict:
            return {"source": source, "kind": "null", **measure(ctx.corpus(kind, seed), ctx, seed, "null")}
        default = list(range(N_NULL))
    seeds = M9.parse_seeds(args.seeds, default)
    params = {"source": source, "variants": {k: {"prefixes": v[0], "suffixes": v[1]} for k, v in variants(ctx).items()},
              "by": BY, "cond_residual": COND, "n_perm": N_PERM, "min_stem": MIN_STEM,
              "min_stem_check": MIN_STEM_CHECK, "sub_frac": M9.SUB_FRAC, "n_sub": N_SUB, "n_sub_ctl": N_SUB_CTL,
              "n_null": N_NULL, "ts2020_config": M9.TS2020_CONFIG, "ac_seed0": M9.AC_SEED0,
              "tmpl_seed0": M9.TMPL_SEED0, "control": ctx.real.attrs.get("control")}
    run_seeds("011_suffix_variation", seeds, one, OUT / "runs" / safe(source), params=params, progress=progress)


def cmd_top(args) -> None:
    """Descriptive: p(affix | stem, place in line / line of paragraph) of the 10 most frequent stems."""
    rows = []
    for tag in ("ZL_split", "GC_split"):
        ctx = M9.Ctx(tag)
        df = V.add_covariates(ctx.real)
        for var, (pre, suf) in variants(ctx).items():
            dfd = decompose_table(df, prefixes=pre, suffixes=suf)
            dfd = dfd[dfd["stem"] != ""]
            for target in TARGETS:
                for cov in ("cv_place", "cv_lip"):
                    t = V.affix_by_covariate(dfd, target, "stem", cov, top=10)
                    t = t.rename(columns={cov: "position", target: "affix"})
                    t.insert(0, "covariate", cov[3:])
                    t.insert(0, "target", target)
                    t.insert(0, "variant", var)
                    t.insert(0, "source", tag)
                    rows.append(t)
    OUT.mkdir(parents=True, exist_ok=True)
    out = pd.concat(rows, ignore_index=True)
    out["position"] = out["position"].astype(str)
    out.round(4).to_csv(OUT / "top_stems.csv", index=False)
    print(f"top_stems.csv: {len(out)} rows")


def all_sources() -> list[str]:
    return ([*REAL] + [f"ctl:{c}" for c in CONTROLS]
            + [f"{n}:{t}" for n, tags in NULLS.items() for t in tags])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("measure")
    p.add_argument("--source", required=True)
    p.add_argument("--seeds", default=None)
    sub.add_parser("top")
    sub.add_parser("summarize")
    sub.add_parser("sources")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from summarize import summarize  # noqa: E402  (experiments/011_suffix_variation/summarize.py)

        summarize()
        return
    if args.cmd == "sources":
        print("\n".join(all_sources()))
        return
    {"measure": cmd_measure, "top": cmd_top}[args.cmd](args)


if __name__ == "__main__":
    main()
