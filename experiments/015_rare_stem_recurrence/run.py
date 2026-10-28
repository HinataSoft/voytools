"""Experiment 015: recurrence of rare stems at a distance.

Usage (from the repository root; see ``spec.md``)::

    python experiments/015_rare_stem_recurrence/run.py run --combo ZL_split                 # all sources
    python experiments/015_rare_stem_recurrence/run.py run --combo GC_join --sources real,perm_tokens
    python experiments/015_rare_stem_recurrence/run.py run --combo ctl:pliny_prose             # a control (009)
    python experiments/015_rare_stem_recurrence/run.py run --combo ctl:gibberish_pooled         # 013
    python experiments/015_rare_stem_recurrence/run.py summarize

Combinations: ``{ZL, GC}_{split, join}`` (real text, ``P`` clean) and ``ctl:<name>`` for the controls of
experiment 009 (``pliny_prose``, ``celsus_prose``, ``linnaeus_items`` ... in the ZL split layout) and 013
(``gibberish_pooled``, ``voynich_crop``: 60 evenly spaced pages, so distances only up to ~8 pages of the
cropped order).  Sources (``seed`` = repetition):

* ``real`` (seed 0); ``real_boot`` - bootstrap over stems (stems resampled with replacement, seeds 0..199);
* ``perm_tokens`` - the spec's null: every occurrence of a rare stem moved to a random page with
  probability proportional to page size (``recurrence.permute_pages``), 200 seeds; the exact expectation
  of the same null is computed for every sample (columns ``e_*``, ``*_exp``), the sampled version gives
  the interval and checks the expectation;
* ``perm_global`` / ``perm_section`` - whole pages permuted among page positions globally / inside the
  section (006), 100 seeds each; metadata (rank, section, language, hand, quire) stay with the position,
  the page size travels with the content and the exact expectations of the sample use the permuted sizes;
* ``selfcitation_21`` - ts2020 config 21 (004), seeds 2100 + rep, 30 reps, Eva only;
* ``template_005`` / ``template_005_notopic`` - the 005 template (configuration ``best`` of the file) with
  its topic and with ``topic_concentration = inf`` (variant ``no_topic``), seeds 100000 + rep, pool_seed
  = rep, 30 reps each;
* ``autocopy_86`` - auto-copy generator config 86 of 003 (``edit_model='trigram'``), seeds 8600 + rep, 30.

Every run reports the metrics of ``voytools.recurrence.recurrence_metrics`` for the decompositions
(a) learned affixes and (b) learned + positional glyphs (011) and for the bands 2-5 (``r``), 6-20 (``c``)
and label stems with 2-5 text occurrences (``l``): column ``<var>__<band>__<metric>``.  Controls use
their own learned affixes, decomposition (a) only, no label band.  Generated samples are decomposed
with the lists of the real text; their label band uses the real file's label stems.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_runs, load_tables, run_seeds, seed_range  # noqa: E402
from voytools import recurrence as R  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.nulls import autocopy_like  # noqa: E402
from voytools.selfcitation import selfcitation_like  # noqa: E402
from voytools.templates import ParadigmTemplate  # noqa: E402
from voytools.wordgrammar import learn_word_grammar  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M6 = _load("exp006", "006_topic_at_distance")

DATA = ROOT / "data"
OUT = ROOT / "results" / "015_rare_stem_recurrence"
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
COMBOS = ["ZL_split", "ZL_join", "GC_split", "GC_join"]
CONTROLS_009 = ("pliny_prose", "celsus_prose", "linnaeus_items")
CONTROLS_013 = ("gibberish_pooled", "voynich_crop")
BANDS = {"r": (2, 5), "c": (6, 20)}
LABEL_BAND = "l"
POS_AFFIXES = {   # decomposition (b): positional glyphs added to the learned lists (008 / 011)
    "ZL": {"prefixes": ["p", "f", "t", "k", "cph", "cfh", "cth", "ckh"], "suffixes": ["m", "g", "am", "om"]},
    "GC": {"prefixes": ["g", "j", "f", "u", "k", "h", "G", "J", "F", "U", "K", "H"],
           "suffixes": ["p", "P", "*", "ap", "op"]},
}
SOURCES = {"real": range(1), "real_boot": range(200), "perm_tokens": range(200), "perm_global": range(100),
           "perm_section": range(100), "selfcitation_21": range(30), "template_005": range(30),
           "template_005_notopic": range(30), "autocopy_86": range(30)}
CONTROL_SOURCES = ("real", "real_boot", "perm_tokens")
TS_SEED0, TPL_SEED0, AC_SEED0 = 2100, 100_000, 8600


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


class Context:
    """One corpus (real combination or control): token table, page table, affix lists, label stems."""

    def __init__(self, combo: str, with_grammar: bool = False):
        self.combo = combo
        self.is_control = combo.startswith("ctl:")
        if self.is_control:
            name = combo[4:]
            self.key, self.policy = "ZL", "split"
            if name in CONTROLS_013:
                M13 = _load("exp013", "013_gibberish_and_cipher_controls")
                self.real = M13.build(name).reset_index(drop=True)
            else:
                M9 = _load("exp009", "009_register_profile")
                self.real = M9.control_table(name).reset_index(drop=True)
            self.file = name
            self.rank = R.page_order(self.real)
            self.real_d = decompose_table(self.real)
            self.labels = None
        else:
            self.key, self.policy = combo.split("_")
            self.file = FILES[self.key]
            tables = load_tables(DATA / self.file, space_policy=self.policy)
            self.rank = R.page_order(tables.tokens)
            self.real = filter_tokens(tables.tokens, generic_type="P", clean=True).reset_index(drop=True)
            self.real_d = decompose_table(self.real)
            self.labels = filter_tokens(tables.tokens, generic_type="L", clean=True).reset_index(drop=True)
        self.prefixes, self.suffixes = self.real_d.attrs["prefixes"], self.real_d.attrs["suffixes"]
        self.variants = {"a": (list(self.prefixes), list(self.suffixes))}
        if not self.is_control:
            add = POS_AFFIXES[self.key]
            self.variants["b"] = (self.prefixes + [a for a in add["prefixes"] if a not in self.prefixes],
                                  self.suffixes + [a for a in add["suffixes"] if a not in self.suffixes])
        self.pages = R.page_table(self.real, self.rank)
        self.expect_same = R.same_by_distance(self.pages)
        self.expect_profile = R.expected_profile(self.pages)
        self.label_stems = {}
        for var, (pre, suf) in self.variants.items():
            if self.labels is not None:
                self.label_stems[var] = R.label_stems(decompose_table(self.labels, prefixes=pre, suffixes=suf))
        self.grammar = None
        if with_grammar:
            self.grammar = learn_word_grammar(self.real["token"].tolist(), self.real.attrs.get("glyph_preset"))
        self._tpl = {}

    def columns(self) -> dict:
        return {"combo": self.combo, "file": self.file, "space_policy": self.policy, "n_pages": len(self.pages)}

    def decompositions(self, df: pd.DataFrame) -> dict:
        return {var: decompose_table(df, prefixes=pre, suffixes=suf) for var, (pre, suf) in self.variants.items()}

    def occurrences(self, df: pd.DataFrame) -> dict:
        """``{(variant, band): occurrence table}`` of a sample laid out on ``self.pages``."""
        out = {}
        for var, dfd in self.decompositions(df).items():
            for band, (lo, hi) in BANDS.items():
                out[(var, band)] = R.stem_occurrences(dfd, self.pages, lo, hi)
            if var in self.label_stems:
                out[(var, LABEL_BAND)] = R.stem_occurrences(dfd, self.pages, *BANDS["r"], within=self.label_stems[var])
        return out

    def template(self, notopic: bool) -> ParadigmTemplate:
        if notopic not in self._tpl:
            cfg = M6.template_config(self.key)
            if notopic:
                cfg["topic_concentration"] = math.inf
            self._tpl[notopic] = ParadigmTemplate.learn(self.real_d, **cfg)
        return self._tpl[notopic]


def metrics(ctx: Context, occs: dict, pages: pd.DataFrame | None = None) -> dict:
    """Metrics of all (variant, band) occurrence tables; ``pages`` overrides the page table (page permutations)."""
    if pages is None:
        pages, es, ep = ctx.pages, ctx.expect_same, ctx.expect_profile
    else:
        es, ep = R.same_by_distance(pages), R.expected_profile(pages)
    out = {}
    for (var, band), occ in occs.items():
        for k, v in R.recurrence_metrics(occ, pages, es, ep).items():
            out[f"{var}__{band}__{k}"] = round(v, 6) if isinstance(v, float) and np.isfinite(v) else v
    return out


def bootstrap_stems(occ: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Stems resampled with replacement; a stem drawn twice becomes two distinct stems (``stem#k``)."""
    rng = np.random.default_rng(seed)
    stems = pd.unique(occ["stem"].to_numpy())
    draw = rng.choice(stems, size=len(stems), replace=True)
    parts = []
    for k, (s, cnt) in enumerate(pd.Series(draw).value_counts().items()):
        block = occ[occ["stem"] == s]
        for r in range(int(cnt)):
            parts.append(block.assign(stem=f"{s}#{r}") if r else block)
    out = pd.concat(parts, ignore_index=True) if parts else occ.iloc[:0]
    return out.sort_values(["stem", "folio"], kind="stable").reset_index(drop=True)


def page_permutation(ctx: Context, within: str | None, seed: int) -> np.ndarray:
    """Content of page ``p`` goes to the position of page ``perm[p]`` (inside the groups ``within``)."""
    rng = np.random.default_rng(seed)
    folios = ctx.pages.index.to_numpy()
    n = len(folios)
    if within is None:
        return folios[rng.permutation(n)]
    groups = ctx.pages[within].to_numpy()
    perm = np.arange(n)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        perm[idx] = idx[rng.permutation(len(idx))]
    return folios[perm]


def make_runner(ctx: Context, source: str, real_occs: dict, ac_params: dict | None):
    def one(seed: int) -> dict:
        base = {**ctx.columns(), "source": source}
        if source == "real":
            return {**base, **metrics(ctx, real_occs)}
        if source == "real_boot":
            return {**base, **metrics(ctx, {k: bootstrap_stems(o, seed) for k, o in real_occs.items()})}
        if source == "perm_tokens":
            return {**base, **metrics(ctx, {k: R.permute_pages(o, ctx.pages, seed) for k, o in real_occs.items()})}
        if source.startswith("perm_"):
            within = {"perm_global": None, "perm_section": "section"}[source]
            target = page_permutation(ctx, within, seed)
            move = dict(zip(ctx.pages.index.to_numpy(), target))
            # page sizes travel with the content: the expectations of this sample use the permuted sizes
            pages = ctx.pages.copy()
            pages.loc[target, "size"] = ctx.pages["size"].to_numpy()
            pages["weight"] = pages["size"] / pages["size"].sum()
            occs = {k: R.relocate(o, pages, o["folio"].map(move).to_numpy()) for k, o in real_occs.items()}
            return {**base, **metrics(ctx, occs, pages)}
        if source == "selfcitation_21":
            if ctx.key != "ZL" or ctx.is_control:
                raise ValueError("selfcitation_like generates basic Eva; no null for a non-Eva file")
            sample = selfcitation_like(ctx.real, seed=TS_SEED0 + seed, **M6.selfcitation_params())
        elif source in ("template_005", "template_005_notopic"):
            sample = ctx.template(source.endswith("notopic")).sample(ctx.real, seed=TPL_SEED0 + seed, pool_seed=seed)
        elif source == "autocopy_86":
            sample = autocopy_like(ctx.real, allowed_from=ctx.real["token"].tolist(), word_grammar=ctx.grammar,
                                   seed=AC_SEED0 + seed, **ac_params)
        else:
            raise ValueError(source)
        return {**base, **metrics(ctx, ctx.occurrences(sample))}

    return one


def progress(i: int, n: int, s: int) -> None:
    if i % 25 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_run(args) -> None:
    default = CONTROL_SOURCES if args.combo.startswith("ctl:") else list(SOURCES)
    sources = list(default) if args.sources == "all" else args.sources.split(",")
    ctx = Context(args.combo, with_grammar="autocopy_86" in sources)
    real_occs = ctx.occurrences(ctx.real)
    ac = M6.autocopy_params() if "autocopy_86" in sources else None
    parts = OUT / "runs" / args.combo.replace(":", "_")
    parts.mkdir(parents=True, exist_ok=True)
    ctx.pages.reset_index().to_csv(parts / "pages.csv", index=False)
    stems = pd.concat([o.groupby("stem", sort=False).agg(count=("count", "first"), n_pages=("folio", "nunique"),
                                                          n_quires=("quire", "nunique"), n_forms=("token", "nunique"))
                       .assign(variant=v, band=b) for (v, b), o in real_occs.items()]).reset_index()
    stems.to_csv(parts / "stems.csv", index=False)
    for source in sources:
        seeds = M6.parse_seeds(args.seeds) if args.seeds else list(SOURCES[source])
        print(f"{args.combo} {source}: {len(seeds)} seeds", flush=True)
        params = {"combo": args.combo, "source": source, "bands": BANDS, "variants": ctx.variants,
                  "n_label_stems": {v: len(s) for v, s in ctx.label_stems.items()}, "n_pages": len(ctx.pages)}
        if source == "autocopy_86":
            params["autocopy"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in ac.items()}
            params["seed_offset"] = AC_SEED0
        if source == "selfcitation_21":
            params["selfcitation_004"] = {"row": M6.TS_ROW, **{k: (list(v) if isinstance(v, tuple) else v)
                                                               for k, v in M6.selfcitation_params().items()}}
            params["seed_offset"] = TS_SEED0
        if source.startswith("template_005"):
            cfg = M6.template_config(ctx.key)
            if source.endswith("notopic"):
                cfg["topic_concentration"] = "inf"
            params["template_005"] = cfg
            params["seed_offset"], params["pool_seed"] = TPL_SEED0, "repetition"
        if ctx.is_control:
            params["control"] = {k: v for k, v in ctx.real.attrs.get("control", {}).items()}
        run_seeds(f"015_{args.combo}_{source}", seeds, make_runner(ctx, source, real_occs, ac),
                  parts / source, params=params, progress=progress)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--combo", required=True)
    p.add_argument("--sources", default="all")
    p.add_argument("--seeds", default=None, help="'0:30', '0:200/4:1' or '3,5'")
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("summarize")
    p.set_defaults(func=lambda a: __import__("summarize").summarize())
    args = ap.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    args.func(args)


if __name__ == "__main__":
    main()
