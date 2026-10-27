"""Experiment 005: structured-template null ("artificial language") and robustness of experiment 002.

Usage (from the repository root; see ``spec.md``)::

    # step 1, phase A: 216 configurations x 3 reps on one text, in 4 parallel parts
    python experiments/005_structured_template_null/run.py grid --tag ZL_split --part 0 --parts 4
    python experiments/005_structured_template_null/run.py rank-grid          # grid_summary.csv
    # step 1, phase B: 5 best configurations x 50 seeds (+ the best one on the join policy)
    python experiments/005_structured_template_null/run.py phase-b --tag ZL_split --part 0 --parts 2
    # step 2: ablations / additions of the best configuration, 50 seeds each
    python experiments/005_structured_template_null/run.py ablation --tag ZL_split --variant no_topic
    # step 3: robustness of 002 (real 80 % page subsamples, slot null, autocopy_536)
    python experiments/005_structured_template_null/run.py robust --file ZL --source real_sub --seeds 0:100
    # step 4: real A/B x scribe contrasts with page-permutation baselines
    python experiments/005_structured_template_null/run.py ab-real
    # tables, plot, checks
    python experiments/005_structured_template_null/run.py summarize

Texts ("tags"): ``ZL_split``, ``GC_split`` (main), ``ZL_join``, ``GC_join`` (control), paragraph text,
clean tokens, learned affix lists (``decompose_table`` defaults), exactly as in 002 (its ``Context`` is
imported).  Every template sample is decomposed with the real lists of its text and goes through the
metric function of 002 (copied unchanged as :func:`metrics_decomposed`; ``check-002`` verifies that the
copy reproduces 002's ``runs.csv`` for the real text), plus vocabulary targets and the A/B x scribe
suffix contrasts of step 4.

Run ids: phase A ``seed = config_id * 100 + rep`` (rep 0..2), phase B and ablations
``seed = config_id * 100 + rep`` (rep 0..49).  The artificial stem inventory of a sample is chosen by
``pool_seed = rep`` (common random numbers: configurations with the same rep share their stems).
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import load_runs, run_seeds  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402
from voytools.morphology import (  # noqa: E402
    apply_bpe,
    decompose_table,
    inventory_jaccard,
    learn_bpe,
    paradigm_metrics,
    productivity_slope,
    slot_independent_sample,
    stem_dispersion,
)
from voytools.nulls import autocopy_like  # noqa: E402
from voytools.stats import jsd, vocabulary_summary  # noqa: E402
from voytools.templates import ParadigmTemplate  # noqa: E402

OUT = ROOT / "results" / "005_structured_template_null"
R002 = ROOT / "results" / "002_paradigm_coherence"

_spec = importlib.util.spec_from_file_location("exp002", ROOT / "experiments" / "002_paradigm_coherence" / "run.py")
E2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E2)

MIN_STEM, MIN_INV = E2.MIN_STEM, E2.MIN_INV
TAGS = ("ZL_split", "GC_split", "ZL_join", "GC_join")
MAIN_TAGS = ("ZL_split", "GC_split")

# --------------------------------------------------------------------------- #
# grid (step 1)
# --------------------------------------------------------------------------- #

AXES = {
    "n_stems": (2000, 4000, 8000),
    "zipf_s": (0.9, 1.1, 1.3),
    "paradigm": ("fixed",),
    "prefix_dependence": ("none", "stem"),
    "topic_concentration": (math.inf, 5.0, 1.0, 0.2),
    "page_locality": (0.0, 0.3, 0.6),
}
CONFIGS = [dict(zip(AXES, v)) for v in itertools.product(*AXES.values())]
REPS_A, REPS_B, TOP_B = 3, 50, 5
# loss as in 001: sum of w * ((x - t) / t) ** 2 over the target metrics
TARGETS = {"n_types": 1.0, "hapax_share": 1.0, "mean_len": 1.0, "n_stems_k": 2.0, "mean_suffixes_per_stem": 2.0,
           "frac_stems_multi_suffix": 2.0}
LOSS_LIMIT = 0.15
# ablations (remove one ingredient) and additions (one more hand-made rule) of the best configuration;
# ids 1000+ so that they never collide with grid configurations
VARIANTS = {
    "best": {},
    "no_topic": {"topic_concentration": math.inf},
    "no_locality": {"page_locality": 0.0},
    "no_prefix_link": {"prefix_dependence": "none"},
    "free_paradigm": {"paradigm": "free"},
    "add_prefix_given_suffix": {"prefix_dependence": "suffix"},
    "add_prefix_stem_suffix": {"prefix_dependence": "stem_suffix"},
    "add_suffix_by_language": {"suffix_by": "language"},
    "add_all": {"prefix_dependence": "stem_suffix", "suffix_by": "language"},
}
VARIANT_ID = {v: 1000 + i for i, v in enumerate(VARIANTS)}

# --------------------------------------------------------------------------- #
# robustness (step 3)
# --------------------------------------------------------------------------- #

ROBUST = {  # variant -> decomposition settings (learned on the full real text, split policy)
    "base": {},
    "suf8": {"n_suffixes": 8},
    "suf20": {"n_suffixes": 20},
    "len3": {"max_len": 3},
    "suffix_first": {"order": "suffix_first"},
    "bpe64": {"bpe": 64},
}
ROBUST_METRICS = ["cmi_excess", "mi_prefix_suffix", "mean_suffixes_per_stem", "frac_stems_multi_suffix",
                  "frac_multi_suffix_5to9", "H_suffix_given_stem", "n_stems_k", "jaccard_suffix_sections",
                  "jaccard_stem_sections", "jaccard_gap_sections", "jsd_suffix_sections", "jsd_stem_sections",
                  "jsd_suffix_AB", "stem_dispersion"]
ROBUST_SOURCES = {"real": range(1), "real_sub": range(200), "slot_indep": range(50), "autocopy_536": range(50)}
SUB_FRAC = 0.8

# --------------------------------------------------------------------------- #
# A/B x scribe contrasts (step 4)
# --------------------------------------------------------------------------- #

# (name, kind, group a, group b); groups are language+hand labels such as 'A1', 'B3'
AB_PAIRS = [
    ("AB_all", "language", "A", "B"),
    ("A1_vs_B2", "between_language", "A1", "B2"),
    ("A1_vs_B3", "between_language", "A1", "B3"),
    ("B2_vs_B3", "within_language_B", "B2", "B3"),
    ("B2_vs_B5", "within_language_B", "B2", "B5"),
    ("B3_vs_B5", "within_language_B", "B3", "B5"),
    ("H3_A_vs_B", "within_hand_3", "A3", "B3"),
]
SECTION_HANDS = ("A1", "B2", "B3")  # within-scribe, between-section contrasts
MIN_SECTION_TOKENS = 500
N_PERM = 200


def parse_seeds(spec: str) -> list[int]:
    """``'0:20'`` -> seeds 0..19; ``'3,5,8'`` -> those seeds."""
    if ":" in spec:
        a, b = spec.split(":")
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(",") if s.strip()]


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #


def metrics_decomposed(dfd: pd.DataFrame, seed: int) -> dict:
    """The metrics of 002 (``E2.metrics``) for an already decomposed table; body copied unchanged."""
    pm = paradigm_metrics(dfd, min_stem_count=MIN_STEM, seed=seed)
    s = pm.summary
    low = pm.per_stem[pm.per_stem["count"] < 2 * MIN_STEM]
    out = {k: float(s[k]) for k in ("n_tokens", "n_stems", "n_stems_k", "token_coverage_k", "mean_suffixes_per_stem",
                                     "frac_stems_multi_suffix", "H_suffix_given_stem", "mi_prefix_suffix",
                                     "cmi_prefix_suffix_given_stem", "cmi_shuffled", "cmi_excess")}
    # frequency-controlled: stems with 5..9 tokens only
    out["frac_multi_suffix_5to9"] = float((low["n_suffixes"] > 1).mean()) if len(low) else float("nan")
    out["jaccard_suffix_sections"], out["jsd_suffix_sections"] = E2._mean_pairs(dfd, "illustration", "suffix")
    out["jaccard_stem_sections"], out["jsd_stem_sections"] = E2._mean_pairs(dfd, "illustration", "stem")
    out["jaccard_gap_sections"] = out["jaccard_suffix_sections"] - out["jaccard_stem_sections"]
    out["jaccard_suffix_AB"], out["jsd_suffix_AB"] = E2._mean_pairs(dfd, "language", "suffix", ["A", "B"])
    out["jaccard_stem_AB"], out["jsd_stem_AB"] = E2._mean_pairs(dfd, "language", "stem", ["A", "B"])
    out["jaccard_gap_AB"] = out["jaccard_suffix_AB"] - out["jaccard_stem_AB"]
    ps = productivity_slope(dfd, min_stem_count=MIN_STEM)
    out["productivity_slope"], out["productivity_intercept"], out["productivity_r2"] = ps["slope"], ps["intercept"], ps["r2"]
    out["stem_dispersion"] = stem_dispersion(dfd, min_stem_count=MIN_STEM)
    return out


def lang_hand(df: pd.DataFrame) -> np.ndarray:
    return (df["language"].fillna("?").astype(str) + df["hand"].fillna("?").astype(str)).to_numpy()


def _jsd_of(values: np.ndarray, mask_a: np.ndarray, mask_b: np.ndarray, support: list[str]) -> float:
    ca = pd.Series(values[mask_a]).value_counts()
    cb = pd.Series(values[mask_b]).value_counts()
    return jsd([int(ca.get(v, 0)) for v in support], [int(cb.get(v, 0)) for v in support])


def ab_contrasts(dfd: pd.DataFrame, column: str = "suffix") -> dict:
    """JSD of the ``column`` distribution for the A/B x scribe contrasts (step 4)."""
    lh = lang_hand(dfd)
    lang = dfd["language"].fillna("?").astype(str).to_numpy()
    sec = dfd["illustration"].fillna("?").astype(str).to_numpy()
    vals = dfd[column].astype(str).to_numpy()
    support = sorted(set(vals.tolist()))
    out = {}
    for name, kind, a, b in AB_PAIRS:
        key = lang if kind == "language" else lh
        out[f"ab_{column}_{name}"] = _jsd_of(vals, key == a, key == b, support)
    for h in SECTION_HANDS:
        m = lh == h
        secs = [s for s, c in pd.Series(sec[m]).value_counts().items() if c >= MIN_SECTION_TOKENS]
        pairs = [_jsd_of(vals, m & (sec == x), m & (sec == y), support) for i, x in enumerate(secs) for y in secs[i + 1:]]
        out[f"ab_{column}_sections_within_{h}"] = float(np.mean(pairs)) if pairs else float("nan")
    return out


def vocab(df: pd.DataFrame) -> dict:
    v = vocabulary_summary(df).iloc[0]
    return {"n_types": float(v["n_types"]), "hapax_share": float(v["hapax_share"]), "mean_len": float(v["mean_len"])}


def all_metrics(df: pd.DataFrame, ctx, seed: int) -> dict:
    """Vocabulary targets + 002 metrics + A/B x scribe contrasts, real affix lists of ``ctx``."""
    dfd = decompose_table(df, prefixes=ctx.prefixes, suffixes=ctx.suffixes)
    return {**vocab(df), **metrics_decomposed(dfd, seed), **ab_contrasts(dfd)}


# --------------------------------------------------------------------------- #
# contexts and templates
# --------------------------------------------------------------------------- #

_CTX: dict = {}


def context(tag: str):
    if tag not in _CTX:
        _CTX[tag] = E2.Context(f"{tag}_learned")
    return _CTX[tag]


_TPL: dict = {}


def template(tag: str) -> ParadigmTemplate:
    if tag not in _TPL:
        _TPL[tag] = ParadigmTemplate.learn(context(tag).real_d)
    return _TPL[tag]


def config_columns(cfg: dict) -> dict:
    return {k: (str(v) if k == "topic_concentration" else v) for k, v in cfg.items()}


def template_run(tag: str, cfg: dict, config_id: int, rep: int, seed: int) -> dict:
    ctx = context(tag)
    tpl = template(tag).with_params(**cfg)
    x = tpl.sample(ctx.real, seed=seed, pool_seed=rep)
    info = x.attrs["template_info"]
    return {"tag": tag, "config_id": config_id, "rep": rep, **config_columns(cfg),
            "n_stems_used": info["n_stems_used"], "share_stems_real": info["share_stems_real"],
            "share_tokens_local": info["share_tokens_local"], **all_metrics(x, ctx, seed)}


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def _warm(tag: str, reps) -> None:
    for r in reps:
        template(tag).stem_pool(r, max(AXES["n_stems"]))


# --------------------------------------------------------------------------- #
# step 1
# --------------------------------------------------------------------------- #


def cmd_grid(args) -> None:
    ids = [i for i in range(len(CONFIGS)) if i % args.parts == args.part]
    seeds = [i * 100 + r for i in ids for r in range(REPS_A)]
    _warm(args.tag, range(REPS_A))

    def one_run(seed: int) -> dict:
        cid, rep = divmod(seed, 100)
        return template_run(args.tag, CONFIGS[cid], cid, rep, seed)

    run_seeds("005_structured_template_null", seeds, one_run, OUT / "grid" / args.tag / f"part{args.part}",
              params={"step": "grid", "tag": args.tag, "axes": {k: [str(x) for x in v] for k, v in AXES.items()},
                      "reps": REPS_A, "targets": TARGETS}, progress=progress)


def load_dir(path: Path) -> pd.DataFrame:
    frames = [load_runs(p.parent) for p in sorted(path.glob("**/runs.csv"))]
    frames = [f for f in frames if not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def real_point(tag: str) -> pd.Series:
    path = OUT / "real" / tag / "point.csv"
    if not path.exists():
        ctx = context(tag)
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{"tag": tag, **all_metrics(ctx.real, ctx, 0)}]).to_csv(path, index=False, float_format="%.6g")
    return pd.read_csv(path).iloc[0]


def add_loss(means: pd.DataFrame, target: pd.Series) -> pd.DataFrame:
    out = means.copy()
    out["loss"] = 0.0
    for f, w in TARGETS.items():
        comp = ((out[f] - float(target[f])) / float(target[f])) ** 2
        out[f"comp_{f}"] = comp
        out["loss"] += w * comp
    return out


def config_means(runs: pd.DataFrame) -> pd.DataFrame:
    num = [c for c in runs.columns if c not in ("tag", "seed", "rep", "config_id", *AXES, "suffix_by", "variant")
           and pd.api.types.is_numeric_dtype(runs[c])]
    g = runs.groupby(["tag", "config_id"], sort=False)
    means = g[num].mean()
    means.insert(0, "n_seeds", g.size())
    cfg = g[list(AXES)].first()
    return cfg.join(means).reset_index()


def cmd_rank_grid(args) -> None:
    runs = load_dir(OUT / "grid")
    rows = []
    for tag, sub in runs.groupby("tag", sort=False):
        m = add_loss(config_means(sub), real_point(tag)).sort_values("loss", kind="stable").reset_index(drop=True)
        m.insert(2, "rank", np.arange(1, len(m) + 1))
        rows.append(m)
    grid = pd.concat(rows, ignore_index=True)
    lead = ["tag", "rank", "config_id", *AXES, "n_seeds", "loss", *[f"comp_{f}" for f in TARGETS]]
    grid = grid[lead + [c for c in grid.columns if c not in lead]]
    grid.to_csv(OUT / "grid_summary.csv", index=False, float_format="%.6g")
    for tag, sub in grid.groupby("tag", sort=False):
        print(tag, f"configs {len(sub)}, runs {int(sub['n_seeds'].sum())}")
        print(sub.head(8)[["rank", "config_id", *AXES, "loss"] + list(TARGETS)].to_string(index=False))


def top_configs(tag: str, n: int = TOP_B) -> list[int]:
    grid = pd.read_csv(OUT / "grid_summary.csv")
    return grid[grid["tag"] == tag].sort_values("rank").head(n)["config_id"].astype(int).tolist()


def cmd_phase_b(args) -> None:
    """5 best grid configurations of ``tag`` x 50 seeds; ``--on`` evaluates them on another text."""
    on = args.on or args.tag
    ids = top_configs(args.tag)[: args.top]
    seeds = [cid * 100 + r for cid in ids for r in range(REPS_B)]
    seeds = [s for k, s in enumerate(seeds) if k % args.parts == args.part]
    _warm(on, range(REPS_B))

    def one_run(seed: int) -> dict:
        cid, rep = divmod(seed, 100)
        return {"calibrated_on": args.tag, **template_run(on, CONFIGS[cid], cid, rep, seed)}

    run_seeds("005_structured_template_null", seeds, one_run, OUT / "phaseB" / on / f"part{args.part}",
              params={"step": "phaseB", "tag": on, "calibrated_on": args.tag, "configs": ids, "reps": REPS_B},
              progress=progress)


# --------------------------------------------------------------------------- #
# step 2
# --------------------------------------------------------------------------- #


def best_config(tag: str) -> tuple[int, dict]:
    """Rank 1 of phase B (loss over 50 seeds) for ``tag``."""
    pb = pd.read_csv(OUT / "best_configs.csv")
    r = pb[pb["tag"] == tag].sort_values("rank").iloc[0]
    cfg = {k: r[k] for k in AXES}
    cfg["n_stems"] = int(cfg["n_stems"])
    cfg["zipf_s"] = float(cfg["zipf_s"])
    cfg["topic_concentration"] = float(cfg["topic_concentration"])
    cfg["page_locality"] = float(cfg["page_locality"])
    return int(r["config_id"]), cfg


def cmd_ablation(args) -> None:
    cid0, cfg0 = best_config(args.tag)
    variants = args.variant.split(",") if args.variant else list(VARIANTS)
    _warm(args.tag, range(REPS_B))
    for v in variants:
        cfg = {**cfg0, "suffix_by": None, **VARIANTS[v]}
        vid = VARIANT_ID[v]
        seeds = [vid * 100 + r for r in range(REPS_B)]
        print(f"[{args.tag}] {v}: {cfg}", flush=True)

        def one_run(seed: int, cfg=cfg, v=v) -> dict:
            _, rep = divmod(seed, 100)
            row = template_run(args.tag, cfg, cid0, rep, seed)
            return {"variant": v, **row}

        run_seeds("005_structured_template_null", seeds, one_run, OUT / "ablation" / args.tag / v,
                  params={"step": "ablation", "tag": args.tag, "variant": v, "best_config_id": cid0,
                          "config": config_columns(cfg)}, progress=progress)


# --------------------------------------------------------------------------- #
# step 3
# --------------------------------------------------------------------------- #


def _bpe_decompose(df: pd.DataFrame, merges, preset, prefixes: list[str], suffixes: list[str]) -> pd.DataFrame:
    """Prefix = first BPE unit, suffix = last unit, if they are in the lists and a stem unit remains."""
    pre, suf = set(prefixes), set(suffixes)
    units = apply_bpe(df["token"].tolist(), merges, preset=preset)
    P, S, X = [], [], []
    for u in units:
        p = u[0] if len(u) >= 2 and u[0] in pre else ""
        rest = u[1:] if p else u
        s = rest[-1] if len(rest) >= 2 and rest[-1] in suf else ""
        stem = rest[:-1] if s else rest
        P.append(p)
        S.append("".join(stem))
        X.append(s)
    out = df.copy()
    out["prefix"], out["stem"], out["suffix"] = P, S, X
    out.attrs.update(df.attrs)
    out.attrs["prefixes"], out.attrs["suffixes"] = list(prefixes), list(suffixes)
    return out


class RobustContext:
    """Real text (split policy) of one file and the decomposition of every robustness variant."""

    def __init__(self, key: str):
        self.ctx = context(f"{key}_split")
        self.real = self.ctx.real
        preset = self.real.attrs.get("glyph_preset")
        self.preset = preset
        self.decomp = {}
        for v, opt in ROBUST.items():
            if "bpe" in opt:
                merges = learn_bpe(self.real["token"].tolist(), merges=opt["bpe"], preset=preset)
                units = apply_bpe(self.real["token"].tolist(), merges, preset=preset)
                first = pd.Series([u[0] for u in units if len(u) >= 2]).value_counts()
                last = pd.Series([u[-1] for u in units if len(u) >= 2]).value_counts()
                # as many affixes as the learned lists of 002 (12 prefixes, 14 suffixes), most frequent units
                self.decomp[v] = {"bpe": merges, "prefixes": first.head(12).index.tolist(),
                                  "suffixes": last.head(14).index.tolist(), "order": "prefix_first"}
            else:
                kw = {k: opt[k] for k in ("n_suffixes", "max_len") if k in opt}
                d = decompose_table(self.real, **kw)
                self.decomp[v] = {"prefixes": d.attrs["prefixes"], "suffixes": d.attrs["suffixes"],
                                  "order": opt.get("order", "prefix_first")}

    def decompose(self, df: pd.DataFrame, v: str) -> pd.DataFrame:
        d = self.decomp[v]
        if "bpe" in d:
            return _bpe_decompose(df, d["bpe"], self.preset, d["prefixes"], d["suffixes"])
        return decompose_table(df, prefixes=d["prefixes"], suffixes=d["suffixes"], order=d["order"])


def page_subsample(df: pd.DataFrame, seed: int, frac: float = SUB_FRAC) -> pd.DataFrame:
    """``frac`` of the pages drawn without replacement, text order kept."""
    rng = np.random.default_rng(seed)
    folios = df["folio"].drop_duplicates().to_numpy()
    keep = set(rng.choice(folios, size=int(round(frac * len(folios))), replace=False).tolist())
    out = df[df["folio"].isin(keep)].reset_index(drop=True)
    out.attrs.update(df.attrs)
    return out


def robust_row(rc: RobustContext, source: str, seed: int) -> dict:
    if source == "real":
        base = rc.real
    elif source == "real_sub":
        base = page_subsample(rc.real, seed)
    elif source == "autocopy_536":
        base = autocopy_like(rc.real, allowed_from=rc.ctx.allowed_from, seed=seed, **E2.ACFG["autocopy_536"])
    elif source == "slot_indep":
        base = None
    else:
        raise ValueError(source)
    row = {}
    for v in ROBUST:
        if base is None:  # slot null of this decomposition, re-parsed with the same lists
            df = slot_independent_sample(rc.decompose(rc.real, v), seed=seed)
        else:
            df = base
        m = metrics_decomposed(rc.decompose(df, v), seed)
        for k in ROBUST_METRICS:
            row[f"{v}__{k}"] = m[k]
    return row


def cmd_robust(args) -> None:
    rc = RobustContext(args.file)
    sources = args.source.split(",") if args.source else list(ROBUST_SOURCES)
    for source in sources:
        seeds = list(ROBUST_SOURCES[source])
        if args.seeds:
            keep = set(parse_seeds(args.seeds))
            seeds = [s for s in seeds if s in keep]

        def one_run(seed: int, source=source) -> dict:
            return {"file": args.file, "source": source, **robust_row(rc, source, seed)}

        run_seeds("005_structured_template_null", seeds, one_run, OUT / "robust" / args.file / source,
                  params={"step": "robustness", "file": args.file, "source": source, "subsample_frac": SUB_FRAC,
                          "variants": {v: {k: (x if k != "bpe" else len(x)) for k, x in d.items()}
                                       for v, d in rc.decomp.items()}},
                  progress=progress)


# --------------------------------------------------------------------------- #
# step 4 (real text): page-subsample intervals and page-label permutation baselines
# --------------------------------------------------------------------------- #


def _perm_null(dfd: pd.DataFrame, n_perm: int, seed: int) -> dict:
    """Page-label permutation baseline of every contrast: labels shuffled across the pages involved."""
    rng = np.random.default_rng(seed)
    lh = lang_hand(dfd)
    lang = dfd["language"].fillna("?").astype(str).to_numpy()
    sec = dfd["illustration"].fillna("?").astype(str).to_numpy()
    page = dfd["folio"].astype(str).to_numpy()
    vals = dfd["suffix"].astype(str).to_numpy()
    support = sorted(set(vals.tolist()))
    out: dict[str, list] = {}

    def perm_pair(key, a, b, name, restrict=None):
        m = (key == a) | (key == b)
        if restrict is not None:
            m &= restrict
        pages = pd.unique(page[m])
        lab = pd.Series(key[m]).groupby(page[m]).first().reindex(pages).to_numpy()
        res = []
        for _ in range(n_perm):
            perm = dict(zip(pages, rng.permutation(lab)))
            k2 = np.array([perm.get(p, "") for p in page])
            res.append(_jsd_of(vals, m & (k2 == a), m & (k2 == b), support))
        out[name] = res

    for name, kind, a, b in AB_PAIRS:
        perm_pair(lang if kind == "language" else lh, a, b, f"ab_suffix_{name}")
    for h in SECTION_HANDS:
        m = lh == h
        secs = [s for s, c in pd.Series(sec[m]).value_counts().items() if c >= MIN_SECTION_TOKENS]
        pages = pd.unique(page[m])
        lab = pd.Series(sec[m]).groupby(page[m]).first().reindex(pages).to_numpy()
        res = []
        for _ in range(n_perm):
            perm = dict(zip(pages, rng.permutation(lab)))
            s2 = np.array([perm.get(p, "") for p in page])
            pairs = [_jsd_of(vals, m & (s2 == x), m & (s2 == y), support) for i, x in enumerate(secs) for y in secs[i + 1:]]
            res.append(float(np.mean(pairs)) if pairs else float("nan"))
        out[f"ab_suffix_sections_within_{h}"] = res
    return out


def cmd_ab_real(args) -> None:
    for tag in args.tags.split(","):
        ctx = context(tag)
        d = ctx.real_d

        def one_run(seed: int) -> dict:
            if seed == 0:
                x = d
            else:
                x = decompose_table(page_subsample(ctx.real, seed), prefixes=ctx.prefixes, suffixes=ctx.suffixes)
            return {"tag": tag, "source": "real" if seed == 0 else "real_sub", **ab_contrasts(x)}

        run_seeds("005_structured_template_null", list(range(201)), one_run, OUT / "ab" / tag / "real",
                  params={"step": "ab_real", "tag": tag, "subsample_frac": SUB_FRAC, "note": "seed 0 = full text"},
                  progress=progress)
        perm = _perm_null(d, N_PERM, seed=12345)
        rows = [{"tag": tag, "contrast": k, "perm_mean": float(np.nanmean(v)), "perm_sd": float(np.nanstd(v, ddof=1)),
                 "perm_q975": float(np.nanquantile(v, 0.975)), "n_perm": len(v)} for k, v in perm.items()]
        pd.DataFrame(rows).to_csv(OUT / "ab" / tag / "perm.csv", index=False, float_format="%.6g")
        print(pd.DataFrame(rows).to_string(index=False))


# --------------------------------------------------------------------------- #
# check of the copied metric function against 002
# --------------------------------------------------------------------------- #


def cmd_check_002(args) -> None:
    r2 = pd.read_csv(R002 / "runs.csv")
    res = {}
    for tag in TAGS:
        ctx = context(tag)
        mine = metrics_decomposed(ctx.real_d, 0)
        ref = r2[(r2["combo"] == f"{tag}_learned") & (r2["source"] == "real")].iloc[0]
        res[tag] = max(abs(float(mine[k]) - float(ref[k])) for k in mine if not (pd.isna(mine[k]) and pd.isna(ref[k])))
    print(json.dumps(res, indent=2))
    (OUT / "check_002.json").write_text(json.dumps({"max_abs_diff_vs_002_real": res}, indent=2), encoding="utf-8")


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #

MAIN_METRICS = ["cmi_excess", "mi_prefix_suffix", "mean_suffixes_per_stem", "frac_stems_multi_suffix",
                "frac_multi_suffix_5to9", "H_suffix_given_stem", "stem_dispersion", "jsd_stem_sections",
                "jsd_suffix_sections", "jsd_suffix_AB", "jaccard_stem_sections", "n_stems_k", "productivity_slope"]


def z_table(runs: pd.DataFrame, by: list[str], point: dict[str, pd.Series], metrics: list[str]) -> pd.DataFrame:
    rows = []
    for key, sub in runs.groupby(by, sort=False):
        key = key if isinstance(key, tuple) else (key,)
        tag = dict(zip(by, key))["tag"]
        for m in metrics:
            v = sub[m].astype(float)
            mu, sd = float(v.mean()), float(v.std(ddof=1))
            real = float(point[tag][m])
            rows.append({**dict(zip(by, key)), "metric": m, "real": real, "mean": mu, "sd": sd,
                         "q025": float(v.quantile(0.025)), "q975": float(v.quantile(0.975)),
                         "diff": real - mu, "z": (real - mu) / sd if sd > 0 else float("nan"), "n": len(v)})
    return pd.DataFrame(rows)


def cmd_summarize(args) -> None:
    points = {t: real_point(t) for t in TAGS}
    runs_all = []
    for step in ("grid", "phaseB", "ablation"):
        r = load_dir(OUT / step)
        if not r.empty:
            r.insert(0, "step", step)
            runs_all.append(r)
    runs = pd.concat(runs_all, ignore_index=True)
    lead = ["step", "tag", "variant", "calibrated_on", "config_id", "rep", "seed"]
    lead = [c for c in lead if c in runs.columns]
    runs = runs[lead + [c for c in runs.columns if c not in lead]]
    runs.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")

    # phase B ranking -> best_configs.csv
    pb = runs[runs["step"] == "phaseB"]
    best = []
    for tag, sub in pb.groupby("tag", sort=False):
        m = add_loss(config_means(sub), points[tag]).sort_values("loss", kind="stable").reset_index(drop=True)
        m.insert(2, "rank", np.arange(1, len(m) + 1))
        best.append(m)
    best = pd.concat(best, ignore_index=True)
    best["meets_criterion1"] = best["loss"] < LOSS_LIMIT
    lead = ["tag", "rank", "config_id", *AXES, "n_seeds", "loss", "meets_criterion1", *[f"comp_{f}" for f in TARGETS]]
    best = best[lead + [c for c in best.columns if c not in lead]]
    best.to_csv(OUT / "best_configs.csv", index=False, float_format="%.6g")

    # ablation table (z of the real text vs every variant)
    ab = runs[runs["step"] == "ablation"]
    abl = pd.DataFrame()
    if not ab.empty:
        metrics = list(TARGETS) + [m for m in MAIN_METRICS if m not in TARGETS]
        abl = z_table(ab, ["tag", "variant"], points, metrics)
        abl["reproduced"] = abl["z"].abs() < 3
        abl.to_csv(OUT / "ablation.csv", index=False, float_format="%.6g")
        ingredient_table(abl).to_csv(OUT / "ingredients.csv", index=False)

    robustness_table()
    ab_by_hand_table(runs)
    meta = {"tags": TAGS, "axes": {k: [str(x) for x in v] for k, v in AXES.items()}, "targets": TARGETS,
            "variants": {k: {a: str(b) for a, b in v.items()} for k, v in VARIANTS.items()},
            "robust_variants": ROBUST, "robust_sources": {k: len(v) for k, v in ROBUST_SOURCES.items()},
            "ab_pairs": AB_PAIRS, "section_hands": SECTION_HANDS, "min_section_tokens": MIN_SECTION_TOKENS}
    write_meta(OUT, "005_structured_template_null", params=meta, n_runs=len(runs),
               runs_per_step=runs.groupby("step").size().to_dict())
    info = checks(runs)
    (OUT / "summary_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps(info, indent=2, ensure_ascii=False, default=str))
    if not abl.empty:
        plot(abl)


def ingredient_table(abl: pd.DataFrame) -> pd.DataFrame:
    """Metric x ingredient: reproduced by the best template (ZL and GC) and which ablation breaks it."""
    rows = []
    for m, sub in abl.groupby("metric", sort=False):
        piv = sub.pivot(index="variant", columns="tag", values="z")
        row = {"metric": m}
        for tag in MAIN_TAGS:
            if tag in piv.columns:
                row[f"z_best_{tag}"] = round(float(piv.loc["best", tag]), 2)
        both = lambda v: all(abs(piv.loc[v, t]) < 3 for t in MAIN_TAGS if t in piv.columns)  # noqa: E731
        row["best_reproduces"] = both("best")
        row["broken_by"] = ",".join(v for v in VARIANTS if v.startswith(("no_", "free_")) and both("best") and not both(v))
        row["fixed_by"] = ",".join(v for v in VARIANTS if v.startswith("add_") and not both("best") and both(v))
        for v in VARIANTS:
            row[f"absz_max_{v}"] = round(float(piv.loc[v].abs().max()), 2) if v in piv.index else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def robustness_table() -> None:
    """Per file x variant: real point, subsample SD (scaled to full-text SD), nulls, verdicts (i)-(iii)."""
    rows = []
    for key in ("ZL", "GC"):
        parts = {s: load_runs(OUT / "robust" / key / s) for s in ROBUST_SOURCES}
        if any(p.empty for p in parts.values()):
            continue
        scale = math.sqrt(SUB_FRAC / (1 - SUB_FRAC))  # delete-d jackknife: SD(full) ~ SD(80 % subsample) * 2
        for v in ROBUST:
            for m in ROBUST_METRICS:
                c = f"{v}__{m}"
                real = float(parts["real"][c].iloc[0])
                sd_r = float(parts["real_sub"][c].std(ddof=1)) * scale
                row = {"file": key, "variant": v, "metric": m, "real": real, "real_sd": sd_r}
                for null in ("slot_indep", "autocopy_536"):
                    x = parts[null][c].astype(float)
                    sd = math.sqrt(sd_r ** 2 + float(x.std(ddof=1)) ** 2)
                    d = real - float(x.mean())
                    row.update({f"{null}_mean": float(x.mean()), f"{null}_sd": float(x.std(ddof=1)),
                                f"diff_{null}": d, f"lo_{null}": d - 1.96 * sd, f"hi_{null}": d + 1.96 * sd})
                rows.append(row)
    if not rows:
        return
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "robustness.csv", index=False, float_format="%.6g")
    verdicts = []
    for (key, v), sub in t.groupby(["file", "variant"], sort=False):
        s = sub.set_index("metric")
        # (i) cmi_excess above slot (lower bound > 0) and below both autocopy (upper bound < 0);
        #     paradigm richer than autocopy, narrower than slot
        i_cmi = s.loc["cmi_excess", "lo_slot_indep"] > 0 and s.loc["cmi_excess", "hi_autocopy_536"] < 0
        i_par = (s.loc["mean_suffixes_per_stem", "lo_autocopy_536"] > 0 and s.loc["mean_suffixes_per_stem", "hi_slot_indep"] < 0
                 and s.loc["frac_multi_suffix_5to9", "lo_autocopy_536"] > 0)
        i_mi = s.loc["mi_prefix_suffix", "lo_slot_indep"] > 0 and s.loc["mi_prefix_suffix", "lo_autocopy_536"] > 0
        # (ii) affix inventory more stable than stem inventory, and stem JSD between slot and autocopy
        ii_gap = s.loc["jaccard_gap_sections", "real"] > 0
        ii_jsd = s.loc["jsd_stem_sections", "lo_slot_indep"] > 0 and s.loc["jsd_stem_sections", "hi_autocopy_536"] < 0
        ii_suf = s.loc["jsd_suffix_sections", "lo_slot_indep"] > 0
        # (iii) dispersion above autocopy, below slot
        iii = s.loc["stem_dispersion", "lo_autocopy_536"] > 0 and s.loc["stem_dispersion", "hi_slot_indep"] < 0
        verdicts.append({"file": key, "variant": v, "i_cmi_between_slot_and_autocopy": bool(i_cmi),
                         "i_paradigm_between": bool(i_par), "i_mi_above_both": bool(i_mi),
                         "ii_jaccard_gap_positive": bool(ii_gap), "ii_jsd_stem_between": bool(ii_jsd),
                         "ii_jsd_suffix_above_slot": bool(ii_suf), "iii_dispersion_between": bool(iii)})
    pd.DataFrame(verdicts).to_csv(OUT / "robustness_verdicts.csv", index=False)


def ab_by_hand_table(runs: pd.DataFrame) -> None:
    rows = []
    cols = [f"ab_suffix_{n}" for n, *_ in AB_PAIRS] + [f"ab_suffix_sections_within_{h}" for h in SECTION_HANDS]
    kinds = {f"ab_suffix_{n}": k for n, k, *_ in AB_PAIRS}
    kinds.update({f"ab_suffix_sections_within_{h}": f"between_sections_within_{h}" for h in SECTION_HANDS})
    scale = math.sqrt(SUB_FRAC / (1 - SUB_FRAC))
    for tag in TAGS:
        real = load_runs(OUT / "ab" / tag / "real")
        perm_path = OUT / "ab" / tag / "perm.csv"
        if real.empty or not perm_path.exists():
            continue
        perm = pd.read_csv(perm_path).set_index("contrast")
        point = real[real["seed"] == 0].iloc[0]
        sub = real[real["seed"] > 0]
        n_tok = token_counts(tag)
        for c in cols:
            sd = float(sub[c].std(ddof=1)) * scale
            x = float(point[c])
            pm = float(perm.loc[c, "perm_mean"])
            row = {"tag": tag, "contrast": c.replace("ab_suffix_", ""), "kind": kinds[c], "n_tokens": n_tok.get(c, ""),
                   "source": "real", "jsd": x, "ci_lo": x - 1.96 * sd, "ci_hi": x + 1.96 * sd,
                   "perm_mean": pm, "perm_q975": float(perm.loc[c, "perm_q975"]), "excess": x - pm,
                   "excess_lo": x - pm - 1.96 * sd, "excess_hi": x - pm + 1.96 * sd}
            rows.append(row)
        ab = runs[(runs["step"] == "ablation") & (runs["tag"] == tag)]
        for v in ("best", "no_topic", "add_suffix_by_language"):
            t = ab[ab["variant"] == v]
            if t.empty:
                continue
            for c in cols:
                x = t[c].astype(float)
                rows.append({"tag": tag, "contrast": c.replace("ab_suffix_", ""), "kind": kinds[c],
                             "n_tokens": n_tok.get(c, ""), "source": f"template_{v}", "jsd": float(x.mean()),
                             "ci_lo": float(x.quantile(0.025)), "ci_hi": float(x.quantile(0.975))})
    if rows:
        pd.DataFrame(rows).to_csv(OUT / "ab_by_hand.csv", index=False, float_format="%.6g")


def token_counts(tag: str) -> dict:
    d = context(tag).real_d
    lh = lang_hand(d)
    lang = d["language"].fillna("?").astype(str).to_numpy()
    out = {}
    for name, kind, a, b in AB_PAIRS:
        key = lang if kind == "language" else lh
        out[f"ab_suffix_{name}"] = f"{int((key == a).sum())}/{int((key == b).sum())}"
    sec = d["illustration"].fillna("?").astype(str).to_numpy()
    for h in SECTION_HANDS:
        vc = pd.Series(sec[lh == h]).value_counts()
        out[f"ab_suffix_sections_within_{h}"] = "/".join(f"{s}:{c}" for s, c in vc.items() if c >= MIN_SECTION_TOKENS)
    return out


def checks(runs: pd.DataFrame) -> dict:
    flags = {"cmi_excess > 1": runs["cmi_excess"] > 1, "stem_dispersion > 10": runs["stem_dispersion"] > 10}
    for m in ("jaccard_suffix_sections", "jaccard_stem_sections", "jaccard_suffix_AB", "jaccard_stem_AB"):
        flags[f"{m} == 1"] = runs[m] == 1.0
        flags[f"{m} == 0"] = runs[m] == 0.0
    sus = {k: runs[v].groupby(runs["step"]).size().to_dict() for k, v in flags.items() if v.any()}
    counts = runs.groupby(["step", "tag"]).size().to_dict()
    return {"runs_per_step_tag": {f"{a}/{b}": int(c) for (a, b), c in counts.items()},
            "suspicious_runs_by_step": sus,
            "share_stems_real_mean": float(runs["share_stems_real"].mean()),
            "check_002": json.loads((OUT / "check_002.json").read_text()) if (OUT / "check_002.json").exists() else None}


def plot(abl: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; template_vs_real.png skipped")
        return
    metrics = [m for m in MAIN_METRICS if m in set(abl["metric"])] + ["n_types", "hapax_share", "mean_len"]
    variants = list(VARIANTS)
    tags = [t for t in MAIN_TAGS if t in set(abl["tag"])]
    fig, axes = plt.subplots(1, len(tags), figsize=(7 * len(tags), 7), sharey=True)
    axes = np.atleast_1d(axes)
    cmap = plt.get_cmap("tab10")
    for ax, tag in zip(axes, tags):
        sub = abl[abl["tag"] == tag]
        y = np.arange(len(metrics))
        for j, v in enumerate(variants):
            z = sub[sub["variant"] == v].set_index("metric")["z"].reindex(metrics).clip(-30, 30)
            ax.scatter(z, y + (j - len(variants) / 2) * 0.08, s=14 if v != "best" else 40, color=cmap(j % 10),
                       marker="o" if v != "best" else "D", label=v, zorder=3 if v == "best" else 2)
        ax.axvspan(-3, 3, color="0.9", zorder=0)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(metrics, fontsize=8)
        ax.set_xlabel("z = (real - template mean) / template SD (clipped at +-30)")
        ax.set_title(tag)
        ax.grid(axis="x", alpha=0.3)
    axes[0].legend(fontsize=7, loc="lower left")
    fig.suptitle("005: real text vs structured template (best config and ablations/additions, 50 seeds each)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "template_vs_real.png", dpi=110)
    print(f"wrote {OUT / 'template_vs_real.png'}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("grid")
    p.add_argument("--tag", required=True, choices=TAGS)
    p.add_argument("--part", type=int, default=0)
    p.add_argument("--parts", type=int, default=1)
    p.set_defaults(fn=cmd_grid)
    sub.add_parser("rank-grid").set_defaults(fn=cmd_rank_grid)
    p = sub.add_parser("phase-b")
    p.add_argument("--tag", required=True, choices=TAGS, help="text whose grid ranking is used")
    p.add_argument("--on", default="", help="evaluate on this text instead (e.g. the join policy)")
    p.add_argument("--top", type=int, default=TOP_B)
    p.add_argument("--part", type=int, default=0)
    p.add_argument("--parts", type=int, default=1)
    p.set_defaults(fn=cmd_phase_b)
    p = sub.add_parser("ablation")
    p.add_argument("--tag", required=True, choices=TAGS)
    p.add_argument("--variant", default="", help="comma-separated subset of " + ",".join(VARIANTS))
    p.set_defaults(fn=cmd_ablation)
    p = sub.add_parser("robust")
    p.add_argument("--file", required=True, choices=("ZL", "GC"))
    p.add_argument("--source", default="", help="comma-separated subset of " + ",".join(ROBUST_SOURCES))
    p.add_argument("--seeds", default="")
    p.set_defaults(fn=cmd_robust)
    p = sub.add_parser("ab-real")
    p.add_argument("--tags", default=",".join(TAGS))
    p.set_defaults(fn=cmd_ab_real)
    sub.add_parser("check-002").set_defaults(fn=cmd_check_002)
    sub.add_parser("summarize").set_defaults(fn=cmd_summarize)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
