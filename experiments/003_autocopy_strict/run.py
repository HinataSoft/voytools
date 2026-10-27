"""Experiment 003: auto-copy with edits bound by a word grammar (trigram / slot).

Usage (from the repository root; see ``spec.md``)::

    # real text: point estimates (seed 0), grammaticality of the real text
    python experiments/003_autocopy_strict/run.py real --file ZL3b-n.txt --space-policy split
    # phase A: 432 configs x 3 seeds on ZL split, target metrics + edge/token MI; 4 parallel chunks
    python experiments/003_autocopy_strict/run.py grid --seeds 0:43200/4:0 --out results/003_autocopy_strict/phaseA/part0
    python experiments/003_autocopy_strict/run.py rankA          # -> grid_summary.csv, phaseA_top10.csv
    # phase B: 10 configs (5 best per edit_model) x 50 seeds, all metrics, per data set
    python experiments/003_autocopy_strict/run.py phaseB --file ZL3b-n.txt --space-policy split
    # grammaticality of the 001 bigram generator (same seeds as 001, so the same samples)
    python experiments/003_autocopy_strict/run.py bigramgram --file ZL3b-n.txt --space-policy split
    # outside the spec: decomposition of the H2 excess (ZL split)
    python experiments/003_autocopy_strict/run.py diag
    # tables, profile, plot
    python experiments/003_autocopy_strict/run.py summarize

Run identification: ``seed = config_id * 100 + rep`` with ``config_id`` 0..431 in the order of
``itertools.product`` over :data:`AXES`.  A ``--seeds`` range is a range of encoded seeds; only
seeds with ``rep < --reps`` are run.  Every command is resumable (``run_seeds``).

The bigram reference (001) is not recomputed: its phase-B runs, the p_new = 0.5 diagnostic point,
the within-line shuffles and the page bootstrap of the real text are read from
``results/001_autocopy_calibration/``.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_runs, load_tables, merge_runs, run_seeds, seed_range  # noqa: E402
from voytools.morphology import bpe_dependence_curve, decompose_table, paradigm_metrics  # noqa: E402
from voytools.nulls import autocopy_like  # noqa: E402
from voytools.stats import (  # noqa: E402
    adjacent_token_mi,
    edge_glyph_mi,
    glyph_entropy_report,
    neighbor_similarity,
    positional_profile,
    vocabulary_summary,
)
from voytools.wordgrammar import grammaticality, learn_word_grammar  # noqa: E402

DATA = ROOT / "data"
OUT = ROOT / "results" / "003_autocopy_strict"
REF001 = ROOT / "results" / "001_autocopy_calibration"

AXES = {
    "edit_model": ("trigram", "slot"),
    "new_token_source": ("seed", "grammar"),
    "p_new": (0.05, 0.15, 0.3),
    "window": (16, 32, 64),
    "p_copy_exact": (0.35, 0.5, 0.65),
    "p_two_edits": (0.25, 0.5),
    "recency": (0.0, 0.3),
}
FIXED = {"edit_type_probs": (0.6, 0.2, 0.2)}
CONFIGS = list(itertools.product(*AXES.values()))  # 432
TARGETS = {"share_near": 2.0, "hapax_share": 2.0, "H2": 2.0, "n_types": 1.0, "mean_len": 1.0}
COMPONENT_LIMIT = 0.05
LOSS_LIMIT = 0.15
P_NEW_MAX_SEED = 0.3
NONTARGET = ["edge_mi_excess", "token_mi_excess_share", "bpe_gap64", "pos_jsd_first_mid",
             "cmi_excess", "frac_stems_multi_suffix", "mean_suffixes_per_stem"]
EXTRA = ["near_ratio"]
CHECKS = ["gram_slot", "gram_trigram", "gram_bigram", "p_unchanged", "p_edit_failed"]
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
TOP_PER_MODEL = 5

# 001 reference: best phase-B config per data set and the p_new = 0.5 diagnostic point (upper bound)
DIAG001 = {"base_config": 573, "p_new": 0.5, "p_copy_exact": 0.5}

# docs/reference_values.md, Rozanova & Temerev (2026) tab. 7 (ZL3b, paragraph text, split-like)
REFERENCE = {
    "ref_RT2026_voynich": {"H2": (2.69, 2.69), "bpe_gap64": (1.05, 1.05), "token_mi_excess_share": (0.0079, 0.0079),
                           "hapax_share": (0.70, 0.70), "edge_mi_excess": (0.197, 0.197), "near_ratio": (1.12, 1.12)},
    "ref_RT2026_selfcitation": {"H2": (2.6, 2.7), "bpe_gap64": (0.65, 0.86), "token_mi_excess_share": (0.005, 0.006),
                                "hapax_share": (0.59, 0.60), "edge_mi_excess": (0.006, 0.009),
                                "near_ratio": (1.06, 1.23)},
}


def parse_seeds(spec: str) -> list[int]:
    """``'0:20'`` -> seeds 0..19; ``'0:100/4:1'`` -> chunk 1 of 4 of 0..99; ``'3,5,8'`` -> those seeds."""
    if "/" in spec:
        rng, part = spec.split("/")
        chunks, k = part.split(":")
        a, b = rng.split(":")
        return seed_range(int(a), int(b), int(chunks))[int(k)]
    if ":" in spec:
        a, b = spec.split(":")
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(",") if s.strip()]


def config_params(cid: int) -> dict:
    """Generator keyword arguments of grid configuration ``cid``."""
    return {**dict(zip(AXES, CONFIGS[cid])), **FIXED}


def config_columns(cid: int) -> dict:
    return {"config_id": cid, **dict(zip(AXES, CONFIGS[cid]))}


def tag_of(file: str, policy: str) -> str:
    return f"{'GC' if file.startswith('GC') else 'ZL'}_{policy}"


# --------------------------------------------------------------------------- #
# data and metrics
# --------------------------------------------------------------------------- #


class Context:
    """Real text of one file/space policy, its word grammar and affix lists (learned once)."""

    def __init__(self, file: str, policy: str):
        self.file, self.policy = file, policy
        self.real = filter_tokens(load_tables(DATA / file, space_policy=policy).tokens, generic_type="P", clean=True)
        self.eva = self.real.attrs.get("glyph_preset") is not None
        self.h2_preset = "extended" if self.eva else None  # as in 001
        self.allowed_from = self.real["token"].tolist()
        # grammar and trigrams learned from the whole real text (generator glyph preset)
        self.grammar = learn_word_grammar(self.allowed_from, self.real.attrs.get("glyph_preset"))
        d = decompose_table(self.real)
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]

    def generate(self, params: dict, seed: int) -> pd.DataFrame:
        return autocopy_like(self.real, allowed_from=self.allowed_from, word_grammar=self.grammar, seed=seed, **params)


def lines_of(df: pd.DataFrame) -> list[list[str]]:
    return df.groupby("locus_id", sort=False)["token"].apply(list).tolist()


def stats_for(df: pd.DataFrame, ctx: Context, seed: int, full: bool) -> dict:
    """Metrics of the spec (same code as 001) plus grammaticality and realised copy counts."""
    near = neighbor_similarity(df, window=5, shuffles=20, seed=seed).iloc[0]
    vs = vocabulary_summary(df).iloc[0]
    ge = glyph_entropy_report(df, preset=ctx.h2_preset, scope="all").iloc[0]
    out = {
        "share_near": float(near["share_near"]),
        "near_ratio": float(near["ratio"]),
        "n_types": float(vs["n_types"]),
        "hapax_share": float(vs["hapax_share"]),
        "mean_len": float(vs["mean_len"]),
        "H2": float(ge["H2"]),
    }
    for m in ("slot", "trigram", "bigram"):
        out[f"gram_{m}"] = float(grammaticality(df, m, grammar=ctx.grammar))
    out["gram_slot_types"] = float(grammaticality(df, "slot", grammar=ctx.grammar, weight="type"))
    info = df.attrs.get("autocopy_info")
    if info:
        out["p_unchanged"] = float(info["p_unchanged"])
        out["p_edit_failed"] = info["n_edit_failed"] / max(info["n_edit_attempted"], 1)
        out["p_new_realised"] = float(info["p_new_realised"])
        out["n_grammar_fallback"] = int(info["n_grammar_fallback"])
    out["edge_mi_excess"] = float(edge_glyph_mi(df, shuffles=20, seed=seed).iloc[0]["excess_bits"])
    out["token_mi_excess_share"] = float(adjacent_token_mi(df, cap=2000, shuffles=20, seed=seed).iloc[0]["excess_share"])
    if full:
        curve = bpe_dependence_curve(lines_of(df), merge_checkpoints=(0, 32, 64, 128)).set_index("merges")["gap"]
        out["bpe_gap64"] = float(curve.get(64, np.nan))
        out["pos_jsd_first_mid"] = float(positional_profile(df, which="first").jsd.loc["line_first", "mid"])
        pm = paradigm_metrics(decompose_table(df, prefixes=ctx.prefixes, suffixes=ctx.suffixes), seed=seed).summary
        out["cmi_excess"] = float(pm["cmi_excess"])
        out["frac_stems_multi_suffix"] = float(pm["frac_stems_multi_suffix"])
        out["mean_suffixes_per_stem"] = float(pm["mean_suffixes_per_stem"])
    return out


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_real(args) -> None:
    ctx = Context(args.file, args.space_policy)
    tag = tag_of(args.file, args.space_policy)
    out = OUT / "real" / tag
    out.mkdir(parents=True, exist_ok=True)
    point = {"file": args.file, "space_policy": args.space_policy, "source": "real",
             **stats_for(ctx.real, ctx, 0, True),
             "gram_slot_types": float(grammaticality(ctx.real, "slot", grammar=ctx.grammar, weight="type")),
             "slot_coverage_fit": ctx.grammar.coverage, "n_slots": len(ctx.grammar.slots),
             "n_trigrams": len(ctx.grammar.trigrams), "repeatable": " ".join(sorted(ctx.grammar.repeatable)),
             "slots": " ".join(ctx.grammar.slots)}
    pd.DataFrame([point]).to_csv(out / "point.csv", index=False)
    # consistency with 001 (same text, same metric code)
    p001 = pd.read_csv(REF001 / "real" / tag / "point.csv").iloc[0]
    diffs = {m: float(point[m]) - float(p001[m]) for m in list(TARGETS) + EXTRA + NONTARGET}
    (out / "check_vs_001.json").write_text(json.dumps(diffs, indent=2))
    print(pd.Series(point).to_string())
    print("max |diff| vs 001:", max(abs(v) for v in diffs.values()))


def _gen_runner(ctx: Context, full: bool):
    def one_run(seed: int) -> dict:
        cid = seed // 100
        res = stats_for(ctx.generate(config_params(cid), seed), ctx, seed, full)
        return {"file": ctx.file, "space_policy": ctx.policy, "source": "autocopy", **config_columns(cid), **res}

    return one_run


def cmd_grid(args) -> None:
    ctx = Context(args.file, args.space_policy)
    seeds = [s for s in parse_seeds(args.seeds) if s % 100 < args.reps and s // 100 < len(CONFIGS)]
    run_seeds("003_phaseA", seeds, _gen_runner(ctx, False), args.out or OUT / "phaseA" / "part0",
              params={"file": args.file, "space_policy": args.space_policy, "seeds": args.seeds, "reps": args.reps,
                      "metrics": "target + edge/token MI + grammaticality",
                      "axes": {k: list(v) for k, v in AXES.items()}, "fixed": {k: list(v) for k, v in FIXED.items()}},
              progress=progress)


def cmd_phase_b(args) -> None:
    ctx = Context(args.file, args.space_policy)
    top = pd.read_csv(OUT / "phaseA_top10.csv")["config_id"].astype(int).tolist()
    seeds = [cid * 100 + r for cid in top for r in range(args.reps)]
    if args.seeds:
        keep = set(parse_seeds(args.seeds))
        seeds = [s for s in seeds if s in keep]
    out = args.out or OUT / "phaseB" / tag_of(args.file, args.space_policy)
    run_seeds("003_phaseB", seeds, _gen_runner(ctx, True), out,
              params={"file": args.file, "space_policy": args.space_policy, "configs": top, "reps": args.reps},
              progress=progress)


def cmd_bigramgram(args) -> None:
    """Grammaticality of the 001 bigram generator (best 001 config of this data set, same seeds as 001)."""
    ctx = Context(args.file, args.space_policy)
    tag = tag_of(args.file, args.space_policy)
    best = int(pd.read_csv(REF001 / "phaseB_configs.csv").query("dataset == @tag and rank == 1")["config_id"].iloc[0])
    axes001 = {"window": (4, 8, 16, 32), "p_copy_exact": (0.1, 0.2, 0.35, 0.5), "p_new": (0.02, 0.05, 0.15),
               "p_two_edits": (0.1, 0.25, 0.5), "recency": (0.0, 0.3), "edit_type_probs": ((0.6, 0.2, 0.2), (0.4, 0.3, 0.3))}
    params = dict(zip(axes001, list(itertools.product(*axes001.values()))[best]))

    def one_run(seed: int) -> dict:
        df = autocopy_like(ctx.real, allowed_from=ctx.allowed_from, seed=seed, **params)  # 001 call, unchanged
        res = {f"gram_{m}": float(grammaticality(df, m, grammar=ctx.grammar)) for m in ("slot", "trigram", "bigram")}
        res["gram_slot_types"] = float(grammaticality(df, "slot", grammar=ctx.grammar, weight="type"))
        res["p_unchanged"] = float(df.attrs["autocopy_info"]["p_unchanged"])
        info = df.attrs["autocopy_info"]
        res["p_edit_failed"] = info["n_edit_failed"] / max(info["n_edit_attempted"], 1)
        res["H2_check"] = float(glyph_entropy_report(df, preset=ctx.h2_preset, scope="all").iloc[0]["H2"])
        return {"file": ctx.file, "space_policy": ctx.policy, "source": "bigram_001", "config_id_001": best, **res}

    seeds = [best * 100 + r for r in range(args.reps)]
    run_seeds("003_bigram_grammaticality", seeds, one_run, OUT / "bigram001" / tag,
              params={"file": args.file, "space_policy": args.space_policy, "config_id_001": best}, progress=progress)


DIAG = [  # outside the grid: where does the H2 excess come from? (base = best trigram/slot config of phase B)
    ("no_edits", {"p_copy_exact": 1.0}),
    ("grammar_words_only", {"p_new": 1.0, "p_copy_exact": 1.0, "new_token_source": "grammar"}),
    ("seed_vocab_only", {"p_new": 1.0, "p_copy_exact": 1.0, "new_token_source": "seed"}),
    ("base", {}),
    ("edits_only_no_copy", {"window": 1, "p_copy_exact": 0.0, "p_new": 1.0, "new_token_source": "seed"}),
]


def cmd_diag(args) -> None:
    """Out of spec: H2 / types of decomposed generator variants around the phase-B choice (ZL split)."""
    ctx = Context(args.file, args.space_policy)
    pb = pd.read_csv(OUT / "phaseB_configs.csv").query("dataset == 'ZL_split'")
    bases = [int(pb[pb["edit_model"] == em].iloc[0]["config_id"]) for em in AXES["edit_model"]]

    def one_run(seed: int) -> dict:
        bi, di = divmod(seed // 100, len(DIAG))
        name, over = DIAG[di]
        params = {**config_params(bases[bi]), **over}
        res = stats_for(ctx.generate(params, seed), ctx, seed, False)
        return {"file": ctx.file, "space_policy": ctx.policy, "base_config": bases[bi], "variant": name,
                **{k: str(v) for k, v in params.items()}, **res}

    seeds = [i * 100 + r for i in range(len(bases) * len(DIAG)) for r in range(args.reps)]
    run_seeds("003_diag", seeds, one_run, OUT / "diag" / tag_of(args.file, args.space_policy),
              params={"bases": bases, "variants": [d[0] for d in DIAG], "reps": args.reps}, progress=progress)
    runs = load_runs(OUT / "diag" / tag_of(args.file, args.space_policy))
    keys = ["base_config", "edit_model", "variant"]
    means = runs.groupby(keys, sort=False)[list(TARGETS) + ["gram_slot", "gram_trigram", "p_unchanged"]].mean()
    means = add_loss(means, real_point(tag_of(args.file, args.space_policy)))
    means.to_csv(OUT / "diag_summary.csv", float_format="%.6g")
    print(means.to_string(float_format=lambda x: f"{x:.4g}"))


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #


def add_loss(means: pd.DataFrame, target: pd.Series) -> pd.DataFrame:
    out = means.copy()
    out["loss"] = 0.0
    out["n_missed"] = 0
    for f, w in TARGETS.items():
        comp = ((out[f] - target[f]) / target[f]) ** 2
        out[f"comp_{f}"] = comp
        out["loss"] += w * comp
        out["n_missed"] += (comp > COMPONENT_LIMIT).astype(int)
    return out


def config_means(runs: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    keys = ["config_id"] + list(AXES)
    metrics = [m for m in metrics if m in runs.columns]
    g = runs.groupby(keys, sort=True)
    means = g[metrics].mean()
    means.insert(0, "n_seeds", g.size())
    return means.reset_index()


def suspicious(runs: pd.DataFrame) -> pd.DataFrame:
    """Spec criterion 5 for samples: share_near > 0.9 or H2 < 2.0 (plus impossible values)."""
    bad = (runs["share_near"] > 0.9) | (runs["H2"] < 2.0) | (runs["H2"] > 4) | (runs["hapax_share"] < 0)
    if "gram_trigram" in runs and "edit_model" in runs:
        bad |= (runs["edit_model"] == "trigram") & (runs["gram_trigram"] < 0.999)
    return runs[bad]


def real_point(tag: str) -> pd.Series:
    return pd.read_csv(OUT / "real" / tag / "point.csv").iloc[0]


def cmd_rank_a(args) -> None:
    parts = sorted((OUT / "phaseA").glob("part*"))
    runs_a = merge_runs(parts, OUT)
    metrics_a = list(TARGETS) + EXTRA + ["edge_mi_excess", "token_mi_excess_share"] + CHECKS + ["gram_slot_types"]
    grid = add_loss(config_means(runs_a, metrics_a), real_point("ZL_split"))
    grid = grid.sort_values("loss", kind="stable").reset_index(drop=True)
    grid.insert(0, "rank", np.arange(1, len(grid) + 1))
    grid.to_csv(OUT / "grid_summary.csv", index=False, float_format="%.6g")
    top = grid.groupby("edit_model", sort=False).head(TOP_PER_MODEL).sort_values("loss", kind="stable")
    top.to_csv(OUT / "phaseA_top10.csv", index=False, float_format="%.6g")
    print(grid.head(20).to_string(float_format=lambda x: f"{x:.4g}"))
    print(top[["rank", "config_id"] + list(AXES) + ["loss"]].to_string())


def _ref001_sources(tag: str) -> list[tuple[str, pd.DataFrame, str]]:
    """001 bigram generator: best config of this data set and the p_new = 0.5 point (upper bound)."""
    out = []
    pb = pd.read_csv(REF001 / "phaseB_configs.csv")
    best = int(pb[(pb["dataset"] == tag) & (pb["rank"] == 1)]["config_id"].iloc[0])
    runs = load_runs(REF001 / "phaseB" / tag)
    r = runs[runs["config_id"] == best].copy()
    gp = OUT / "bigram001" / tag / "runs.csv"
    if gp.exists():  # grammaticality of the same samples (recomputed here, metrics unchanged)
        g = pd.read_csv(gp).set_index("seed")
        for c in ("gram_slot", "gram_trigram", "gram_bigram", "p_unchanged", "p_edit_failed"):
            r[c] = r["seed"].map(g[c])
    out.append(("bigram_001", r, f"001:{best}"))
    d = REF001 / "diag_full" / tag
    if (d / "runs.csv").exists():
        out.append(("bigram_001_pnew05_upper", load_runs(d), f"001:{DIAG001['base_config']}/p_new=0.5"))
    return out


def _source_row(row: dict, name: str, r: pd.DataFrame, cid, m: str, x: float) -> dict:
    if m not in r.columns or r[m].isna().all():
        return {**row, "source": name, "config_id": cid, "n": 0, "mean": np.nan, "lo": np.nan, "hi": np.nan,
                "sd": np.nan, "interval": "not computed", "diff_real": np.nan, "z_real": np.nan}
    v = r[m].dropna().astype(float)
    sd = v.std(ddof=1)
    return {**row, "source": name, "config_id": cid, "n": len(v), "mean": v.mean(),
            "lo": v.quantile(0.025), "hi": v.quantile(0.975), "sd": sd, "interval": "2.5/97.5 % over seeds",
            "diff_real": x - v.mean(), "z_real": (x - v.mean()) / sd if sd > 1e-9 * max(abs(x), 1) else np.nan}


def cmd_summarize(args) -> None:
    grid = pd.read_csv(OUT / "grid_summary.csv")
    runs_a = load_runs(OUT)
    info = {"phaseA_runs": len(runs_a), "phaseA_configs": int(runs_a["config_id"].nunique()),
            "phaseA_suspicious_runs": len(suspicious(runs_a)),
            "phaseA_min_gram_trigram_trigram_model": float(runs_a.query("edit_model == 'trigram'")["gram_trigram"].min()),
            "phaseA_min_gram_slot_slot_model": float(runs_a.query("edit_model == 'slot'")["gram_slot"].min()),
            "phaseA_grammar_fallbacks": int(runs_a["n_grammar_fallback"].sum())}
    tags = sorted([p.name for p in (OUT / "phaseB").glob("*_*") if (p / "runs.csv").exists()],
                  key=lambda t: (t != "ZL_split", t))
    metrics_b = list(TARGETS) + EXTRA + NONTARGET + CHECKS
    best_rows, profile_rows, summary_rows = [], [], []
    choice: dict[str, int] = {}
    for tag in tags:
        target = real_point(tag)
        runs_b = load_runs(OUT / "phaseB" / tag)
        info[f"phaseB_{tag}_runs"] = len(runs_b)
        info[f"phaseB_{tag}_suspicious_runs"] = len(suspicious(runs_b))
        info[f"phaseB_{tag}_grammar_fallbacks"] = int(runs_b["n_grammar_fallback"].sum())
        mb = add_loss(config_means(runs_b, metrics_b + ["gram_slot_types", "p_new_realised"]), target)
        mb = mb.sort_values("loss", kind="stable").reset_index(drop=True)
        mb.insert(0, "dataset", tag)
        mb.insert(1, "rank", np.arange(1, len(mb) + 1))
        best_rows.append(mb)
        if tag == "ZL_split":  # the null model is chosen on the primary data set and kept fixed
            for em in AXES["edit_model"]:
                choice[em] = int(mb[mb["edit_model"] == em].iloc[0]["config_id"])
            choice["recommended"] = int(mb.iloc[0]["config_id"])
        boot = load_runs(REF001 / "real" / tag / "bootstrap")
        shuf = load_runs(REF001 / "shuffle" / tag)
        sources = [(f"{em}_best", runs_b[runs_b["config_id"] == choice[em]], choice[em]) for em in AXES["edit_model"]]
        sources += _ref001_sources(tag) + [("shuffle_lines_001", shuf, "")]
        for m in metrics_b:
            x = float(target[m]) if m in target and pd.notna(target[m]) else np.nan
            kind = "target" if m in TARGETS else "extra" if m in EXTRA else "nontarget" if m in NONTARGET else "check"
            row = {"dataset": tag, "metric": m, "kind": kind}
            has_boot = m in boot.columns
            sd_b = boot[m].std(ddof=1) if has_boot else np.nan
            profile_rows.append({**row, "source": "real", "config_id": "", "n": 1, "mean": x,
                                 "lo": x - 1.96 * sd_b if has_boot else np.nan, "hi": x + 1.96 * sd_b if has_boot else np.nan,
                                 "sd": sd_b, "interval": "001 page bootstrap, normal (x +- 1.96 SD)" if has_boot else "",
                                 "z_real": np.nan, "diff_real": np.nan})
            for name, r, cid in sources:
                if name.startswith("shuffle") and m in CHECKS and m.startswith("p_"):
                    continue
                profile_rows.append(_source_row(row, name, r, cid, m, x))
            if tag == "ZL_split":
                for ref, vals in REFERENCE.items():
                    if m in vals:
                        lo, hi = vals[m]
                        profile_rows.append({**row, "source": ref, "config_id": "", "n": np.nan, "mean": (lo + hi) / 2,
                                             "lo": lo, "hi": hi, "sd": np.nan,
                                             "interval": "docs/reference_values.md (range if given)",
                                             "z_real": np.nan, "diff_real": x - (lo + hi) / 2})
        best = choice["recommended"]
        summ = runs_b[runs_b["config_id"] == best][metrics_b].describe(percentiles=[0.025, 0.5, 0.975]).T
        summ.insert(0, "dataset", tag)
        summ.insert(1, "config_id", best)
        summary_rows.append(summ.reset_index().rename(columns={"index": "metric"}))
    best_all = pd.concat(best_rows, ignore_index=True)
    best_all.to_csv(OUT / "phaseB_configs.csv", index=False, float_format="%.6g")
    profile = pd.DataFrame(profile_rows)
    profile.to_csv(OUT / "profile.csv", index=False, float_format="%.6g")
    pd.concat(summary_rows, ignore_index=True).to_csv(OUT / "summary.csv", index=False, float_format="%.6g")

    # best_configs.csv: phase-B ranking on ZL split; criterion 1 decides `provisional`
    zl = best_all[best_all["dataset"] == "ZL_split"].copy()
    pa = grid.set_index("config_id")[["rank", "loss"]].rename(columns={"rank": "phaseA_rank", "loss": "phaseA_loss"})
    zl = zl.join(pa, on="config_id")
    comps = [f"comp_{f}" for f in TARGETS]
    ok = (zl["loss"] < LOSS_LIMIT) & (zl[comps] < COMPONENT_LIMIT).all(axis=1) & (
        (zl["new_token_source"] != "seed") | (zl["p_new"] <= P_NEW_MAX_SEED))
    zl["meets_criterion1"] = ok.astype(int)
    met = bool(ok.any())
    zl["provisional"] = 0 if met else 1
    zl["note"] = "" if met else "criterion 1 not met (see report.md)"
    if met and not bool(ok.iloc[0]):
        zl = pd.concat([zl[ok].head(1), zl.drop(zl[ok].head(1).index)])
        zl["rank"] = np.arange(1, len(zl) + 1)
    cols = (["rank"] + ["config_id"] + list(AXES) + ["loss", "n_missed", "provisional", "meets_criterion1", "note",
            "phaseA_rank", "phaseA_loss", "n_seeds"] + comps + list(TARGETS) + EXTRA + NONTARGET + CHECKS)
    zl = zl[cols]
    # 001 p_new = 0.5 point: upper bound only (bigram, re-samples the seed vocabulary), not a candidate
    diag = load_runs(REF001 / "diag_full" / "ZL_split")
    dm = add_loss(diag[list(TARGETS) + EXTRA + NONTARGET].mean().to_frame().T, real_point("ZL_split"))
    up = {c: np.nan for c in cols}
    up.update({k: dm.iloc[0][k] for k in dm.columns if k in cols})
    up.update({"rank": np.nan, "config_id": f"001:{DIAG001['base_config']}", "edit_model": "bigram",
               "new_token_source": "seed", "p_new": 0.5, "window": 32, "p_copy_exact": 0.5, "p_two_edits": 0.5,
               "recency": 0.0, "provisional": 1, "meets_criterion1": 0, "n_seeds": len(diag),
               "note": "upper bound only: 001 bigram, p_new=0.5 > 0.3 (not a generative explanation)"})
    zl = pd.concat([zl, pd.DataFrame([up])], ignore_index=True)
    zl.to_csv(OUT / "best_configs.csv", index=False, float_format="%.6g")
    info["choice_on_ZL_split"] = choice
    info["criterion1_met"] = met
    (OUT / "summary_info.json").write_text(json.dumps(info, indent=2))
    plot_profile(profile)
    print(zl.head(12).to_string(float_format=lambda x: f"{x:.4g}"))
    print(json.dumps(info, indent=2))


def plot_profile(profile: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; profile.png skipped")
        return
    profile = profile[profile["kind"] != "check"]
    metrics = list(dict.fromkeys(profile["metric"]))
    datasets = list(dict.fromkeys(profile["dataset"]))
    colors = {"real": "#222222", "trigram_best": "#d0622c", "slot_best": "#2c8d5a", "bigram_001": "#8a8a8a",
              "bigram_001_pnew05_upper": "#b9a0cc", "shuffle_lines_001": "#3a78b5",
              "ref_RT2026_voynich": "#000000", "ref_RT2026_selfcitation": "#c9302c"}
    markers = {"ref_RT2026_voynich": "*", "ref_RT2026_selfcitation": "s"}
    ncol = 4
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.8 * nrow))
    for ax, m in zip(axes.flat, metrics):
        sub = profile[profile["metric"] == m]
        kind = sub["kind"].iloc[0]
        for di, ds in enumerate(datasets):
            d = sub[(sub["dataset"] == ds) & sub["mean"].notna()]
            for si, (_, r) in enumerate(d.iterrows()):
                y = di + 0.11 * si - 0.38
                lo, hi = r["lo"], r["hi"]
                xerr = None if pd.isna(lo) else [[max(r["mean"] - lo, 0)], [max(hi - r["mean"], 0)]]
                ax.errorbar(r["mean"], y, xerr=xerr, fmt=markers.get(r["source"], "o"), ms=4,
                            color=colors.get(r["source"], "grey"), label=r["source"], capsize=2)
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets, fontsize=8)
        ax.set_title(f"{m} ({kind})", fontsize=9, color="#000000" if kind == "target" else "#7a2a00")
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.3)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    handles, labels = [], []
    for ax in axes.flat:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    fig.legend(handles, labels, loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("003: real text vs grammar-bound auto-copy (best trigram / slot config, chosen on ZL split), "
                 "001 bigram generator, within-line shuffle, literature values (95 % intervals)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "profile.png", dpi=110)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["real", "grid", "rankA", "phaseB", "bigramgram", "diag", "summarize"])
    ap.add_argument("--file", default="ZL3b-n.txt")
    ap.add_argument("--space-policy", default="split", choices=["split", "join", "strict"])
    ap.add_argument("--seeds", default=None, help="encoded seeds config_id*100+rep: 'a:b', 'a:b/chunks:k' or list")
    ap.add_argument("--reps", type=int, default=None, help="seeds per config (A: 3, B: 50)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.command == "grid":
        args.seeds = args.seeds or f"0:{len(CONFIGS) * 100}"
        args.reps = args.reps or 3
        cmd_grid(args)
    elif args.command == "phaseB":
        args.reps = args.reps or 50
        cmd_phase_b(args)
    elif args.command == "bigramgram":
        args.reps = args.reps or 5
        cmd_bigramgram(args)
    elif args.command == "diag":
        args.reps = args.reps or 3
        cmd_diag(args)
    elif args.command == "real":
        cmd_real(args)
    elif args.command == "rankA":
        cmd_rank_a(args)
    else:
        cmd_summarize(args)


if __name__ == "__main__":
    main()
