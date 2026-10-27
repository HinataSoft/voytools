"""Experiment 002: do stems and affixes behave like a paradigm or like a slot template?

Usage (from the repository root; see ``spec.md``)::

    # one combination (file x space policy x decomposition), all sources, default seeds
    python experiments/002_paradigm_coherence/run.py run --combo ZL_split_learned
    # only some sources / seeds (resumable; finished seeds are skipped)
    python experiments/002_paradigm_coherence/run.py run --combo GC_join_learned --sources real_boot --seeds 0:100
    # list the combinations
    python experiments/002_paradigm_coherence/run.py combos
    # merge parts, summary.csv, effects.csv, inventories.csv, paradigm.png
    python experiments/002_paradigm_coherence/run.py summarize

Combinations: ZL3b x {split, join, strict} x {learned, reference} affix lists and GC2a x {split, join,
strict} x learned (9).  Sources per combination (``seed`` = repetition number inside a source):

* ``real`` – the real text, seed 0 (point estimate);
* ``real_boot`` – page bootstrap of the real text, seeds 0..199;
* ``slot_indep`` – prefix, stem and suffix drawn independently from the real marginals, 0..49;
* ``autocopy_573`` / ``autocopy_536`` – auto-copy generator with rows 1 and 2 of
  ``results/001_autocopy_calibration/best_configs.csv``, 0..49 each;
* ``shuffle_lines`` – within-line shuffle, 0..4.

The stem-conditional shuffle (null 2 of the spec) is ``cmi_shuffled`` inside every run.  Affix lists
are learned once on the real text of the combination (or the reference lists) and reused for every
sample.  Runs go to ``results/002_paradigm_coherence/parts/<combo>/<source>/``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_runs, load_tables, run_seeds  # noqa: E402
from voytools.ivtff import UNCERTAIN_SPACE_MARK  # noqa: E402
from voytools.morphology import (  # noqa: E402
    decompose_table,
    inventory_jaccard,
    paradigm_metrics,
    productivity_slope,
    slot_independent_sample,
    stem_dispersion,
)
from voytools.nulls import autocopy_like, shuffle_within_lines  # noqa: E402

DATA = ROOT / "data"
OUT = ROOT / "results" / "002_paradigm_coherence"
BEST = ROOT / "results" / "001_autocopy_calibration" / "best_configs.csv"

FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
POLICIES = ("split", "join", "strict")
COMBOS = [f"ZL_{p}_{d}" for p in POLICIES for d in ("learned", "reference")] + [f"GC_{p}_learned" for p in POLICIES]
GEN_PARAMS = ("window", "p_copy_exact", "p_new", "p_two_edits", "recency", "edit_type_probs")
MIN_STEM = 5  # paradigm_metrics, productivity_slope, stem_dispersion
MIN_INV = 3  # inventory threshold for the Jaccard index


def autocopy_configs() -> dict[str, dict]:
    """Rows 1 and 2 of 001's best_configs.csv as generator keyword arguments."""
    best = pd.read_csv(BEST).sort_values("rank").head(2)
    out = {}
    for _, r in best.iterrows():
        p = {k: r[k] for k in GEN_PARAMS}
        p["window"] = int(p["window"])
        p["edit_type_probs"] = tuple(float(x) for x in str(p["edit_type_probs"]).split("/"))
        for k in ("p_copy_exact", "p_new", "p_two_edits", "recency"):
            p[k] = float(p[k])
        out[f"autocopy_{int(r['config_id'])}"] = p
    return out


ACFG = autocopy_configs()
SOURCES = {"real": range(1), "real_boot": range(200), "slot_indep": range(50),
           **{name: range(50) for name in ACFG}, "shuffle_lines": range(5)}


def parse_seeds(spec: str) -> list[int]:
    """``'0:20'`` -> seeds 0..19; ``'3,5,8'`` -> those seeds."""
    if ":" in spec:
        a, b = spec.split(":")
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(",") if s.strip()]


# --------------------------------------------------------------------------- #
# data and metrics
# --------------------------------------------------------------------------- #


class Context:
    """Real text of one combination plus the affix lists every sample is decomposed with."""

    def __init__(self, combo: str):
        self.combo = combo
        key, self.policy, self.decomposition = combo.split("_")
        self.file = FILES[key]
        real = filter_tokens(load_tables(DATA / self.file, space_policy=self.policy).tokens, generic_type="P", clean=True)
        comma = real["token"].str.contains(UNCERTAIN_SPACE_MARK, regex=False)
        self.n_dropped_comma = int(comma.sum())  # strict: tokens with an uncertain space inside are dropped
        self.real = real[~comma].reset_index(drop=True) if self.n_dropped_comma else real
        self.real.attrs.update(real.attrs)
        spec = "reference" if self.decomposition == "reference" else None
        d = decompose_table(self.real, prefixes=spec, suffixes=spec)
        self.real_d = d
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]
        self.allowed_from = self.real["token"].tolist()

    def columns(self) -> dict:
        return {"file": self.file, "space_policy": self.policy, "decomposition": self.decomposition}

    def sample(self, source: str, seed: int) -> pd.DataFrame:
        if source == "real":
            return self.real
        if source == "real_boot":
            return page_bootstrap(self.real, seed)
        if source == "slot_indep":
            return slot_independent_sample(self.real_d, seed=seed)
        if source == "shuffle_lines":
            return shuffle_within_lines(self.real, seed=seed)
        if source in ACFG:
            return autocopy_like(self.real, allowed_from=self.allowed_from, seed=seed, **ACFG[source])
        raise ValueError(f"unknown source {source!r}")


def page_bootstrap(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Pages resampled with replacement; duplicated pages get distinct folio labels (as in 001)."""
    rng = np.random.default_rng(seed)
    folios = df["folio"].drop_duplicates().tolist()
    groups = {f: g for f, g in df.groupby("folio", sort=False)}
    parts = []
    for i, f in enumerate(rng.choice(len(folios), size=len(folios), replace=True)):
        g = groups[folios[f]].copy()
        g["folio"] = f"{folios[f]}#{i}"
        parts.append(g)
    out = pd.concat(parts, ignore_index=True)
    out.attrs.update(df.attrs)
    return out


def _mean_pairs(dfd: pd.DataFrame, by: str, column: str, groups=None) -> tuple[float, float]:
    t = inventory_jaccard(dfd, by, column=column, min_count=MIN_INV, groups=groups)
    return float(t["jaccard"].mean()), float(t["jsd"].mean())


def metrics(df: pd.DataFrame, ctx: Context, seed: int) -> dict:
    """All metrics of the spec for one table (real, bootstrap or null sample), identical code."""
    dfd = decompose_table(df, prefixes=ctx.prefixes, suffixes=ctx.suffixes)
    pm = paradigm_metrics(dfd, min_stem_count=MIN_STEM, seed=seed)
    s = pm.summary
    low = pm.per_stem[pm.per_stem["count"] < 2 * MIN_STEM]
    out = {k: float(s[k]) for k in ("n_tokens", "n_stems", "n_stems_k", "token_coverage_k", "mean_suffixes_per_stem",
                                     "frac_stems_multi_suffix", "H_suffix_given_stem", "mi_prefix_suffix",
                                     "cmi_prefix_suffix_given_stem", "cmi_shuffled", "cmi_excess")}
    # frequency-controlled: stems with 5..9 tokens only
    out["frac_multi_suffix_5to9"] = float((low["n_suffixes"] > 1).mean()) if len(low) else float("nan")
    out["jaccard_suffix_sections"], out["jsd_suffix_sections"] = _mean_pairs(dfd, "illustration", "suffix")
    out["jaccard_stem_sections"], out["jsd_stem_sections"] = _mean_pairs(dfd, "illustration", "stem")
    out["jaccard_gap_sections"] = out["jaccard_suffix_sections"] - out["jaccard_stem_sections"]
    out["jaccard_suffix_AB"], out["jsd_suffix_AB"] = _mean_pairs(dfd, "language", "suffix", ["A", "B"])
    out["jaccard_stem_AB"], out["jsd_stem_AB"] = _mean_pairs(dfd, "language", "stem", ["A", "B"])
    out["jaccard_gap_AB"] = out["jaccard_suffix_AB"] - out["jaccard_stem_AB"]
    ps = productivity_slope(dfd, min_stem_count=MIN_STEM)
    out["productivity_slope"], out["productivity_intercept"], out["productivity_r2"] = ps["slope"], ps["intercept"], ps["r2"]
    out["stem_dispersion"] = stem_dispersion(dfd, min_stem_count=MIN_STEM)
    return out


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def progress(i: int, n: int, s: int) -> None:
    if i % 25 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_run(args) -> None:
    ctx = Context(args.combo)
    sources = args.sources.split(",") if args.sources else list(SOURCES)
    print(f"{args.combo}: {len(ctx.real)} tokens, dropped (comma) {ctx.n_dropped_comma}, "
          f"prefixes {ctx.prefixes}, suffixes {ctx.suffixes}", flush=True)
    for source in sources:
        seeds = list(SOURCES[source])
        if args.seeds:
            keep = set(parse_seeds(args.seeds))
            seeds = [s for s in seeds if s in keep]
        print(f"[{args.combo}] {source}: {len(seeds)} seeds", flush=True)

        def one_run(seed: int, source=source) -> dict:
            return {"source": source, **ctx.columns(), **metrics(ctx.sample(source, seed), ctx, seed)}

        params = {"combo": args.combo, **ctx.columns(), "source": source, "n_tokens_real": len(ctx.real),
                  "n_dropped_comma": ctx.n_dropped_comma, "prefixes": ctx.prefixes, "suffixes": ctx.suffixes,
                  "min_stem_count": MIN_STEM, "min_inventory_count": MIN_INV,
                  "generator": {k: str(v) for k, v in ACFG.get(source, {}).items()}}
        run_seeds("002_paradigm_coherence", seeds, one_run, OUT / "parts" / args.combo / source,
                  params=params, progress=progress)


def cmd_combos(args) -> None:
    print("\n".join(COMBOS))


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #

ID_COLS = ["combo", "source", "file", "space_policy", "decomposition", "seed"]
NULLS = ["slot_indep", *ACFG, "shuffle_lines"]
# metrics that depend only on the bag of tokens (identical under a within-line shuffle)
TYPE_LEVEL = ["n_stems_k", "mean_suffixes_per_stem", "frac_stems_multi_suffix", "H_suffix_given_stem",
              "mi_prefix_suffix", "cmi_prefix_suffix_given_stem", "jaccard_suffix_sections",
              "jaccard_stem_sections", "jaccard_suffix_AB", "jaccard_stem_AB", "productivity_slope"]
PLOT = ["cmi_excess", "frac_stems_multi_suffix", "jaccard_stem_sections", "jaccard_gap_sections", "stem_dispersion"]
COLORS = {"real": "#222222", "slot_indep": "#4c78a8", "autocopy_573": "#e45756", "autocopy_536": "#f58518",
          "shuffle_lines": "#72b7b2"}


def load_all() -> pd.DataFrame:
    frames = []
    for combo in COMBOS:
        for source in SOURCES:
            r = load_runs(OUT / "parts" / combo / source)
            if not r.empty:
                r.insert(0, "combo", combo)
                frames.append(r)
    runs = pd.concat(frames, ignore_index=True)
    return runs[ID_COLS + [c for c in runs.columns if c not in ID_COLS and c != "elapsed_s"] + ["elapsed_s"]]


def metric_columns(runs: pd.DataFrame) -> list[str]:
    return [c for c in runs.columns if c not in ID_COLS and c != "elapsed_s"]


def cmd_summarize(args) -> None:
    runs = load_all()
    runs.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    mcols = metric_columns(runs)
    keys = ["combo", "file", "space_policy", "decomposition", "source"]

    long = runs.melt(id_vars=keys + ["seed"], value_vars=mcols, var_name="metric")
    g = long.groupby(keys + ["metric"], sort=False)["value"]
    summary = pd.concat({"n": g.count(), "mean": g.mean(), "sd": g.std(ddof=1),
                         "q025": g.quantile(0.025), "q975": g.quantile(0.975)}, axis=1).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False, float_format="%.6g")

    # effects: real point estimate minus null mean; interval from bootstrap SD and null SD
    rows = []
    for combo, sub in runs.groupby("combo", sort=False):
        point = sub[sub["source"] == "real"].iloc[0]
        boot = sub[sub["source"] == "real_boot"]
        for m in mcols:
            sd_r = float(boot[m].std(ddof=1))
            for null in NULLS:
                v = sub.loc[sub["source"] == null, m].astype(float)
                if v.empty:
                    continue
                diff = float(point[m]) - v.mean()
                sd_n = float(v.std(ddof=1))
                sd = float(np.sqrt(np.nansum([sd_r**2, sd_n**2])))
                rows.append({"combo": combo, "metric": m, "null": null, "real": float(point[m]), "real_boot_sd": sd_r,
                             "null_mean": v.mean(), "null_sd": sd_n, "null_q025": v.quantile(0.025),
                             "null_q975": v.quantile(0.975), "diff": diff, "ci_lo": diff - 1.96 * sd,
                             "ci_hi": diff + 1.96 * sd, "z_null": diff / sd_n if sd_n > 0 else float("nan"),
                             "n_null": len(v)})
        # null 2: stem-conditional shuffle of the prefix labels (inside paradigm_metrics)
        sd_r = float(boot["cmi_excess"].std(ddof=1))
        diff = float(point["cmi_excess"])
        rows.append({"combo": combo, "metric": "cmi_prefix_suffix_given_stem", "null": "stem_cond_shuffle",
                     "real": float(point["cmi_prefix_suffix_given_stem"]), "real_boot_sd": float(boot["cmi_prefix_suffix_given_stem"].std(ddof=1)),
                     "null_mean": float(point["cmi_shuffled"]), "null_sd": float("nan"), "null_q025": float("nan"),
                     "null_q975": float("nan"), "diff": diff, "ci_lo": diff - 1.96 * sd_r, "ci_hi": diff + 1.96 * sd_r,
                     "z_null": diff / sd_r if sd_r > 0 else float("nan"), "n_null": 20})
    effects = pd.DataFrame(rows)
    effects.insert(1, "file", effects["combo"].map(lambda c: FILES[c.split("_")[0]]))
    effects.insert(2, "space_policy", effects["combo"].map(lambda c: c.split("_")[1]))
    effects.insert(3, "decomposition", effects["combo"].map(lambda c: c.split("_")[2]))
    effects.to_csv(OUT / "effects.csv", index=False, float_format="%.6g")

    # inventories of the real text, every pair of sections and A vs B
    inv = []
    for combo in COMBOS:
        ctx = Context(combo)
        dfd = ctx.real_d
        for by, groups in (("illustration", None), ("language", ["A", "B"])):
            for column in ("suffix", "stem"):
                t = inventory_jaccard(dfd, by, column=column, min_count=MIN_INV, groups=groups)
                t.insert(0, "combo", combo)
                t.insert(1, "by", by)
                t.insert(2, "column", column)
                inv.append(t)
    pd.concat(inv, ignore_index=True).to_csv(OUT / "inventories.csv", index=False, float_format="%.6g")

    info = checks(runs, effects)
    (OUT / "summary_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(info, indent=2, ensure_ascii=False))
    plot(summary)


def checks(runs: pd.DataFrame, effects: pd.DataFrame) -> dict:
    """Completeness, suspicious values and the shuffle-invariance check."""
    counts = runs.groupby(["combo", "source"]).size().unstack(fill_value=0)
    expected = {s: len(r) for s, r in SOURCES.items()}
    complete = {c: bool(all(counts.loc[c].get(s, 0) == n for s, n in expected.items())) if c in counts.index else False
                for c in COMBOS}
    sus = []
    for _, r in runs.iterrows():
        tag = f"{r['combo']}/{r['source']}/{r['seed']}"
        if r["cmi_excess"] > 1:
            sus.append(f"{tag}: cmi_excess {r['cmi_excess']:.3f} > 1")
        for m in ("jaccard_suffix_sections", "jaccard_stem_sections", "jaccard_suffix_AB", "jaccard_stem_AB"):
            if r[m] in (0.0, 1.0):
                sus.append(f"{tag}: {m} = {r[m]:g}")
        if r["stem_dispersion"] > 10:
            sus.append(f"{tag}: stem_dispersion {r['stem_dispersion']:.2f} > 10")
    real = runs[runs["source"] == "real"].set_index("combo")
    sh = runs[runs["source"] == "shuffle_lines"]
    inv = {}
    for m in TYPE_LEVEL + ["stem_dispersion"]:
        d = (sh.set_index("combo")[m] - real[m].reindex(sh["combo"]).to_numpy()).abs()
        inv[m] = float(d.max()) if len(d) else float("nan")
    return {"complete": complete, "runs_per_combo_source": counts.to_dict(orient="index"),
            "n_suspicious": len(sus), "suspicious_examples": sus[:40],
            "max_abs_diff_shuffle_vs_real": inv}


def plot(summary: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; paradigm.png skipped")
        return
    src = ["real", *NULLS]
    fig, axes = plt.subplots(len(PLOT), 1, figsize=(11, 2.6 * len(PLOT)), sharex=True)
    x = np.arange(len(COMBOS))
    width = 0.14
    for ax, m in zip(axes, PLOT):
        for j, s in enumerate(src):
            sub = summary[(summary["metric"] == m)]
            if s == "real":  # point estimate with the bootstrap interval
                pt = sub[sub["source"] == "real"].set_index("combo")["mean"].reindex(COMBOS)
                bt = sub[sub["source"] == "real_boot"].set_index("combo").reindex(COMBOS)
                lo, hi = bt["q025"], bt["q975"]
                y = pt
            else:
                t = sub[sub["source"] == s].set_index("combo").reindex(COMBOS)
                y, lo, hi = t["mean"], t["q025"], t["q975"]
            xx = x + (j - (len(src) - 1) / 2) * width
            ax.errorbar(xx, y, yerr=[(y - lo).clip(lower=0), (hi - y).clip(lower=0)], fmt="o", ms=4, capsize=2,
                        color=COLORS.get(s, "grey"), label=s)
        ax.set_ylabel(m, fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    axes[0].legend(ncol=len(src), fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.35))
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(COMBOS, rotation=30, ha="right", fontsize=8)
    fig.suptitle("002: real (bootstrap 95 %) vs null models (2.5-97.5 %)", y=0.995, fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "paradigm.png", dpi=110)
    print(f"wrote {OUT / 'paradigm.png'}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="all (or some) sources of one combination")
    p.add_argument("--combo", required=True, choices=COMBOS)
    p.add_argument("--sources", default="", help="comma-separated subset of " + ",".join(SOURCES))
    p.add_argument("--seeds", default="", help="restrict seeds, e.g. 0:100 or 3,5,8")
    p.set_defaults(fn=cmd_run)
    sub.add_parser("combos").set_defaults(fn=cmd_combos)
    sub.add_parser("summarize").set_defaults(fn=cmd_summarize)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
