"""Experiment 001: calibrate the auto-copy generator on ZL3b and list what it leaves unexplained.

Usage (from the repository root; see ``spec.md``)::

    # real text: point estimates (seed 0) + page bootstrap
    python experiments/001_autocopy_calibration/run.py real --file ZL3b-n.txt --space-policy split
    # phase A: 576 configs x 3 seeds, target metrics + edge/token MI; 4 parallel chunks
    python experiments/001_autocopy_calibration/run.py grid --seeds 0:57600/4:0 --out results/001_autocopy_calibration/phaseA/part0
    # phase B: 10 best configs x 50 seeds, all metrics (after `summarize` wrote phaseA_top10.csv)
    python experiments/001_autocopy_calibration/run.py phaseB --file ZL3b-n.txt --space-policy split
    # pure within-line shuffle, all metrics
    python experiments/001_autocopy_calibration/run.py shuffle --file ZL3b-n.txt --space-policy split
    # outside the grid: p_new / seed vocabulary / no-edit diagnostic for the best config
    python experiments/001_autocopy_calibration/run.py diag
    # tables, profile, plot
    python experiments/001_autocopy_calibration/run.py summarize

Run identification: ``seed = config_id * 100 + rep`` with ``config_id`` 0..575 in the order of
``itertools.product`` over the grid axes (see :data:`AXES`).  A ``--seeds`` range is read as a
range of encoded seeds; only seeds with ``rep < --reps`` are run, so ``0:57600/4:k`` gives chunk
``k`` of the configs.  Every command is resumable (``run_seeds`` skips finished seeds).
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
from voytools.nulls import autocopy_like, shuffle_within_lines  # noqa: E402
from voytools.stats import (  # noqa: E402
    adjacent_token_mi,
    edge_glyph_mi,
    glyph_entropy_report,
    neighbor_similarity,
    positional_profile,
    vocabulary_summary,
)

DATA = ROOT / "data"
OUT = ROOT / "results" / "001_autocopy_calibration"

AXES = {
    "window": (4, 8, 16, 32),
    "p_copy_exact": (0.1, 0.2, 0.35, 0.5),
    "p_new": (0.02, 0.05, 0.15),
    "p_two_edits": (0.1, 0.25, 0.5),
    "recency": (0.0, 0.3),
    "edit_type_probs": ((0.6, 0.2, 0.2), (0.4, 0.3, 0.3)),
}
CONFIGS = list(itertools.product(*AXES.values()))  # 576
TARGETS = {"share_near": 2.0, "hapax_share": 2.0, "H2": 2.0, "n_types": 1.0, "mean_len": 1.0}
COMPONENT_LIMIT = 0.05  # ((x - t) / t)^2 > 0.05, i.e. relative deviation > 22 %, is "missed"
NONTARGET = ["edge_mi_excess", "token_mi_excess_share", "bpe_gap64", "pos_jsd_first_mid",
             "cmi_excess", "frac_stems_multi_suffix", "mean_suffixes_per_stem"]
EXTRA = ["near_ratio"]  # reported, not in the loss
BOOT_METRICS = ["share_near", "near_ratio", "hapax_share", "n_types", "mean_len", "H2"]
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}


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
    return dict(zip(AXES, CONFIGS[cid]))


def config_columns(cid: int) -> dict:
    p = config_params(cid)
    p["edit_type_probs"] = "/".join(f"{x:g}" for x in p["edit_type_probs"])
    return {"config_id": cid, **p}


def tag_of(file: str, policy: str) -> str:
    return f"{'GC' if file.startswith('GC') else 'ZL'}_{policy}"


# --------------------------------------------------------------------------- #
# data and metrics
# --------------------------------------------------------------------------- #


class Context:
    """Real text of one file/space policy plus everything the metrics need to be comparable."""

    def __init__(self, file: str, policy: str):
        self.file, self.policy = file, policy
        self.real = filter_tokens(load_tables(DATA / file, space_policy=policy).tokens, generic_type="P", clean=True)
        self.eva = self.real.attrs.get("glyph_preset") is not None
        # 'extended' Eva composites for H2 (spec); v101 glyphs are already single characters and
        # the Eva composites would wrongly merge some of them, so GC uses the plain characters.
        self.h2_preset = "extended" if self.eva else None
        self.allowed_from = self.real["token"].tolist()
        d = decompose_table(self.real)  # affix lists learned on the real text, reused for every sample
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]

    def generate(self, cid: int, seed: int) -> pd.DataFrame:
        return autocopy_like(self.real, allowed_from=self.allowed_from, seed=seed, **config_params(cid))


def lines_of(df: pd.DataFrame) -> list[list[str]]:
    return df.groupby("locus_id", sort=False)["token"].apply(list).tolist()


def stats_for(df: pd.DataFrame, ctx: Context, seed: int, full: bool) -> dict:
    """All metrics of the spec for one table (real, bootstrap, shuffle or generated), identical code."""
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
    if full is None:  # bootstrap: target metrics only
        return out
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


def page_bootstrap(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Pages resampled with replacement; duplicated pages get distinct folio labels."""
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


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_real(args) -> None:
    ctx = Context(args.file, args.space_policy)
    out = OUT / "real" / tag_of(args.file, args.space_policy)
    out.mkdir(parents=True, exist_ok=True)
    point = {"file": args.file, "space_policy": args.space_policy, "source": "real",
             **stats_for(ctx.real, ctx, 0, True)}
    pd.DataFrame([point]).to_csv(out / "point.csv", index=False)
    print(pd.Series(point).to_string())
    run_seeds("001_real_bootstrap", range(args.boot), lambda s: stats_for(page_bootstrap(ctx.real, s), ctx, s, None),
              out / "bootstrap", params={"file": args.file, "space_policy": args.space_policy, "unit": "page"},
              progress=progress)


def _gen_runner(ctx: Context, full: bool):
    def one_run(seed: int) -> dict:
        cid = seed // 100
        res = stats_for(ctx.generate(cid, seed), ctx, seed, full)
        return {"file": ctx.file, "space_policy": ctx.policy, "source": "autocopy", **config_columns(cid), **res}

    return one_run


def cmd_grid(args) -> None:
    ctx = Context(args.file, args.space_policy)
    seeds = [s for s in parse_seeds(args.seeds) if s % 100 < args.reps and s // 100 < len(CONFIGS)]
    run_seeds("001_phaseA", seeds, _gen_runner(ctx, False), args.out or OUT / "phaseA" / "part0",
              params={"file": args.file, "space_policy": args.space_policy, "seeds": args.seeds,
                      "reps": args.reps, "metrics": "target+edge/token MI", "axes": {k: list(v) for k, v in AXES.items()}},
              progress=progress)


def cmd_phase_b(args) -> None:
    ctx = Context(args.file, args.space_policy)
    top = pd.read_csv(OUT / "phaseA_top10.csv")["config_id"].astype(int).tolist()
    seeds = [cid * 100 + r for cid in top for r in range(args.reps)]
    if args.seeds:
        keep = set(parse_seeds(args.seeds))
        seeds = [s for s in seeds if s in keep]
    out = args.out or OUT / "phaseB" / tag_of(args.file, args.space_policy)
    run_seeds("001_phaseB", seeds, _gen_runner(ctx, True), out,
              params={"file": args.file, "space_policy": args.space_policy, "configs": top, "reps": args.reps},
              progress=progress)


DIAG_BASE = 573  # best phase-A config; the diagnostic varies p_new, p_copy_exact and the seed vocabulary
DIAG = list(itertools.product(("first100", "all"), (0.15, 0.3, 0.5, 0.8), (0.5, 1.0)))


def cmd_diag(args) -> None:
    """Outside the spec grid: which knob brings H2 / n_types / share_near down together?"""
    ctx = Context(args.file, args.space_policy)

    def one_run(seed: int) -> dict:
        vocab, p_new, pce = DIAG[seed // 100]
        params = {**config_params(DIAG_BASE), "p_new": p_new, "p_copy_exact": pce}
        seeds = ctx.allowed_from if vocab == "all" else None
        df = autocopy_like(ctx.real, seed_tokens=seeds, allowed_from=ctx.allowed_from, seed=seed, **params)
        return {"file": ctx.file, "space_policy": ctx.policy, "base_config": DIAG_BASE, "seed_vocab": vocab,
                "p_new": p_new, "p_copy_exact": pce, **stats_for(df, ctx, seed, False)}

    seeds = [i * 100 + r for i in range(len(DIAG)) for r in range(args.reps)]
    run_seeds("001_diag", seeds, one_run, OUT / "diag", params={"base_config": DIAG_BASE, "reps": args.reps},
              progress=progress)
    runs = load_runs(OUT / "diag")
    means = runs.groupby(["seed_vocab", "p_new", "p_copy_exact"])[list(TARGETS) + EXTRA + ["edge_mi_excess", "token_mi_excess_share"]].mean()
    means = add_loss(means, real_point(tag_of(args.file, args.space_policy)))
    means.to_csv(OUT / "diag_summary.csv", float_format="%.6g")
    print(means.to_string(float_format=lambda x: f"{x:.4g}"))


def cmd_shuffle(args) -> None:
    ctx = Context(args.file, args.space_policy)

    def one_run(seed: int) -> dict:
        res = stats_for(shuffle_within_lines(ctx.real, seed=seed), ctx, seed, True)
        return {"file": ctx.file, "space_policy": ctx.policy, "source": "shuffle_lines", **res}

    run_seeds("001_shuffle", range(args.reps), one_run, OUT / "shuffle" / tag_of(args.file, args.space_policy),
              params={"file": args.file, "space_policy": args.space_policy}, progress=progress)


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
    g = runs.groupby(keys, sort=True)
    means = g[metrics].mean()
    means.insert(0, "n_seeds", g.size())
    return means.reset_index()


def suspicious(runs: pd.DataFrame) -> pd.DataFrame:
    bad = (runs["share_near"] < 0.05) | (runs["share_near"] > 0.9) | (runs["H2"] < 1.5) | (runs["H2"] > 4) | (runs["hapax_share"] < 0)
    return runs[bad]


def real_point(tag: str) -> pd.Series:
    return pd.read_csv(OUT / "real" / tag / "point.csv").iloc[0]


def cmd_summarize(args) -> None:
    parts = sorted((OUT / "phaseA").glob("part*"))
    runs_a = merge_runs(parts, OUT)
    metrics_a = list(TARGETS) + EXTRA + ["edge_mi_excess", "token_mi_excess_share"]
    t_zl = real_point("ZL_split")
    grid = add_loss(config_means(runs_a, metrics_a), t_zl).sort_values("loss", kind="stable").reset_index(drop=True)
    grid.insert(0, "rank", np.arange(1, len(grid) + 1))
    grid.to_csv(OUT / "grid_summary.csv", index=False, float_format="%.6g")
    grid.head(10).to_csv(OUT / "phaseA_top10.csv", index=False, float_format="%.6g")
    sus = suspicious(runs_a)
    info = {"phaseA_runs": len(runs_a), "phaseA_configs": int(runs_a["config_id"].nunique()),
            "phaseA_suspicious_runs": len(sus)}
    print(grid.head(15).to_string(float_format=lambda x: f"{x:.4g}"))

    tags = [p.name for p in (OUT / "phaseB").glob("*_*") if (p / "runs.csv").exists()]
    tags = sorted(tags, key=lambda t: (t != "ZL_split", t))  # ZL split first: its best config is the reference
    if not tags:
        (OUT / "summary_info.json").write_text(json.dumps(info, indent=2))
        return
    metrics_b = list(TARGETS) + EXTRA + NONTARGET
    profile_rows, best_rows, summary_rows = [], [], []
    for tag in tags:
        target = real_point(tag)
        runs_b = load_runs(OUT / "phaseB" / tag)
        info[f"phaseB_{tag}_runs"] = len(runs_b)
        info[f"phaseB_{tag}_suspicious_runs"] = len(suspicious(runs_b))
        mb = add_loss(config_means(runs_b, metrics_b), target).sort_values("loss", kind="stable").reset_index(drop=True)
        mb.insert(0, "dataset", tag)
        mb.insert(1, "rank", np.arange(1, len(mb) + 1))
        best_rows.append(mb)
        best = int(mb.iloc[0]["config_id"])
        zl_best = int(best_rows[0].iloc[0]["config_id"])
        boot = load_runs(OUT / "real" / tag / "bootstrap")
        shuf = load_runs(OUT / "shuffle" / tag) if (OUT / "shuffle" / tag / "runs.csv").exists() else pd.DataFrame()
        for m in metrics_b:
            x = float(target[m])
            row = {"dataset": tag, "metric": m, "kind": "target" if m in TARGETS else ("extra" if m in EXTRA else "nontarget")}
            prof = [{**row, "source": "real", "config_id": "", "n": 1, "mean": x,
                     "lo": x - 1.96 * boot[m].std(ddof=1) if m in boot else np.nan,
                     "hi": x + 1.96 * boot[m].std(ddof=1) if m in boot else np.nan,
                     "sd": boot[m].std(ddof=1) if m in boot else np.nan,
                     "interval": "page bootstrap, normal (x +- 1.96 SD)" if m in boot else "",
                     "z_real": np.nan, "diff_real": np.nan}]
            sources = [("autocopy_best", runs_b[runs_b["config_id"] == best], best)]
            if best != zl_best:
                sources.append(("autocopy_zlbest", runs_b[runs_b["config_id"] == zl_best], zl_best))
            if not shuf.empty:
                sources.append(("shuffle_lines", shuf, ""))
            for name, r, cid in sources:
                v = r[m].astype(float)
                sd = v.std(ddof=1)
                prof.append({**row, "source": name, "config_id": cid, "n": len(v), "mean": v.mean(),
                             "lo": v.quantile(0.025), "hi": v.quantile(0.975), "sd": sd,
                             "interval": "2.5/97.5 % over seeds", "diff_real": x - v.mean(),
                             "z_real": (x - v.mean()) / sd if sd > 0 else np.nan})
            profile_rows += prof
        summ = runs_b[runs_b["config_id"] == best][metrics_b].describe(percentiles=[0.025, 0.5, 0.975]).T
        summ.insert(0, "dataset", tag)
        summ.insert(1, "config_id", best)
        summary_rows.append(summ.reset_index().rename(columns={"index": "metric"}))
    best_all = pd.concat(best_rows, ignore_index=True)
    params_cols = ["config_id"] + list(AXES)
    zl = best_all[best_all["dataset"] == "ZL_split"]
    pa = grid.set_index("config_id")[["rank", "loss"]].rename(columns={"rank": "phaseA_rank", "loss": "phaseA_loss"})
    best_cfg = zl.join(pa, on="config_id")
    best_cfg = best_cfg[["rank"] + params_cols + ["loss", "n_missed", "phaseA_rank", "phaseA_loss", "n_seeds"]
                        + [f"comp_{f}" for f in TARGETS] + list(TARGETS) + EXTRA + NONTARGET]
    best_cfg.to_csv(OUT / "best_configs.csv", index=False, float_format="%.6g")
    best_all.to_csv(OUT / "phaseB_configs.csv", index=False, float_format="%.6g")
    profile = pd.DataFrame(profile_rows)
    profile.to_csv(OUT / "profile.csv", index=False, float_format="%.6g")
    pd.concat(summary_rows, ignore_index=True).to_csv(OUT / "summary.csv", index=False, float_format="%.6g")
    (OUT / "summary_info.json").write_text(json.dumps(info, indent=2))
    plot_profile(profile)
    print(profile.to_string(float_format=lambda x: f"{x:.4g}"))


def plot_profile(profile: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; profile.png skipped")
        return
    metrics = list(dict.fromkeys(profile["metric"]))
    datasets = list(dict.fromkeys(profile["dataset"]))
    colors = {"real": "#222222", "autocopy_best": "#d0622c", "autocopy_zlbest": "#e8a33d", "shuffle_lines": "#3a78b5"}
    ncol = 4
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.6 * nrow))
    for ax, m in zip(axes.flat, metrics):
        sub = profile[profile["metric"] == m]
        kind = sub["kind"].iloc[0]
        for di, ds in enumerate(datasets):
            d = sub[sub["dataset"] == ds]
            for si, (_, r) in enumerate(d.iterrows()):
                y = di + 0.18 * si - 0.25
                lo, hi = r["lo"], r["hi"]
                xerr = None if pd.isna(lo) else [[r["mean"] - lo], [hi - r["mean"]]]
                ax.errorbar(r["mean"], y, xerr=xerr, fmt="o", ms=4, color=colors.get(r["source"], "grey"),
                            label=r["source"] if (di == 0) else None, capsize=2)
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets, fontsize=8)
        ax.set_title(f"{m} ({kind})", fontsize=9, color="#000000" if kind == "target" else "#7a2a00")
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.3)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("001: real text vs best auto-copy config vs within-line shuffle (95 % intervals)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "profile.png", dpi=110)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["real", "grid", "phaseB", "shuffle", "diag", "summarize"])
    ap.add_argument("--file", default="ZL3b-n.txt")
    ap.add_argument("--space-policy", default="split", choices=["split", "join", "strict"])
    ap.add_argument("--seeds", default=None, help="encoded seeds config_id*100+rep: 'a:b', 'a:b/chunks:k' or list")
    ap.add_argument("--reps", type=int, default=None, help="seeds per config (A: 3, B: 50) or shuffles")
    ap.add_argument("--boot", type=int, default=200, help="page bootstrap replicates for the real text")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.command == "grid":
        args.seeds = args.seeds or f"0:{len(CONFIGS) * 100}"
        args.reps = args.reps or 3
        cmd_grid(args)
    elif args.command == "phaseB":
        args.reps = args.reps or 50
        cmd_phase_b(args)
    elif args.command == "shuffle":
        args.reps = args.reps or 50
        cmd_shuffle(args)
    elif args.command == "diag":
        args.reps = args.reps or 3
        cmd_diag(args)
    elif args.command == "real":
        cmd_real(args)
    else:
        cmd_summarize(args)


if __name__ == "__main__":
    main()
