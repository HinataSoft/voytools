"""Experiment 011: tables and plot from ``results/011_suffix_variation/runs/*``.

Writes ``runs.csv`` (all rows), ``decomposition.csv``, ``residual_structure.csv``, ``consistency.csv``,
``decomposition.png`` and ``meta.json``.  Intervals as in 007/009: real text and controls = point
estimate + 2 x the 2.5 / 97.5 % quantiles of the subsample deviations (80 % of pages, scale
sqrt(0.8 / 0.2) = 2); nulls = mean and 2.5 / 97.5 % quantiles over seeds.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from voytools import load_runs  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402

import run as R  # noqa: E402

OUT = R.OUT
SCALE = float(np.sqrt(R.M9.SUB_FRAC / (1 - R.M9.SUB_FRAC)))  # = 2
Q = (0.025, 0.975)


def load_all() -> pd.DataFrame:
    frames = [load_runs(d) for d in sorted((OUT / "runs").iterdir()) if (d / "runs.csv").exists()]
    runs = pd.concat(frames, ignore_index=True)
    runs["kind"] = runs["kind"].fillna("null")  # read_csv parses the label 'null' as NA
    return runs


def stat(g: pd.DataFrame, col: str) -> dict:
    """Point/interval of one column for one source (real/control: point + subsamples; null: seeds)."""
    if col not in g.columns:
        return {"value": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
    point = g[g["kind"] == "point"]
    if len(point):
        v = float(point[col].iloc[0])
        s = g.loc[g["kind"] == "sub", col].dropna()
        if len(s) < 2:
            return {"value": v, "lo": np.nan, "hi": np.nan, "n": len(s)}
        dev = s - s.mean()
        return {"value": v, "lo": v + SCALE * float(dev.quantile(Q[0])), "hi": v + SCALE * float(dev.quantile(Q[1])),
                "n": len(s)}
    s = g.loc[g["kind"] == "null", col].dropna()
    if not len(s):
        return {"value": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
    return {"value": float(s.mean()), "lo": float(s.quantile(Q[0])), "hi": float(s.quantile(Q[1])), "n": len(s)}


def point_or_mean(g: pd.DataFrame, col: str) -> float:
    if col not in g.columns:
        return np.nan
    point = g[g["kind"] == "point"]
    return float(point[col].iloc[0]) if len(point) else float(g[col].mean())


def source_class(src: str) -> str:
    if src in R.REAL:
        return "voynich"
    return "control" if src.startswith("ctl:") else "null"


def decomposition_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src, g in runs.groupby("source", sort=False):
        for var in ("a", "b"):
            if f"{var}__share_empty_stem" not in g.columns or g[f"{var}__share_empty_stem"].isna().all():
                continue
            for target in R.TARGETS:
                for x in R.BY:
                    pre = f"{var}__{target}__{x}__"
                    e = stat(g, pre + "explained")
                    rows.append({
                        "source": src, "class": source_class(src), "variant": var, "target": target, "X": x,
                        "n_stems": point_or_mean(g, f"{var}__{target}__n_stems"),
                        "n_tokens": point_or_mean(g, f"{var}__{target}__n_tok"),
                        "H_given": point_or_mean(g, pre + "H_given"), "H_cond": point_or_mean(g, pre + "H_cond"),
                        "bias": point_or_mean(g, pre + "bias"), "reduction": point_or_mean(g, pre + "reduction"),
                        "explained": e["value"], "lo": e["lo"], "hi": e["hi"], "n_rep": e["n"],
                        "explained_raw": point_or_mean(g, pre + "explained_raw"),
                        "explained_cv": point_or_mean(g, pre + "explained_cv"),
                        "explained_min5": point_or_mean(g, pre + "explained_min5"),
                        "H_given_min5": point_or_mean(g, pre + "H_given_min5"),
                        "z_point": point_or_mean(g, pre + "z"),
                    })
    return pd.DataFrame(rows)


def residual_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src, g in runs.groupby("source", sort=False):
        for var in ("a", "b"):
            if f"{var}__share_empty_stem" not in g.columns or g[f"{var}__share_empty_stem"].isna().all():
                continue
            for name in ("nb_prev_suffix", "nb_next_prefix"):
                pre = f"{var}__{name}__"
                e = stat(g, pre + "mi")
                rows.append({"source": src, "class": source_class(src), "variant": var, "metric": name,
                             "raw": point_or_mean(g, pre + "mi_raw"), "perm_mean": point_or_mean(g, pre + "mi_perm"),
                             "value": e["value"], "lo": e["lo"], "hi": e["hi"],
                             "perm_sd": point_or_mean(g, pre + "perm_sd"), "p": point_or_mean(g, pre + "p"),
                             "share_of_H": np.nan, "n": point_or_mean(g, pre + "n_pairs")})
            for x in R.RESID_BY:
                pre = f"{var}__resid_{x}__"
                e = stat(g, pre + "reduction")
                h = point_or_mean(g, pre + "H_given")
                rows.append({"source": src, "class": source_class(src), "variant": var,
                             "metric": f"I(suffix;{x}|stem,pos,lang,hand)",
                             "raw": point_or_mean(g, pre + "reduction_raw"), "perm_mean": point_or_mean(g, pre + "bias"),
                             "value": e["value"], "lo": e["lo"], "hi": e["hi"], "perm_sd": np.nan,
                             "p": np.nan, "share_of_H": e["value"] / h if h else np.nan,
                             "n": point_or_mean(g, f"{var}__suffix__n_tok")})
    return pd.DataFrame(rows)


def consistency_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src, g in runs.groupby("source", sort=False):
        for var in ("a", "b"):
            if f"{var}__share_empty_stem" not in g.columns or g[f"{var}__share_empty_stem"].isna().all():
                continue
            for unit in ("page", "par"):
                pre = f"{var}__cons_{unit}__"
                e = stat(g, pre + "excess")
                rows.append({"source": src, "class": source_class(src), "variant": var, "unit": unit,
                             "within": point_or_mean(g, pre + "within"), "expected": point_or_mean(g, pre + "expected"),
                             "excess": e["value"], "lo": e["lo"], "hi": e["hi"], "z_point": point_or_mean(g, pre + "z"),
                             "global_mode_share": point_or_mean(g, pre + "global"),
                             "n_cells": point_or_mean(g, pre + "n_cells"), "n_tokens": point_or_mean(g, pre + "n_tokens")})
    return pd.DataFrame(rows)


def plot(dec: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    xs = ["pos", "place", "lang", "hand", "sec", "page", "pos+lang", "pos+hand", "all-page", "all"]
    srcs = ["ZL_split", "GC_split", "ctl:caesar_prose", "ctl:celsus_prose", "ctl:pliny_prose",
            "ctl:linnaeus_items", "ts2020:ZL_split", "tmpl:ZL_split", "ac86:ZL_split"]
    colors = ["#1f4e79", "#4f81bd", "#7f7f7f", "#a5a5a5", "#c9c9c9", "#9e480e", "#c0504d", "#e6a23c", "#9bbb59"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for ax, target in zip(axes, R.TARGETS):
        d = dec[(dec["variant"] == "a") & (dec["target"] == target)]
        w = 0.8 / len(srcs)
        for i, (src, c) in enumerate(zip(srcs, colors)):
            s = d[d["source"] == src].set_index("X").reindex(xs)
            if s["explained"].isna().all():
                continue
            pos = np.arange(len(xs)) + (i - len(srcs) / 2) * w + w / 2
            err = np.vstack([s["explained"] - s["lo"], s["hi"] - s["explained"]]).clip(min=0)
            ax.bar(pos, s["explained"], w, color=c, label=src, yerr=np.nan_to_num(err), error_kw={"lw": 0.6})
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(np.arange(len(xs)))
        ax.set_xticklabels(xs, rotation=45, ha="right")
        ax.set_title(f"explained share of H({target} | stem), permutation-corrected")
    axes[0].set_ylabel("(H(S|Z) - H(S|Z,X)) / H(S|Z)")
    axes[1].legend(fontsize=7, ncol=1, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "decomposition.png", dpi=110)
    plt.close(fig)


def summarize() -> None:
    runs = load_all()
    runs.to_csv(OUT / "runs.csv", index=False)
    dec = decomposition_table(runs)
    dec.round(5).to_csv(OUT / "decomposition.csv", index=False)
    residual_table(runs).round(5).to_csv(OUT / "residual_structure.csv", index=False)
    consistency_table(runs).round(5).to_csv(OUT / "consistency.csv", index=False)
    plot(dec)
    write_meta(OUT, "011_suffix_variation", {"sources": sorted(runs["source"].unique().tolist()),
                     "n_rows": int(len(runs)), "by": R.BY, "n_perm": R.N_PERM, "min_stem": R.MIN_STEM,
                     "interval": "point + 2 x quantiles of 80 % page-subsample deviations; nulls: seed quantiles"})
    print("written:", ", ".join(p.name for p in sorted(OUT.glob("*.*"))))


if __name__ == "__main__":
    summarize()
