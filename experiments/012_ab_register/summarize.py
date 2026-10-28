"""Experiment 012: tables and plot from ``results/012_ab_register/runs/*`` and ``extras/*``.

Writes ``runs.csv``, ``stem_sharing.csv``, ``mapping_fit.csv``, ``suffix_mapping_ZL.csv`` /
``suffix_mapping_GC.csv``, ``transfers_top.csv``, ``page_scores.csv``, ``gradient.csv``, ``gradient.png``,
``contributors.csv``, ``sections_within.csv`` and ``meta.json``.  Intervals: 2.5 / 97.5 % quantiles over
the equal-size draws (real) or over the generator seeds (nulls).  ``d_halves`` = per draw
(A vs B) - mean(A|A, B|B), i.e. how far the A/B difference exceeds the difference between two halves of
the same language drawn at the same size.
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
Q = (0.025, 0.975)
CMP_DRAW = ("AB", "AA", "BB", "PERM")
CMP_POINT = ("ABfull", "AAcontig", "BBcontig")


def load_all() -> pd.DataFrame:
    frames = [load_runs(d) for d in sorted((OUT / "runs").iterdir()) if (d / "runs.csv").exists()]
    runs = pd.concat(frames, ignore_index=True)
    runs["kind"] = runs["kind"].fillna("null")
    return runs


def _q(s: pd.Series) -> dict:
    s = s.dropna()
    if not len(s):
        return {"value": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
    return {"value": float(s.mean()), "lo": float(s.quantile(Q[0])), "hi": float(s.quantile(Q[1])), "n": len(s)}


def long_table(runs: pd.DataFrame, block: str, fields) -> pd.DataFrame:
    """``block`` = 'share' (sharing) or a target ('suffix' / 'prefix')."""
    rows = []
    for (tag, source), g in runs.groupby(["set", "source"], sort=False):
        draws = g[g["kind"].isin(["draw", "null"])]
        point = g[g["kind"] == "point"]
        for scope in R.SCOPES:
            for var in ("a", "b"):
                for f in fields:
                    def col(cmp):
                        return f"{scope}__{cmp}__{var}__{block}__{f}"
                    if col("AB") not in g.columns:
                        continue
                    for cmp in CMP_DRAW:
                        if col(cmp) in draws.columns and draws[col(cmp)].notna().any():
                            rows.append({"set": tag, "source": source, "scope": scope, "variant": var, "metric": f,
                                         "comparison": cmp, **_q(draws[col(cmp)])})
                    halves = draws[[col("AA"), col("BB")]].mean(axis=1)
                    rows.append({"set": tag, "source": source, "scope": scope, "variant": var, "metric": f,
                                 "comparison": "AB-halves", **_q(draws[col("AB")] - halves)})
                    for cmp in CMP_POINT:
                        if len(point) and col(cmp) in point.columns:
                            v = float(point[col(cmp)].iloc[0])
                            rows.append({"set": tag, "source": source, "scope": scope, "variant": var, "metric": f,
                                         "comparison": cmp, "value": v, "lo": np.nan, "hi": np.nan, "n": 1})
    return pd.DataFrame(rows)


def extras(name: str) -> pd.DataFrame:
    frames = []
    for d in sorted((OUT / "extras").iterdir()):
        p = d / name
        if p.exists():
            f = pd.read_csv(p, keep_default_na=False, na_values=[""])
            if "set" not in f.columns:
                f.insert(0, "set", d.name)
            frames.append(f)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def null_gradient(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (tag, source), g in runs[runs["kind"] == "null"].groupby(["set", "source"], sort=False):
        for c in [c for c in g.columns if c.startswith("grad__")]:
            rows.append({"set": tag, "source": source, "stat": c[6:], **_q(g[c])})
    return pd.DataFrame(rows)


def plot(scores: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    sets = [s for s in ("ZL_split", "GC_split") if s in set(scores["set"])]
    fig, axes = plt.subplots(len(sets), 1, figsize=(13, 3.6 * len(sets)), squeeze=False)
    markers = {"1": "o", "2": "s", "3": "^", "4": "D", "5": "v"}
    secs = sorted(scores["sec"].astype(str).unique())
    cmap = plt.get_cmap("tab10")
    for ax, tag in zip(axes[:, 0], sets):
        s = scores[scores["set"] == tag].sort_values("order")
        for (h, sec), g in s.groupby([s["hand"].astype(str), s["sec"].astype(str)]):
            ax.scatter(g["order"], g["score"], s=18, marker=markers.get(h, "x"), color=cmap(secs.index(sec) % 10),
                       edgecolor=np.where(g["label"] == "A", "k", "none"), linewidths=0.6,
                       label=f"hand {h}, sec {sec}")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(f"{tag}: page score on the A-B suffix axis (leave-one-page-out; >0 = B-like); black edge = A")
        ax.set_xlabel("page order in the file")
        ax.set_ylabel("score (bits)")
    h, lab = axes[0, 0].get_legend_handles_labels()
    axes[0, 0].legend(h, lab, fontsize=6, ncol=4, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "gradient.png", dpi=110)
    plt.close(fig)


def summarize() -> None:
    runs = load_all()
    runs.to_csv(OUT / "runs.csv", index=False)
    long_table(runs, "share", R.SHARE_F).round(5).to_csv(OUT / "stem_sharing.csv", index=False)
    fits = pd.concat([long_table(runs, t, R.MAP_F).assign(target=t) for t in R.TARGETS], ignore_index=True)
    fits.round(5).to_csv(OUT / "mapping_fit.csv", index=False)
    mats = extras("mapping.csv")
    if len(mats):
        for key in ("ZL", "GC"):
            m = mats[mats["set"] == f"{key}_split"]
            m.round(5).to_csv(OUT / f"suffix_mapping_{key}.csv", index=False)
        off = mats[(mats["affix_A"] != mats["affix_B"])]
        top = (off.sort_values("excess", ascending=False).groupby(["set", "variant", "target"], sort=False).head(10))
        top.round(5).to_csv(OUT / "transfers_top.csv", index=False)
    sc = extras("page_scores.csv")
    if len(sc):
        sc.round(5).to_csv(OUT / "page_scores.csv", index=False)
        plot(sc)
    gr = extras("gradient.csv")
    ng = null_gradient(runs)
    with open(OUT / "gradient.csv", "w", encoding="utf-8") as fh:
        gr.round(5).to_csv(fh, index=False)
    ng.round(5).to_csv(OUT / "gradient_nulls.csv", index=False)
    extras("contributors.csv").round(6).to_csv(OUT / "contributors.csv", index=False)
    extras("sections_within.csv").round(6).to_csv(OUT / "sections_within.csv", index=False)
    write_meta(OUT, "012_ab_register", {"sets": sorted(runs["set"].unique().tolist()),
                                        "sources": sorted(runs["source"].unique().tolist()), "n_rows": int(len(runs)),
                                        "n_draw": R.N_DRAW, "n_null": R.N_NULL, "min_share": R.MIN_SHARE,
                                        "min_map": R.MIN_MAP, "sec_n": R.SEC_N, "sec_reps": R.SEC_REPS})
    print("written:", ", ".join(p.name for p in sorted(OUT.glob("*.*"))))


if __name__ == "__main__":
    summarize()
