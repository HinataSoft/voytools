"""Aggregation of experiment 015 runs into the spec's tables and figures (``run.py summarize``).

Reads ``results/015_rare_stem_recurrence/runs/<combo>/<source>/runs.csv`` and writes, in the results
folder: ``runs.csv`` (all runs), ``summary.csv`` (long table: combo x source x variant x band x metric,
value with 2.5 / 97.5 % interval), ``distance_profile.csv``, ``far_pairs_sections.csv``,
``form_by_distance.csv``, ``far_clusters.csv``, ``label_stems.csv``, ``meta.json`` and the figures
``distance_profile.png`` / ``form_by_distance.png``.

Value of a source = the real run (``real``, interval = stem bootstrap ``real_boot``) or the mean over
seeds (nulls, generators; interval = quantiles over seeds).  Expectations (``e_*``, ``*_exp``) are the
exact size-weighted page-permutation expectations computed per sample.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from voytools.experiment import write_meta  # noqa: E402
from voytools.recurrence import DISTANCE_LABELS, FORM_LABELS, SAME_COLUMNS  # noqa: E402

OUT = ROOT / "results" / "015_rare_stem_recurrence"
RUNS = OUT / "runs"
Q = (0.025, 0.975)
COL = re.compile(r"^([ab])__([rcl])__(.+)$")
SOURCE_ORDER = ["real", "perm_tokens", "perm_global", "perm_section", "selfcitation_21", "template_005",
                "template_005_notopic", "autocopy_86"]
PLOT_SOURCES = ["real", "perm_global", "selfcitation_21", "template_005", "template_005_notopic", "autocopy_86",
                "ctl:pliny_prose", "ctl:celsus_prose", "ctl:linnaeus_items", "ctl:gibberish_pooled"]


def load_all() -> pd.DataFrame:
    frames = []
    for p in sorted(RUNS.glob("*/*/runs.csv")):
        df = pd.read_csv(p)
        frames.append(df)
    runs = pd.concat(frames, ignore_index=True)
    return runs


def long_summary(runs: pd.DataFrame) -> pd.DataFrame:
    """One row per combo x source x variant x band x metric with value, lo, hi, n_seeds."""
    metric_cols = [c for c in runs.columns if COL.match(c)]
    rows = []
    for (combo, source), g in runs.groupby(["combo", "source"], sort=False):
        if source == "real_boot":
            continue
        boot = runs[(runs["combo"] == combo) & (runs["source"] == "real_boot")] if source == "real" else None
        vals = g[metric_cols].apply(pd.to_numeric, errors="coerce")
        mean = vals.iloc[0] if source == "real" else vals.mean()
        spread = boot[metric_cols].apply(pd.to_numeric, errors="coerce") if boot is not None and len(boot) else vals
        lo, hi = spread.quantile(Q[0]), spread.quantile(Q[1])
        for c in metric_cols:
            var, band, metric = COL.match(c).groups()
            rows.append({"combo": combo, "source": source, "variant": var, "band": band, "metric": metric,
                         "value": mean[c], "lo": lo[c], "hi": hi[c],
                         "n_seeds": int(len(spread)) if source != "real" else int(len(spread))})
    return pd.DataFrame(rows)


def _pick(s: pd.DataFrame, metric: str) -> pd.DataFrame:
    return s[s["metric"] == metric].set_index(["combo", "source", "variant", "band"])


def distance_profile(s: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = s[["combo", "source", "variant", "band"]].drop_duplicates()
    for _, k in keys.iterrows():
        idx = tuple(k)
        n = _pick(s, "n_pairs").loc[idx, "value"]
        ns = _pick(s, "n_stems").loc[idx, "value"]
        for b in list(DISTANCE_LABELS) + ["far", "far_nl"]:
            p = _pick(s, f"p_{b}").loc[idx]
            e = _pick(s, f"e_{b}").loc[idx]
            x = _pick(s, f"x_{b}").loc[idx] if b.startswith("far") else None
            rows.append({**k.to_dict(), "bin": b, "n_stems": ns, "n_pairs": n, "share": p["value"],
                         "share_lo": p["lo"], "share_hi": p["hi"], "expected": e["value"],
                         "excess": p["value"] - e["value"] if x is None else x["value"],
                         "excess_lo": p["lo"] - e["value"] if x is None else x["lo"],
                         "excess_hi": p["hi"] - e["value"] if x is None else x["hi"],
                         "ratio": p["value"] / e["value"] if e["value"] else np.nan})
    return pd.DataFrame(rows)


def far_pairs_sections(s: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = s[["combo", "source", "variant", "band"]].drop_duplicates()
    for _, k in keys.iterrows():
        idx = tuple(k)
        for rng in ("far", "mid"):
            n = _pick(s, f"{rng}_n").loc[idx, "value"]
            for c in SAME_COLUMNS:
                p = _pick(s, f"{rng}_{c}").loc[idx]
                e = _pick(s, f"{rng}_{c}_exp").loc[idx]
                rows.append({**k.to_dict(), "range": rng, "attribute": c, "n_pairs": n, "share": p["value"],
                             "share_lo": p["lo"], "share_hi": p["hi"], "expected": e["value"],
                             "excess": p["value"] - e["value"], "excess_lo": p["lo"] - e["value"],
                             "excess_hi": p["hi"] - e["value"]})
    return pd.DataFrame(rows)


def form_by_distance(s: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = s[["combo", "source", "variant", "band"]].drop_duplicates()
    for _, k in keys.iterrows():
        idx = tuple(k)
        for c in FORM_LABELS:
            n = _pick(s, f"form_n_{c}").loc[idx, "value"]
            p = _pick(s, f"form_same_{c}").loc[idx]
            e = _pick(s, f"form_exp_{c}").loc[idx]
            rows.append({**k.to_dict(), "class": c, "n_pairs": n, "same_form": p["value"], "same_lo": p["lo"],
                         "same_hi": p["hi"], "expected": e["value"], "diff": p["value"] - e["value"],
                         "diff_lo": p["lo"] - e["value"], "diff_hi": p["hi"] - e["value"]})
    return pd.DataFrame(rows)


def far_clusters(s: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = s[["combo", "source", "variant", "band"]].drop_duplicates()
    for _, k in keys.iterrows():
        idx = tuple(k)
        for tag, by in (("clq", "quire"), ("cls", "section")):
            n = _pick(s, f"{tag}_n").loc[idx]
            rows.append({**k.to_dict(), "by": by, "n_stems_ge3": _pick(s, f"{tag}_n_stems").loc[idx, "value"],
                         "n_clusters": n["value"], "n_lo": n["lo"], "n_hi": n["hi"],
                         "expected": _pick(s, f"{tag}_exp").loc[idx, "value"],
                         "share": _pick(s, f"{tag}_share").loc[idx, "value"],
                         "share_expected": _pick(s, f"{tag}_share_exp").loc[idx, "value"]})
    out = pd.DataFrame(rows)
    out["ratio"] = out["n_clusters"] / out["expected"]
    return out


def label_stems() -> pd.DataFrame:
    frames = []
    for p in sorted(RUNS.glob("*/stems.csv")):
        df = pd.read_csv(p)
        df = df[df["band"] == "l"].assign(combo=p.parent.name)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out[["combo", "variant", "stem", "count", "n_pages", "n_quires", "n_forms"]]


def plots(prof: pd.DataFrame, form: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, combo in zip(axes, ("ZL_split", "GC_split")):
        d = prof[(prof["variant"] == "a") & (prof["band"] == "r") & prof["bin"].isin(DISTANCE_LABELS)]
        for src in PLOT_SOURCES:
            c = combo if not src.startswith("ctl:") else src
            dd = d[(d["combo"] == c) & (d["source"] == ("real" if src.startswith("ctl:") else src))]
            if dd.empty:
                continue
            dd = dd.set_index("bin").loc[list(DISTANCE_LABELS)]
            ax.plot(range(len(DISTANCE_LABELS)), dd["ratio"], marker="o", label=src,
                    lw=2.2 if src == "real" else 1.2)
        ax.axhline(1, color="k", lw=0.6)
        ax.set_xticks(range(len(DISTANCE_LABELS)))
        ax.set_xticklabels(DISTANCE_LABELS)
        ax.set_yscale("log")
        ax.set_title(f"{combo}: stems 2-5, share / expected per distance bin")
        ax.set_xlabel("page distance")
    axes[0].set_ylabel("observed / expected (page permutation)")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "distance_profile.png", dpi=120)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, combo in zip(axes, ("ZL_split", "GC_split")):
        d = form[(form["variant"] == "a") & (form["band"] == "r")]
        for src in PLOT_SOURCES:
            c = combo if not src.startswith("ctl:") else src
            dd = d[(d["combo"] == c) & (d["source"] == ("real" if src.startswith("ctl:") else src))]
            if dd.empty:
                continue
            dd = dd.set_index("class").loc[list(FORM_LABELS)]
            ax.plot(range(len(FORM_LABELS)), dd["diff"], marker="o", label=src, lw=2.2 if src == "real" else 1.2)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(len(FORM_LABELS)))
        ax.set_xticklabels(FORM_LABELS)
        ax.set_title(f"{combo}: identical form share - paradigm expectation")
        ax.set_xlabel("page distance class")
    axes[0].set_ylabel("same form - expected")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "form_by_distance.png", dpi=120)
    plt.close(fig)


def summarize() -> None:
    runs = load_all()
    runs.to_csv(OUT / "runs.csv", index=False)
    s = long_summary(runs)
    s.to_csv(OUT / "summary.csv", index=False)
    prof = distance_profile(s, runs)
    prof.to_csv(OUT / "distance_profile.csv", index=False)
    far_pairs_sections(s).to_csv(OUT / "far_pairs_sections.csv", index=False)
    form = form_by_distance(s)
    form.to_csv(OUT / "form_by_distance.csv", index=False)
    far_clusters(s).to_csv(OUT / "far_clusters.csv", index=False)
    label_stems().to_csv(OUT / "label_stems.csv", index=False)
    metas = {}
    for p in sorted(RUNS.glob("*/*/meta.json")):
        metas[f"{p.parent.parent.name}/{p.parent.name}"] = json.loads(p.read_text())
    counts = runs.groupby(["combo", "source"]).size().to_dict()
    write_meta(OUT, "015_rare_stem_recurrence", params={"runs": {f"{k[0]}/{k[1]}": int(v) for k, v in counts.items()},
                                                        "parts": metas})
    plots(prof, form)
    print(f"{len(runs)} runs, {len(s)} summary rows")


if __name__ == "__main__":
    summarize()
