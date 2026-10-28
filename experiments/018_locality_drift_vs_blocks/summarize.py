"""Aggregation of experiment 018 runs into the spec's tables and figures (``run.py summarize``).

Reads ``results/018_locality_drift_vs_blocks/runs/<combo>/<source>/runs.csv`` (and the ``extras`` files
of the real text) and writes, in the results folder: ``runs.csv`` (all runs without the 200-seed page
permutations, key columns only; the full rows stay in ``runs/``), ``summary.csv`` (long table of the key
metrics with intervals), ``similarity_profile.csv``, ``segmentation_gain.csv``, ``crossval_blocks.csv``,
``changepoints.csv`` (precision / recall summary) + ``changepoints_list.csv`` (the boundaries of the real
text), ``physical_order.csv``, ``neighbours.csv`` + ``neighbours_list.csv``, ``quire_matrix.csv``,
``meta.json`` and the figures ``distance_profile.png``, ``segmentation_gain.png``,
``similarity_matrix.png``, ``quire_matrix.png``.

Value of a source = mean over its rows (real: over the 50 stem splits; nulls and generators: over
seeds), interval = 2.5 / 97.5 % quantiles over the rows.
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

OUT = ROOT / "results" / "018_locality_drift_vs_blocks"
RUNS = OUT / "runs"
Q = (0.025, 0.975)
COL = re.compile(r"^(rare|cos|jac)__(all|within|herbalA)__(.+)$")
MEASURES = ("rare", "cos", "jac")
CUTS = ("all", "within", "herbalA")
REFS = ("section", "language", "hand", "quire")
KEY_K = (2, 3, 5, 8, 10, 15, 20, 25, 30)
BIG_SOURCES = ("perm_global", "perm_section")
SOURCE_ORDER = ["real", "perm_global", "perm_section", "selfcitation_21", "selfcitation_mem", "template_005"]
# Okabe-Ito (colour-blind safe), fixed per source / control
COLORS = {"real": "#000000", "perm_global": "#999999", "perm_section": "#E69F00", "selfcitation_21": "#56B4E9",
          "selfcitation_mem": "#0072B2", "template_005": "#009E73", "ctl:pliny_prose": "#D55E00",
          "ctl:celsus_prose": "#CC79A7", "ctl:linnaeus_items": "#F0E442", "ctl:gibberish_pooled": "#7f4f24"}
PLOT_SOURCES = ["real", "perm_global", "perm_section", "selfcitation_21", "selfcitation_mem", "template_005",
                "ctl:pliny_prose", "ctl:celsus_prose", "ctl:linnaeus_items", "ctl:gibberish_pooled"]


def load_all() -> pd.DataFrame:
    frames = [pd.read_csv(p) for p in sorted(RUNS.glob("*/*/runs.csv")) if p.parent.parent.name != "memory"]
    return pd.concat(frames, ignore_index=True)


def _agg(g: pd.DataFrame, cols: list[str]) -> dict:
    out = {}
    for c in cols:
        if c not in g.columns:
            continue
        s = pd.to_numeric(g[c], errors="coerce").dropna()
        if s.empty:
            continue
        out[c] = (float(s.mean()), float(s.quantile(Q[0])), float(s.quantile(Q[1])), int(len(s)))
    return out


def profile_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (combo, source), g in runs.groupby(["combo", "source"], sort=False):
        for m in MEASURES:
            for cut in CUTS:
                cols = [c for c in g.columns if c.startswith(f"{m}__{cut}__prof_d")]
                for c, (v, lo, hi, n) in _agg(g, cols).items():
                    d = int(c.rsplit("prof_d", 1)[1])
                    rows.append({"combo": combo, "source": source, "measure": m, "cut": cut, "d": d,
                                 "ratio": v, "lo": lo, "hi": hi, "n": n})
    return pd.DataFrame(rows)


def gain_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (combo, source), g in runs.groupby(["combo", "source"], sort=False):
        for m in MEASURES:
            for cut in CUTS:
                kref = g.get(f"{m}__{cut}__K_ref")
                if kref is None or kref.dropna().empty:
                    continue
                kref = int(kref.dropna().iloc[0])
                cv = _agg(g, [f"{m}__{cut}__cv_K{k}" for k in range(1, 31)])
                fit = _agg(g, [f"{m}__{cut}__fit_K{k}" for k in range(1, 31)])
                for k in range(1, 31):
                    c = f"{m}__{cut}__cv_K{k}"
                    if c not in cv:
                        continue
                    v, lo, hi, n = cv[c]
                    rows.append({"combo": combo, "source": source, "measure": m, "cut": cut, "K": k, "K_ref": kref,
                                 "gain_cv": v, "lo": lo, "hi": hi, "gain_fit": fit.get(f"{m}__{cut}__fit_K{k}", (np.nan,))[0],
                                 "n": n})
    return pd.DataFrame(rows)


def crossval_table(runs: pd.DataFrame) -> pd.DataFrame:
    real = runs[runs["source"] == "real"]
    rows = []
    for _, r in real.iterrows():
        for m in MEASURES:
            for cut in CUTS:
                kref = r.get(f"{m}__{cut}__K_ref")
                if kref is None or not np.isfinite(kref):
                    continue
                ks = sorted(set(KEY_K) | {min(int(kref) * x, 30) for x in (1, 2, 3)})
                for k in ks:
                    rows.append({"combo": r["combo"], "measure": m, "cut": cut, "split_seed": int(r["seed"]), "K": k,
                                 "K_ref": int(kref), "gain_cv": r.get(f"{m}__{cut}__cv_K{k}", np.nan)})
    return pd.DataFrame(rows)


def changepoint_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (combo, source), g in runs.groupby(["combo", "source"], sort=False):
        for m in MEASURES:
            for cut in CUTS:
                if f"{m}__{cut}__K_ref" not in g.columns or g[f"{m}__{cut}__K_ref"].dropna().empty:
                    continue
                for mult in (1, 2, 3):
                    kcol = f"{m}__{cut}__cp{mult}x_K"
                    for ref in REFS:
                        base = f"{m}__{cut}__cp{mult}x_{ref}_"
                        a = _agg(g, [base + x for x in ("prec", "rec", "prec0", "rec0")])
                        if any(base + x not in a for x in ("prec", "rec", "prec0", "rec0")):
                            continue
                        rows.append({"combo": combo, "source": source, "measure": m, "cut": cut, "K_mult": mult,
                                     "K": int(g[kcol].dropna().iloc[0]), "ref": ref,
                                     "n_ref": int(g[base + "n_ref"].dropna().iloc[0]),
                                     "prec": a[base + "prec"][0], "prec_lo": a[base + "prec"][1], "prec_hi": a[base + "prec"][2],
                                     "rec": a[base + "rec"][0], "rec_lo": a[base + "rec"][1], "rec_hi": a[base + "rec"][2],
                                     "prec_random": a[base + "prec0"][0], "rec_random": a[base + "rec0"][0],
                                     "n": a[base + "prec"][3]})
    return pd.DataFrame(rows)


def physical_table(runs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    phys, nb = [], []
    for (combo, source), g in runs.groupby(["combo", "source"], sort=False):
        for m in MEASURES:
            pre = f"{m}__all__"
            for unit in ("bifolio", "quire"):
                a = _agg(g, [f"{pre}phys_{unit}_{x}" for x in ("same", "matched", "excess")])
                if f"{pre}phys_{unit}_excess" not in a:
                    continue
                n_same = g[f"{pre}phys_{unit}_n_same"].dropna()
                phys.append({"combo": combo, "source": source, "measure": m, "unit": unit,
                             "n_same_pairs": int(n_same.iloc[0]) if len(n_same) else np.nan,
                             "same": a[f"{pre}phys_{unit}_same"][0], "matched": a[f"{pre}phys_{unit}_matched"][0],
                             "excess": a[f"{pre}phys_{unit}_excess"][0], "lo": a[f"{pre}phys_{unit}_excess"][1],
                             "hi": a[f"{pre}phys_{unit}_excess"][2], "n": a[f"{pre}phys_{unit}_excess"][3]})
            a = _agg(g, [f"{pre}qo_adjacent_share", f"{pre}qo_rho"])
            for k in ("adjacent_share", "rho"):
                if f"{pre}qo_{k}" in a:
                    v, lo, hi, n = a[f"{pre}qo_{k}"]
                    phys.append({"combo": combo, "source": source, "measure": m, "unit": f"quire_order_{k}",
                                 "n_same_pairs": np.nan, "same": np.nan, "matched": np.nan, "excess": v, "lo": lo,
                                 "hi": hi, "n": n})
            a = _agg(g, [f"{pre}nb_{x}" for x in ("share_le2", "share_le5", "share_gt16", "n_anomalous", "n_pages")])
            if f"{pre}nb_share_le5" in a:
                nb.append({"combo": combo, "source": source, "measure": m,
                           **{k: a[f"{pre}nb_{k}"][0] for k in ("n_pages", "share_le2", "share_le5", "share_gt16", "n_anomalous")},
                           **{f"{k}_lo": a[f"{pre}nb_{k}"][1] for k in ("share_le2", "share_le5", "share_gt16", "n_anomalous")},
                           **{f"{k}_hi": a[f"{pre}nb_{k}"][2] for k in ("share_le2", "share_le5", "share_gt16", "n_anomalous")},
                           "n": a[f"{pre}nb_share_le5"][3]})
    return pd.DataFrame(phys), pd.DataFrame(nb)


def summary_table(gain: pd.DataFrame, cp: pd.DataFrame, phys: pd.DataFrame, nb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in gain.iterrows():
        if r["K"] in {min(r["K_ref"] * x, 30) for x in (1, 2, 3)} or r["K"] in (2, 10, 30):
            rows.append({"combo": r["combo"], "source": r["source"], "measure": r["measure"], "cut": r["cut"],
                         "metric": f"gain_cv_K{int(r['K'])}" + ("_ref" if r["K"] == r["K_ref"] else ""),
                         "value": r["gain_cv"], "lo": r["lo"], "hi": r["hi"], "n": r["n"]})
    for _, r in cp.iterrows():
        if r["K_mult"] == 1:
            rows.append({"combo": r["combo"], "source": r["source"], "measure": r["measure"], "cut": r["cut"],
                         "metric": f"cp_prec_{r['ref']}", "value": r["prec"], "lo": r["prec_lo"], "hi": r["prec_hi"], "n": r["n"]})
            rows.append({"combo": r["combo"], "source": r["source"], "measure": r["measure"], "cut": r["cut"],
                         "metric": f"cp_prec_random_{r['ref']}", "value": r["prec_random"], "lo": np.nan, "hi": np.nan, "n": r["n"]})
    for _, r in phys.iterrows():
        rows.append({"combo": r["combo"], "source": r["source"], "measure": r["measure"], "cut": "all",
                     "metric": f"phys_{r['unit']}", "value": r["excess"], "lo": r["lo"], "hi": r["hi"], "n": r["n"]})
    for _, r in nb.iterrows():
        for k in ("share_le2", "share_le5", "share_gt16", "n_anomalous"):
            rows.append({"combo": r["combo"], "source": r["source"], "measure": r["measure"], "cut": "all",
                         "metric": f"nb_{k}", "value": r[k], "lo": r[f"{k}_lo"], "hi": r[f"{k}_hi"], "n": r["n"]})
    return pd.DataFrame(rows)


def extras_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cps, nbs, qms = [], [], []
    for d in sorted(RUNS.glob("*")):
        if (d / "changepoints.csv").exists():
            cps.append(pd.read_csv(d / "changepoints.csv"))
        if (d / "neighbours.csv").exists():
            nbs.append(pd.read_csv(d / "neighbours.csv"))
        if (d / "quire_matrix.csv").exists():
            qms.append(pd.read_csv(d / "quire_matrix.csv"))
    cat = lambda xs: pd.concat(xs, ignore_index=True) if xs else pd.DataFrame()  # noqa: E731
    return cat(cps), cat(nbs), cat(qms)


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #


def _shared_legend(fig, axes) -> None:
    handles = {}
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(l, h)
    fig.legend(handles.values(), handles.keys(), fontsize=7, loc="lower center", ncol=5, frameon=False)


def _series(table: pd.DataFrame, combo: str, src: str) -> pd.DataFrame:
    if src.startswith("ctl:"):
        return table[(table["combo"] == src) & (table["source"] == "real")]
    return table[(table["combo"] == combo) & (table["source"] == src)]


def plots(prof: pd.DataFrame, gain: pd.DataFrame, qm: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # distance profile: rare stem sharing, pooled ratio per distance (cut all)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, combo in zip(axes, ("ZL_split", "GC_split")):
        for src in PLOT_SOURCES:
            dd = _series(prof[(prof["measure"] == "rare") & (prof["cut"] == "all")], combo, src).sort_values("d")
            dd = dd[dd["d"] <= 40]
            if dd.empty:
                continue
            ax.plot(dd["d"], dd["ratio"], label=src, color=COLORS[src], lw=2.2 if src == "real" else 1.1,
                    marker="o" if src == "real" else None, ms=3)
        ax.axhline(1, color="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_title(f"{combo}: rare stems 2-5, observed / expected pairs per page distance")
        ax.set_xlabel("page distance d")
    axes[0].set_ylabel("observed / expected (size-weighted page permutation)")
    _shared_legend(fig, axes)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(OUT / "distance_profile.png", dpi=120)
    plt.close(fig)

    # segmentation gain G(K)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5))
    for row, measure in enumerate(("rare", "cos")):
        for ax, (combo, cut) in zip(axes[row], (("ZL_split", "all"), ("ZL_split", "herbalA"), ("GC_split", "all"))):
            for src in PLOT_SOURCES:
                dd = _series(gain[(gain["measure"] == measure) & (gain["cut"] == cut)], combo, src).sort_values("K")
                if dd.empty:
                    continue
                ax.plot(dd["K"], dd["gain_cv"], label=src, color=COLORS[src], lw=2.2 if src == "real" else 1.1)
                if src in ("real", "perm_global", "selfcitation_mem"):
                    ax.fill_between(dd["K"], dd["lo"], dd["hi"], color=COLORS[src], alpha=0.12, lw=0)
            ax.axhline(0, color="k", lw=0.6)
            ax.set_title(f"{combo} / {cut}: {measure}, cross-validated block gain", fontsize=9)
            ax.set_xlabel("K blocks")
        axes[row][0].set_ylabel("held-out pooled residual inside blocks")
    _shared_legend(fig, axes.ravel())
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT / "segmentation_gain.png", dpi=120)
    plt.close(fig)

    # residual matrices of the real ZL split text with section and quire boundaries
    pages = pd.read_csv(RUNS / "ZL_split" / "pages.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))
    for ax, measure, vmax in zip(axes, ("rare_smooth", "cos"), (1.0, 0.15)):
        p = RUNS / "ZL_split" / f"residual_{measure}.csv"
        if not p.exists():
            continue
        m = pd.read_csv(p, index_col=0).to_numpy(dtype=float)
        im = ax.imshow(np.clip(np.nan_to_num(m, nan=0.0), -vmax, vmax), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       interpolation="nearest")
        for col, color, ls in (("section", "k", "-"), ("quire", "#444444", ":")):
            b = np.flatnonzero(pages[col].astype(str).to_numpy()[1:] != pages[col].astype(str).to_numpy()[:-1]) + 1
            for x in b:
                ax.axhline(x - 0.5, color=color, lw=0.5, ls=ls)
                ax.axvline(x - 0.5, color=color, lw=0.5, ls=ls)
        ax.set_title(f"ZL split, {measure}: residual vs distance model (solid: section, dotted: quire)"
                     + (" - 7-page window" if measure == "rare_smooth" else ""), fontsize=9)
        ax.set_xlabel("page (file order)")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    fig.tight_layout()
    fig.savefig(OUT / "similarity_matrix.png", dpi=120)
    plt.close(fig)

    # quire x quire matrices
    if not qm.empty:
        sel = qm[(qm["combo"] == "ZL_split")]
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
        for ax, measure, vmax in zip(axes, ("rare", "cos"), (1.0, 0.05)):
            s = sel[sel["measure"] == measure]
            if s.empty:
                continue
            names = s["quire"].tolist()
            mat = s[names].to_numpy(dtype=float)
            im = ax.imshow(np.clip(np.nan_to_num(mat, nan=0.0), -vmax, vmax), cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            ax.set_xticks(range(len(names)))
            ax.set_xticklabels(names, fontsize=7)
            ax.set_yticks(range(len(names)))
            ax.set_yticklabels(names, fontsize=7)
            ax.set_title(f"ZL split, {measure}: quire x quire pooled residual", fontsize=9)
            fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        fig.tight_layout()
        fig.savefig(OUT / "quire_matrix.png", dpi=120)
        plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def summarize() -> None:
    runs = load_all()
    prof = profile_table(runs)
    gain = gain_table(runs)
    cv = crossval_table(runs)
    cp = changepoint_table(runs)
    phys, nb = physical_table(runs)
    summ = summary_table(gain, cp, phys, nb)
    cps, nbs, qms = extras_tables()
    prof.to_csv(OUT / "similarity_profile.csv", index=False)
    gain.to_csv(OUT / "segmentation_gain.csv", index=False)
    cv.to_csv(OUT / "crossval_blocks.csv", index=False)
    cp.to_csv(OUT / "changepoints.csv", index=False)
    cps.to_csv(OUT / "changepoints_list.csv", index=False)
    phys.to_csv(OUT / "physical_order.csv", index=False)
    nb.to_csv(OUT / "neighbours.csv", index=False)
    nbs.to_csv(OUT / "neighbours_list.csv", index=False)
    qms.to_csv(OUT / "quire_matrix.csv", index=False)
    summ.to_csv(OUT / "summary.csv", index=False)
    keep = [c for c in runs.columns if not COL.match(c)]
    keep += [c for c in runs.columns if COL.match(c) and (
        re.search(r"__cv_K(2|5|10|20|30)$", c) or re.search(r"__cp\dx_\w+_(prec|rec)$", c) or "__phys_" in c
        or "__nb_" in c or "__qo_" in c or c.endswith("__K_ref"))]
    runs[~runs["source"].isin(BIG_SOURCES)][keep].to_csv(OUT / "runs.csv", index=False)
    metas = {}
    for p in sorted(RUNS.glob("*/*/meta.json")):
        metas[f"{p.parent.parent.name}/{p.parent.name}"] = json.loads(p.read_text(encoding="utf-8"))
    write_meta(OUT, "018_locality_drift_vs_blocks",
               params={"runs": {k: v.get("n_seeds_done", v.get("n_seeds_requested")) for k, v in metas.items()},
                       "sources": {k: v.get("params", {}) for k, v in metas.items()}},
               n_rows=int(len(runs)))
    plots(prof, gain, qms)
    print(f"summary written: {len(runs)} runs, {len(summ)} summary rows")


if __name__ == "__main__":
    summarize()
