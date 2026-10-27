"""Experiment 010: tables and plots from ``results/010_line_as_item/runs/*``.

Every metric is summarised per source (real text / controls: full-text point estimate with the 007
interval from 80 % page subsamples, SD x 2; nulls: mean and 2.5-97.5 % over seeds).  For every real
tag, z against each null of the same tag and against each control.  Tables: ``boundary_vs_inside.csv``,
``line_header.csv``, ``line_templates.csv``, ``same_length.csv``, ``vertical_copy.csv``, ``line_end.csv``
(long format: metric, source, value, lo, hi, sd, z vs null/controls), ``signatures_top20.csv``,
``summary.csv`` (all), ``runs.csv``; plots ``boundary.png``, ``vertical_copy.png``, ``templates.png``.
"""

from __future__ import annotations

import json
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
SCALE = float(np.sqrt(R.R9.SUB_FRAC / (1 - R.R9.SUB_FRAC)))
Q = (0.025, 0.975)
TABLES = {"boundary_vs_inside": ("bnd_", "n_pairs"), "line_header": ("head_",), "line_templates": ("tpl_",),
          "same_length": ("jac_",), "vertical_copy": ("vert_",), "line_end": ("end_",)}
SKIP = {"source", "kind", "seed", "elapsed_s", "end_enrich_glyph", "n_tokens", "n_lines"}


def source_class(src: str) -> str:
    if src in R.REAL_TAGS:
        return "real"
    if src.startswith("ctl:"):
        return "items" if src.endswith("_items") else "control"
    return "nullmodel"


def load_all() -> pd.DataFrame:
    runs = pd.concat([load_runs(d) for d in sorted((OUT / "runs").iterdir()) if (d / "runs.csv").exists()],
                     ignore_index=True)
    runs["kind"] = runs["kind"].fillna("null")  # read_csv parses the label 'null' as NA
    return runs


def metrics_of(runs: pd.DataFrame) -> list[str]:
    return [c for c in runs.columns if c not in SKIP and pd.api.types.is_numeric_dtype(runs[c])]


def summarize_sources(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src, g in runs.groupby("source", sort=False):
        point, sub, nul = g[g["kind"] == "point"], g[g["kind"] == "sub"], g[g["kind"] == "null"]
        for m in metrics_of(runs):
            if len(point):
                v = float(point[m].iloc[0])
                s = sub[m].dropna()
                dev = s - s.mean()
                r = {"value": v, "lo": v + SCALE * float(dev.quantile(Q[0])) if len(s) else np.nan,
                     "hi": v + SCALE * float(dev.quantile(Q[1])) if len(s) else np.nan,
                     "sd": SCALE * float(s.std(ddof=1)) if len(s) > 1 else np.nan, "n": len(s)}
            else:
                s = nul[m].dropna()
                r = {"value": float(s.mean()) if len(s) else np.nan,
                     "lo": float(s.quantile(Q[0])) if len(s) else np.nan,
                     "hi": float(s.quantile(Q[1])) if len(s) else np.nan,
                     "sd": float(s.std(ddof=1)) if len(s) > 1 else np.nan, "n": len(s)}
            rows.append({"metric": m, "source": src, "class": source_class(src), **r})
    return pd.DataFrame(rows)


def add_z(summ: pd.DataFrame) -> pd.DataFrame:
    """z of every real tag against each null of the same tag and each control (SD of the other source)."""
    P = summ.set_index(["source", "metric"])
    rows = []
    for tag in R.REAL_TAGS:
        if tag not in summ["source"].values:
            continue
        others = [s for s in summ["source"].unique()
                  if s.startswith("ctl:") or (":" in s and not s.startswith("ctl:") and s.split(":")[1] == tag)]
        for m in summ["metric"].unique():
            v = P.loc[(tag, m)]
            row = {"metric": m, "real": tag, "real_value": v["value"], "real_lo": v["lo"], "real_hi": v["hi"]}
            for o in others:
                s = P.loc[(o, m)]
                key = o.replace("ctl:", "ctl_").split(":")[0]
                row[f"{key}_value"] = s["value"]
                row[f"z_{key}"] = (v["value"] - s["value"]) / s["sd"] if s["sd"] and s["sd"] == s["sd"] else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def signatures() -> pd.DataFrame:
    d = OUT / "signatures"
    files = sorted(d.glob("*.csv")) if d.exists() else []
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        return None


def _order(summ: pd.DataFrame) -> list[str]:
    srcs = summ["source"].unique().tolist()
    return ([s for s in R.REAL_TAGS if s in srcs] + sorted(s for s in srcs if s.startswith("ctl:"))
            + sorted(s for s in srcs if source_class(s) == "nullmodel"))


def _bars(ax, summ: pd.DataFrame, metric: str, order: list[str], title: str) -> None:
    sub = summ[summ["metric"] == metric].set_index("source").loc[order]
    colors = {"real": "#1f4e79", "control": "#7f7f7f", "items": "#c55a11", "nullmodel": "#70ad47"}
    y = np.arange(len(order))
    v = sub["value"].to_numpy(dtype=float)
    err = np.vstack([v - sub["lo"].to_numpy(dtype=float), sub["hi"].to_numpy(dtype=float) - v])
    ax.barh(y, v, xerr=np.clip(np.nan_to_num(err), 0, None), color=[colors[source_class(s)] for s in order],
            error_kw={"lw": 0.8})
    ax.set_yticks(y)
    ax.set_yticklabels([s.replace("ctl:", "") for s in order], fontsize=7)
    ax.invert_yaxis()
    ax.axvline(0, color="k", lw=0.6)
    ax.set_title(title, fontsize=9)


def plots(summ: pd.DataFrame) -> None:
    plt = _plt()
    if plt is None:
        return
    order = _order(summ)
    for name, panels in {
        "boundary": [("bnd_edge_ratio_line", "edge glyphs: excess across line / inside"),
                     ("bnd_affix_ratio_line", "affix classes: across line / inside"),
                     ("bnd_stem_ratio_line", "stems: across line / inside")],
        "vertical_copy": [("vert_excess_above", "match (edit <= 1) with token above minus random position"),
                          ("vert_ratio_above", "match rate above / elsewhere in line above")],
        "templates": [("head_second", "I(first; second class), corrected"),
                      ("tpl_top20_cov", "coverage of 20 most frequent line signatures (4-8 tokens)"),
                      ("tpl_zgain_colsh", "zlib ratio gain of same-length column shuffle")],
    }.items():
        fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 0.22 * len(order) + 1.5))
        for ax, (m, t) in zip(np.atleast_1d(axes), panels):
            _bars(ax, summ, m, order, t)
        fig.tight_layout()
        fig.savefig(OUT / f"{name}.png", dpi=110)
        plt.close(fig)


def summarize() -> None:
    runs = load_all()
    runs.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    summ = summarize_sources(runs)
    summ.to_csv(OUT / "summary.csv", index=False, float_format="%.6g")
    z = add_z(summ)
    for name, prefixes in TABLES.items():
        keep = summ["metric"].str.startswith(prefixes)
        summ[keep].to_csv(OUT / f"{name}.csv", index=False, float_format="%.6g")
        zk = z[z["metric"].str.startswith(prefixes)]
        zk.to_csv(OUT / f"{name}_z.csv", index=False, float_format="%.4g")
    sig = signatures()
    if len(sig):
        sig.to_csv(OUT / "signatures_top20.csv", index=False, float_format="%.5g")
    plots(summ)
    counts = runs.groupby(["source", "kind"]).size().reset_index(name="n")
    write_meta(OUT, "010_line_as_item", {
        "n_classes": R.N_CLASSES, "top_stems": R.TOP_STEMS, "shuffles": R.K, "n_sub": R.N_SUB, "n_null": R.N_NULL,
        "nulls": {k: list(v) for k, v in R.NULLS.items()}, "controls": list(R.CONTROLS),
        "interval": "real/controls: point + 2 x subsample quantile deviations; nulls: 2.5-97.5 % over seeds"},
        runs=json.loads(counts.to_json(orient="records")))


if __name__ == "__main__":
    summarize()
