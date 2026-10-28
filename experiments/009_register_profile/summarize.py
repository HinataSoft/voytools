"""Experiment 009: tables and plot from ``results/009_register_profile/runs/*``.

Writes ``runs.csv`` (all rows), ``profile.csv``, ``profile_z.csv``, ``distances.csv``, ``ranking.csv``,
``lists.csv``, ``checks.csv``, ``profile_heatmap.png`` and ``meta.json``.
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
SCALE = float(np.sqrt(R.SUB_FRAC / (1 - R.SUB_FRAC)))  # = 2 (007: SD x 2)
Q = (0.025, 0.975)
VOY = list(R.REAL_TAGS)
PROSE = ["caesar_prose", "celsus_prose", "pliny_prose", "plinyP_prose", "culpeper_prose", "italian_prose"]
ITEMS = ["linnaeus_items", "culpeper_items"]
CATALOGUE_PROSE = ["linnaeus_prose"]
NULL_KINDS = list(R.NULLS)


def source_class(src: str) -> str:
    if src in VOY:
        return "voynich"
    if src.startswith("ctl:"):
        n = src[4:]
        return "prose" if n in PROSE else "items" if n in ITEMS else "catalogue_prose"
    return "shuffle" if src.startswith("line:") else "generator"


def group_of(metric: str) -> str:
    return next(g for g, ms in R.GROUPS.items() if metric in ms)


def load_all() -> pd.DataFrame:
    frames = []
    for d in sorted((OUT / "runs").iterdir()):
        if (d / "runs.csv").exists():
            frames.append(load_runs(d))
    runs = pd.concat(frames, ignore_index=True)
    runs["kind"] = runs["kind"].fillna("null")  # read_csv parses the label 'null' as NA
    return runs


def summarize_source(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src, g in runs.groupby("source", sort=False):
        point = g[g["kind"] == "point"]
        sub = g[g["kind"] == "sub"]
        nul = g[g["kind"] == "null"]
        for m in R.METRICS:
            if len(point):
                v = float(point[m].iloc[0])
                s = sub[m].dropna()
                dev = s - s.mean()
                row = {"value": v, "lo": v + SCALE * float(dev.quantile(Q[0])), "hi": v + SCALE * float(dev.quantile(Q[1])),
                       "sd": SCALE * float(s.std(ddof=1)), "n": len(s)}
            else:
                s = nul[m].dropna()
                row = {"value": float(s.mean()), "lo": float(s.quantile(Q[0])), "hi": float(s.quantile(Q[1])),
                       "sd": float(s.std(ddof=1)), "n": len(s)}
            rows.append({"metric": m, "group": group_of(m), "source": src, "class": source_class(src), **row})
    return pd.DataFrame(rows)


def comparison_set(prof: pd.DataFrame, tag: str) -> list[str]:
    """Sources compared with Voynich ``tag``: all controls and the nulls generated from ``tag``."""
    srcs = prof["source"].unique().tolist()
    return [s for s in srcs if s.startswith("ctl:") or (":" in s and s.split(":")[1] == tag)]


def z_table(prof: pd.DataFrame) -> pd.DataFrame:
    P = prof.set_index(["source", "metric"])
    rows = []
    for tag in VOY:
        if tag not in prof["source"].values:
            continue
        for src in comparison_set(prof, tag):
            for m in R.METRICS:
                v, s = P.loc[(tag, m)], P.loc[(src, m)]
                sd = s["sd"] if s["sd"] > 0 else np.nan
                rows.append({"voynich": tag, "source": src, "class": source_class(src), "metric": m,
                             "group": group_of(m), "voynich_value": v["value"], "source_value": s["value"],
                             "source_sd": s["sd"], "z": (v["value"] - s["value"]) / sd if sd == sd else np.nan})
    return pd.DataFrame(rows)


def distances(prof: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """RMS of standardized differences Voynich - source; all metrics, leave-one-group-out, single group.

    Standardization: per metric the SD of the point values over the comparison set (controls, the
    tag's nulls and the Voynich tag itself).
    """
    W = prof.pivot(index="source", columns="metric", values="value")
    rows = []
    for tag in VOY:
        if tag not in W.index:
            continue
        comp = comparison_set(prof, tag)
        sub = W.loc[comp + [tag], R.METRICS]
        scale = sub.std(ddof=1).replace(0, np.nan)
        Z = (sub - sub.loc[tag]) / scale
        for src in comp:
            z = Z.loc[src]
            rows.append({"voynich": tag, "source": src, "class": source_class(src), "metrics": "all",
                         "distance": float(np.sqrt(np.nanmean(z ** 2)))})
            for g, ms in R.GROUPS.items():
                keep = [m for m in R.METRICS if m not in ms]
                rows.append({"voynich": tag, "source": src, "class": source_class(src), "metrics": f"without_{g}",
                             "distance": float(np.sqrt(np.nanmean(z[keep] ** 2)))})
                rows.append({"voynich": tag, "source": src, "class": source_class(src), "metrics": f"only_{g}",
                             "distance": float(np.sqrt(np.nanmean(z[ms] ** 2)))})
    d = pd.DataFrame(rows)
    d["rank"] = d.groupby(["voynich", "metrics"])["distance"].rank(method="min").astype(int)
    rank = (d.sort_values(["voynich", "metrics", "rank"])
            .groupby(["voynich", "metrics"])
            .agg(order=("source", lambda s: " < ".join(x.replace("ctl:", "") for x in s)),
                 nearest=("source", "first"), nearest_class=("class", "first"))
            .reset_index())
    cls = d.groupby(["voynich", "metrics", "class"])["distance"].min().unstack("class").reset_index()
    rank = rank.merge(cls, on=["voynich", "metrics"])
    return d.sort_values(["voynich", "metrics", "rank"]), rank


def lists(prof: pd.DataFrame) -> pd.DataFrame:
    """(a) inside the range of the real texts, (b) outside real texts and all generators, (c) outside real, inside
    a generator (ts2020, template, auto-copy).

    The within-line shuffle is reported (``in_shuffle``) but does not decide (b)/(c): it is the real text with
    the token order changed, so every type-level or page-level metric equals the real value by construction.
    """
    P = prof.set_index(["source", "metric"])
    ctl = [s for s in prof["source"].unique() if s.startswith("ctl:")]
    rows = []
    for tag in VOY:
        if tag not in prof["source"].values:
            continue
        nulls = [s for s in comparison_set(prof, tag) if not s.startswith("ctl:")]
        for m in R.METRICS:
            v = P.loc[(tag, m)]
            lo = min(P.loc[(c, m)]["lo"] for c in ctl)
            hi = max(P.loc[(c, m)]["hi"] for c in ctl)
            vlo = min(P.loc[(c, m)]["value"] for c in ctl)
            vhi = max(P.loc[(c, m)]["value"] for c in ctl)
            inside_real = lo <= v["value"] <= hi
            in_nulls = [n for n in nulls if P.loc[(n, m)]["lo"] <= v["value"] <= P.loc[(n, m)]["hi"]]
            in_gen = [n for n in in_nulls if not n.startswith("line:")]
            lab = "a" if inside_real else ("c" if in_gen else "b")
            side = "" if inside_real else ("above" if v["value"] > hi else "below")
            rows.append({"voynich": tag, "metric": m, "group": group_of(m), "value": v["value"], "lo": v["lo"],
                         "hi": v["hi"], "real_min": vlo, "real_max": vhi, "real_lo": lo, "real_hi": hi,
                         "list": lab, "side": side, "outside_all_generators": not in_gen,
                         "in_shuffle": any(n.startswith("line:") for n in in_nulls),
                         "generators_containing": " ".join(n.split(":")[0] for n in in_gen)})
    out = pd.DataFrame(rows)
    main = out[out["voynich"].isin(["ZL_split", "GC_split"])].groupby("metric")["list"].agg(
        lambda s: s.iloc[0] if s.nunique() == 1 else "mixed")
    out["list_ZL_GC"] = out["metric"].map(main)
    return out


def checks(prof: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    P = prof.set_index(["source", "metric"])
    rows = []
    for src in [s for s in prof["source"].unique() if s.startswith("ctl:")]:
        h = P.loc[(src, "hapax_share")]["value"]
        rows.append({"check": "hapax_share > 0.9", "source": src, "metric": "hapax_share", "value": h,
                     "voynich_ZL_split": P.loc[("ZL_split", "hapax_share")]["value"], "flag": bool(h > 0.9)})
        if src[4:] in PROSE + CATALOGUE_PROSE:
            for y in R.M7.TARGETS:
                m = f"{y}__I_pos"
                v, ref = P.loc[(src, m)]["value"], P.loc[("ZL_split", m)]["value"]
                rows.append({"check": "I(Y;pos) of mechanically broken prose > Voynich", "source": src, "metric": m,
                             "value": v, "voynich_ZL_split": ref, "flag": bool(v > ref)})
    return pd.DataFrame(rows)


def heatmap(prof: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    tag = "ZL_split"
    comp = comparison_set(prof, tag)
    order = [tag] + [s for s in comp if s.startswith("ctl:")] + [s for s in comp if not s.startswith("ctl:")]
    W = prof.pivot(index="source", columns="metric", values="value").loc[order, R.METRICS]
    Z = (W - W.loc[tag]) / W.std(ddof=1).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(16, 6.5))
    im = ax.imshow(Z.clip(-3, 3).to_numpy(dtype=float), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([s.replace("ctl:", "") for s in order], fontsize=8)
    ax.set_xticks(range(len(R.METRICS)))
    ax.set_xticklabels(R.METRICS, rotation=75, ha="right", fontsize=7)
    x = 0
    for g, ms in R.GROUPS.items():
        ax.axvline(x - 0.5, color="k", lw=0.8)
        ax.text(x + len(ms) / 2 - 0.5, -0.9, g, ha="center", fontsize=8)
        x += len(ms)
    ax.axhline(0.5, color="k", lw=0.8)
    ax.axhline(0.5 + sum(s.startswith("ctl:") for s in comp), color="k", lw=0.8)
    fig.colorbar(im, ax=ax, fraction=0.02, label="(source - ZL split) / SD over sources, clipped at 3")
    ax.set_title("Register profile: difference from Voynich ZL split (red = source higher)", pad=18)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def summarize() -> None:
    runs = load_all()
    runs.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    prof = summarize_source(runs)
    prof.to_csv(OUT / "profile.csv", index=False, float_format="%.6g")
    z_table(prof).to_csv(OUT / "profile_z.csv", index=False, float_format="%.4g")
    d, rank = distances(prof)
    d.to_csv(OUT / "distances.csv", index=False, float_format="%.4g")
    rank.to_csv(OUT / "ranking.csv", index=False, float_format="%.4g")
    lists(prof).to_csv(OUT / "lists.csv", index=False, float_format="%.5g")
    checks(prof, runs).to_csv(OUT / "checks.csv", index=False, float_format="%.5g")
    heatmap(prof, OUT / "profile_heatmap.png")
    counts = runs.groupby(["source", "kind"]).size().reset_index(name="n")
    write_meta(OUT, "009_register_profile", {
        "groups": R.GROUPS, "shuffles": R.SHUF, "n_sub": R.N_SUB, "sub_frac": R.SUB_FRAC, "n_null": R.N_NULL,
        "k_sections": R.K_SECTIONS, "template": R.TEMPLATE_TAG, "nulls": {k: list(v) for k, v in R.NULLS.items()},
        "ts2020_config": R.TS2020_CONFIG, "interval": "real/controls: point + 2 x subsample quantile deviations; "
        "nulls: 2.5-97.5 % over seeds"}, runs=json.loads(counts.to_json(orient="records")))
    print(rank[rank["metrics"] == "all"].to_string())


if __name__ == "__main__":
    summarize()
