"""Aggregation of experiment 020 runs into the spec's tables and figures (``run.py summarize``).

Reads ``results/020_paragraph_as_unit/runs/<combo>/<source>/runs.csv`` and writes, in the results folder:
``runs.csv`` (all runs except the 200-seed sources, which stay in ``runs/``), ``summary.csv`` (long table:
combo x source x variant x metric with value and 2.5 / 97.5 % interval), ``line_pair_sharing.csv``,
``boundary_step.csv``, ``within_par_decay.csv``, ``first_line.csv``, ``page_boundary.csv``,
``form_by_boundary.csv``, ``meta.json`` and the figures ``boundary_step.png`` / ``page_boundary.png``.

Value of a source = the real run (``real``; interval = page bootstrap ``real_boot``) or the mean over seeds
(layout nulls, generators; interval = 2.5 / 97.5 % quantiles over seeds).  ``p`` = share of null seeds at least
as extreme as the real value (one-sided, in the direction of the real value's sign), ``(1 + #) / (1 + n)``.
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

OUT = ROOT / "results" / "020_paragraph_as_unit"
RUNS = OUT / "runs"
Q = (0.025, 0.975)
COL = re.compile(r"^(st|st50|tok|ng)__(.+)$")
SOURCE_ORDER = ["real", "shuffle_bounds", "perm_lines", "selfcitation_21", "template_005"]
COMBO_ORDER = ["ZL_split", "ZL_join", "GC_split", "GC_join", "ctl:celsus_chap", "ctl:pliny_chap", "ctl:celsus_sent40",
               "ctl:pliny_sent40", "ctl:linnaeus_par", "ctl:pliny_prose", "ctl:gibberish_pooled", "ctl:gibberish_native"]
BIG = ("real_boot", "shuffle_bounds", "perm_lines")


def load_all() -> pd.DataFrame:
    frames = [pd.read_csv(p) for p in sorted(RUNS.glob("*/*/runs.csv"))]
    return pd.concat(frames, ignore_index=True)


def _q(s: pd.Series, q: float) -> float:
    s = s.dropna()
    return float(s.quantile(q)) if len(s) else float("nan")


def long_summary(runs: pd.DataFrame) -> pd.DataFrame:
    """One row per combo x source x variant x metric: value, lo, hi, n_seeds, p (nulls vs the real value)."""
    cols = [c for c in runs.columns if COL.match(c)]
    rows = []
    for (combo, source), g in runs.groupby(["combo", "source"], sort=False):
        if source == "real_boot":
            continue
        real = runs[(runs["combo"] == combo) & (runs["source"] == "real")]
        boot = runs[(runs["combo"] == combo) & (runs["source"] == "real_boot")] if source == "real" else None
        for c in cols:
            var, metric = COL.match(c).groups()
            vals = g[c].astype(float)
            if source == "real":
                v = float(vals.iloc[0])
                lo, hi = (_q(boot[c].astype(float), Q[0]), _q(boot[c].astype(float), Q[1])) if boot is not None and len(boot) else (np.nan, np.nan)
                p = np.nan
            else:
                v = float(vals.mean()) if vals.notna().any() else np.nan
                lo, hi = _q(vals, Q[0]), _q(vals, Q[1])
                p = np.nan
                if len(real) and np.isfinite(real[c].astype(float).iloc[0]) and vals.notna().any():
                    rv = float(real[c].astype(float).iloc[0])
                    vv = vals.dropna()
                    ext = (vv >= rv).sum() if rv >= vv.mean() else (vv <= rv).sum()
                    p = (1 + ext) / (1 + len(vv))
            rows.append((combo, source, var, metric, v, lo, hi, int(vals.notna().sum()), p))
    out = pd.DataFrame(rows, columns=["combo", "source", "variant", "metric", "value", "lo", "hi", "n_seeds", "p"])
    out["combo"] = pd.Categorical(out["combo"], [c for c in COMBO_ORDER if c in set(out["combo"])] +
                                  sorted(set(out["combo"]) - set(COMBO_ORDER)))
    out["source"] = pd.Categorical(out["source"], SOURCE_ORDER + sorted(set(out["source"]) - set(SOURCE_ORDER)))
    return out.sort_values(["combo", "source", "variant", "metric"]).reset_index(drop=True)


def _pick(summary: pd.DataFrame, pattern: str) -> pd.DataFrame:
    rx = re.compile(pattern)
    sub = summary[summary["metric"].str.match(rx)].copy()
    parts = sub["metric"].str.extract(rx)
    for k, col in enumerate(parts.columns):
        sub[f"g{k}"] = parts[col].to_numpy()
    return sub


def table_sharing(summary: pd.DataFrame) -> pd.DataFrame:
    """line_pair_sharing.csv: source x variant x expectation x d x same/diff, s with interval and n pairs."""
    sub = _pick(summary, r"^share(nf)?_(p|g|n)__d(\d+)_(same|diff)$")
    sub["set"] = np.where(sub["g0"] == "nf", "no_first_lines", "all")
    sub["expect"] = sub["g1"].map({"p": "page", "g": "global", "n": "n"})
    sub["d"] = sub["g2"].astype(int)
    sub["same"] = sub["g3"]
    n = sub[sub["expect"] == "n"][["combo", "source", "variant", "set", "d", "same", "value"]].rename(columns={"value": "n_pairs"})
    s = sub[sub["expect"] != "n"][["combo", "source", "variant", "set", "expect", "d", "same", "value", "lo", "hi", "n_seeds", "p"]]
    out = s.merge(n, on=["combo", "source", "variant", "set", "d", "same"], how="left").rename(columns={"value": "s"})
    return out.sort_values(["combo", "source", "variant", "set", "expect", "d", "same"]).reset_index(drop=True)


def table_step(summary: pd.DataFrame) -> pd.DataFrame:
    sub = _pick(summary, r"^step(nf)?_(p|g)__d(\d+)$")
    sub["set"] = np.where(sub["g0"] == "nf", "no_first_lines", "all")
    sub["expect"] = sub["g1"].map({"p": "page", "g": "global"})
    sub["d"] = sub["g2"].astype(int)
    out = sub[["combo", "source", "variant", "set", "expect", "d", "value", "lo", "hi", "n_seeds", "p"]].rename(columns={"value": "delta"})
    return out.sort_values(["combo", "source", "variant", "set", "expect", "d"]).reset_index(drop=True)


def table_decay(summary: pd.DataFrame) -> pd.DataFrame:
    sub = _pick(summary, r"^decay_(p|g|n)__(d\d+|slope)$")
    sub["expect"] = sub["g0"].map({"p": "page", "g": "global", "n": "n"})
    n = sub[sub["expect"] == "n"][["combo", "source", "variant", "g1", "value"]].rename(columns={"value": "n_pairs"})
    s = sub[sub["expect"] != "n"][["combo", "source", "variant", "expect", "g1", "value", "lo", "hi", "n_seeds", "p"]]
    out = s.merge(n, on=["combo", "source", "variant", "g1"], how="left").rename(columns={"g1": "d", "value": "s"})
    out["d"] = out["d"].str.replace("d", "", regex=False)
    return out.sort_values(["combo", "source", "variant", "expect", "d"]).reset_index(drop=True)


def table_first(summary: pd.DataFrame) -> pd.DataFrame:
    sub = _pick(summary, r"^first_(p|g|n)__(first_own|first_next|mid_own|mid_next)_d(\d+)$")
    sub["expect"] = sub["g0"].map({"p": "page", "g": "global", "n": "n"})
    sub["d"] = sub["g2"].astype(int)
    n = sub[sub["expect"] == "n"][["combo", "source", "variant", "g1", "d", "value"]].rename(columns={"value": "n_pairs"})
    s = sub[sub["expect"] != "n"][["combo", "source", "variant", "expect", "g1", "d", "value", "lo", "hi", "n_seeds", "p"]]
    out = s.merge(n, on=["combo", "source", "variant", "g1", "d"], how="left").rename(columns={"g1": "kind", "value": "s"})
    return out.sort_values(["combo", "source", "variant", "expect", "kind", "d"]).reset_index(drop=True)


def table_page(summary: pd.DataFrame) -> pd.DataFrame:
    sub = _pick(summary, r"^page(_n|_ratio)?__(across_all|across_same|across_diff|across_open|within|within_cross)_d(\d+)$")
    sub["kind"] = sub["g0"].fillna("").map({"": "s", "_n": "n", "_ratio": "ratio"})
    sub["d"] = sub["g2"].astype(int)
    base = sub[sub["kind"] == "s"][["combo", "source", "variant", "g1", "d", "value", "lo", "hi", "n_seeds", "p"]].rename(columns={"value": "s"})
    n = sub[sub["kind"] == "n"][["combo", "source", "variant", "g1", "d", "value"]].rename(columns={"value": "n_pairs"})
    r = sub[sub["kind"] == "ratio"][["combo", "source", "variant", "g1", "d", "value", "lo", "hi", "p"]].rename(
        columns={"value": "ratio_vs_within", "lo": "ratio_lo", "hi": "ratio_hi", "p": "ratio_p"})
    out = base.merge(n, on=["combo", "source", "variant", "g1", "d"], how="left").merge(r, on=["combo", "source", "variant", "g1", "d"], how="left")
    return out.rename(columns={"g1": "group"}).sort_values(["combo", "source", "variant", "group", "d"]).reset_index(drop=True)


def table_form(summary: pd.DataFrame) -> pd.DataFrame:
    sub = _pick(summary, r"^form(_n)?__(excess|ident)?_?(same|diff)_d(\d+)$")
    sub["kind"] = np.where(sub["g0"].fillna("") == "_n", "n_occ", sub["g1"])
    sub["same"] = sub["g2"]
    sub["d"] = sub["g3"].astype(int)
    keys = ["combo", "source", "variant", "same", "d"]
    ex = sub[sub["kind"] == "excess"][keys + ["value", "lo", "hi", "n_seeds", "p"]].rename(columns={"value": "excess"})
    idn = sub[sub["kind"] == "ident"][keys + ["value"]].rename(columns={"value": "ident"})
    n = sub[sub["kind"] == "n_occ"][keys + ["value"]].rename(columns={"value": "n_occ"})
    out = ex.merge(idn, on=keys, how="left").merge(n, on=keys, how="left")
    return out.sort_values(keys).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

PLOT = [("ZL_split", "real", "#1f77b4", "ZL real"), ("GC_split", "real", "#17becf", "GC real"),
        ("ZL_split", "perm_lines", "#7f7f7f", "ZL perm. řádků"), ("ZL_split", "shuffle_bounds", "#bcbd22", "ZL přelosování hranic"),
        ("ZL_split", "selfcitation_21", "#ff7f0e", "ts2020"), ("ZL_split", "template_005", "#2ca02c", "šablona 005"),
        ("ctl:celsus_chap", "real", "#d62728", "Celsus kapitoly"), ("ctl:pliny_chap", "real", "#e377c2", "Plinius kapitoly"),
        ("ctl:pliny_sent40", "real", "#f7b6d2", "Plinius věty ~40"), ("ctl:linnaeus_par", "real", "#8c564b", "Linnaeus položky"),
        ("ctl:pliny_prose", "real", "#9467bd", "Plinius, hranice ZL"), ("ctl:gibberish_pooled", "real", "#c5b0d5", "blábol, hranice ZL")]


def _val(t: pd.DataFrame, combo: str, source: str, **kw) -> pd.DataFrame:
    sub = t[(t["combo"].astype(str) == combo) & (t["source"].astype(str) == source)]
    for k, v in kw.items():
        sub = sub[sub[k] == v]
    return sub


def plot_step(step: pd.DataFrame, share: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    for k, (combo, source, color, label) in enumerate(PLOT):
        sub = _val(step, combo, source, variant="st", set="all", expect="page").sort_values("d")
        if sub.empty:
            continue
        x = sub["d"].to_numpy() + (k - len(PLOT) / 2) * 0.03
        ax.errorbar(x, sub["delta"], yerr=[sub["delta"] - sub["lo"], sub["hi"] - sub["delta"]], fmt="o-", ms=3, lw=1,
                    color=color, label=label, capsize=2)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("vzdálenost řádků d")
    ax.set_ylabel("Δ(d) = s_same − s_diff (očekávání ze stránky)")
    ax.set_title("Skok sdílení kmenů na hranici odstavce")
    ax.legend(fontsize=6.5, ncol=2)
    ax = axes[1]
    for combo, source, color, label in PLOT[:2] + PLOT[4:8]:
        for same, ls in (("same", "-"), ("diff", "--")):
            sub = _val(share, combo, source, variant="st", set="all", expect="page", same=same).sort_values("d")
            if sub.empty:
                continue
            ax.plot(sub["d"], sub["s"], ls, color=color, lw=1, marker="o", ms=2.5, label=f"{label} {same}")
    ax.set_xlabel("vzdálenost řádků d")
    ax.set_ylabel("s(d) = pozorované / očekávané")
    ax.set_title("Sdílení v témž odstavci (plně) a přes hranici (čárkovaně)")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_page(page: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.2))
    rows = PLOT + [("ctl:gibberish_native", "real", "#636363", "blábol vlastní stránky")]
    for k, (combo, source, color, label) in enumerate(rows):
        sub = _val(page, combo, source, variant="st", group="across_all").sort_values("d")
        if sub.empty or sub["ratio_vs_within"].isna().all():
            continue
        x = sub["d"].to_numpy() + (k - len(rows) / 2) * 0.03
        ax.errorbar(x, sub["ratio_vs_within"], yerr=[sub["ratio_vs_within"] - sub["ratio_lo"], sub["ratio_hi"] - sub["ratio_vs_within"]],
                    fmt="o-", ms=3, lw=1, color=color, label=label, capsize=2)
    ax.axhline(1, color="k", lw=0.6)
    ax.set_xlabel("vzdálenost řádků přes hranici stránky d")
    ax.set_ylabel("sdílení přes hranici / uvnitř stránky (glob. očekávání)")
    ax.set_title("Hranice stránky: pokračování (≈ 1) nebo zlom (< 1)")
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def summarize() -> None:
    runs = load_all()
    OUT.mkdir(parents=True, exist_ok=True)
    small = runs[~runs["source"].isin(BIG)]
    small.to_csv(OUT / "runs.csv", index=False)
    summary = long_summary(runs)
    summary.round(6).to_csv(OUT / "summary.csv", index=False)
    share = table_sharing(summary)
    share.round(5).to_csv(OUT / "line_pair_sharing.csv", index=False)
    step = table_step(summary)
    step.round(5).to_csv(OUT / "boundary_step.csv", index=False)
    table_decay(summary).round(5).to_csv(OUT / "within_par_decay.csv", index=False)
    table_first(summary).round(5).to_csv(OUT / "first_line.csv", index=False)
    page = table_page(summary)
    page.round(5).to_csv(OUT / "page_boundary.csv", index=False)
    table_form(summary).round(5).to_csv(OUT / "form_by_boundary.csv", index=False)
    try:
        plot_step(step, share, OUT / "boundary_step.png")
        plot_page(page, OUT / "page_boundary.png")
    except ImportError:
        pass
    counts = runs.groupby(["combo", "source"]).size()
    metas = {}
    for p in sorted(RUNS.glob("*/*/meta.json")):
        m = json.loads(p.read_text(encoding="utf-8"))
        metas[f"{p.parent.parent.name}/{p.parent.name}"] = m.get("params", {})
    write_meta(OUT, "020_paragraph_as_unit",
               params={"runs": {f"{c}/{s}": int(n) for (c, s), n in counts.items()}, "quantiles": list(Q),
                       "per_source_params": metas},
               n_rows=int(len(runs)))
    print(f"{len(runs)} runs, {len(summary)} summary rows -> {OUT}")


if __name__ == "__main__":
    summarize()
