"""Experiment 007: tables and plots from the run folders (``run.py summarize``).

Intervals of the real text come from the 200 page subsamples (80 % of the pages without
replacement): the deviations of the subsample estimates from their mean are scaled by
``sqrt(m / (n - m)) = 2`` (subsampling without replacement; the spread of an 80 % subsample
around the full-text value is half the sampling spread of the full text) and added to the
full-text point estimate.  Null intervals are the 2.5 / 97.5 % quantiles over seeds.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run import (  # noqa: E402  (same folder)
    CORRECTIONS, DMAX, DRIFT_AXES, FILES, N_NULL, N_SUB, OUT, POLICIES, SUB_FRAC, TARGETS, TERMS, UNITS,
)
from voytools import load_runs  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402

SCALE = float(np.sqrt(SUB_FRAC / (1 - SUB_FRAC)))
TAGS = [f"{k}_{p}" for k in FILES for p in POLICIES]
NULL_SOURCES = {"decomp": ("autocopy",), "mi": ("posclass", "line", "autocopy"), "drift": ("autocopy",)}
Q = (0.025, 0.975)


def real_interval(point: float, sub: pd.Series) -> tuple[float, float]:
    s = sub.dropna()
    if not len(s) or not np.isfinite(point):
        return float("nan"), float("nan")
    dev = s - s.mean()
    return point + SCALE * float(dev.quantile(Q[0])), point + SCALE * float(dev.quantile(Q[1]))


def null_summary(s: pd.Series) -> dict:
    s = s.dropna()
    if not len(s):
        return {"value": float("nan"), "lo": float("nan"), "hi": float("nan"), "sd": float("nan")}
    return {"value": float(s.mean()), "lo": float(s.quantile(Q[0])), "hi": float(s.quantile(Q[1])),
            "sd": float(s.std(ddof=1)) if len(s) > 1 else float("nan")}


# --------------------------------------------------------------------------- #
# metric 1
# --------------------------------------------------------------------------- #


def decomposition_table() -> pd.DataFrame:
    rows = []
    for tag in TAGS:
        base = OUT / "decomp" / tag
        if not (base / "point.csv").exists():
            continue
        point = pd.read_csv(base / "point.csv")
        sub = load_runs(base / "real_sub")
        nulls = {s: load_runs(base / s) for s in NULL_SOURCES["decomp"]}
        for _, r in point.iterrows():
            y, term = r["Y"], r["term"]
            col = f"{y}__{term}__corr"
            lo, hi = real_interval(r["corr"], sub[col])
            hs = sub[f"{y}__H"] - sub[f"{y}__H"].mean()
            row = {"dataset": tag, "Y": y, "term": term, "H": r["H"],
                   "H_lo": r["H"] + SCALE * hs.quantile(Q[0]), "H_hi": r["H"] + SCALE * hs.quantile(Q[1]),
                   "raw": r["raw"], "bias": r["bias"], "corr": r["corr"], "corr_lo": lo, "corr_hi": hi,
                   "share": r["corr"] / r["H"], "share_lo": lo / r["H"], "share_hi": hi / r["H"],
                   "n_sub": int(sub[col].notna().sum())}
            for nl in ("within_page", "page_labels"):
                for c in ("mean", "q975", "p"):
                    row[f"perm_{nl}_{c}"] = r.get(f"{nl}_{c}", np.nan)
            for s, df in nulls.items():
                if df.empty:
                    continue
                ns = null_summary(df[col])
                row.update({f"{s}_corr": ns["value"], f"{s}_corr_lo": ns["lo"], f"{s}_corr_hi": ns["hi"],
                            f"{s}_H": float(df[f"{y}__H"].mean()), f"{s}_share": ns["value"] / df[f"{y}__H"].mean(),
                            f"{s}_n": int(df[col].notna().sum())})
            # the permutation null that can test this term: within-page for positional X, page labels otherwise
            row["p_perm"] = row["perm_within_page_p"] if term.startswith("I_pos") else row["perm_page_labels_p"]
            row["flag_I_gt_H"] = bool(r["raw"] > r["H"] + 1e-9)
            row["flag_corr_lt_-0.01"] = bool(r["corr"] < -0.01)
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# metric 2
# --------------------------------------------------------------------------- #


def corrected(df: pd.DataFrame, u: str, corr: str) -> np.ndarray:
    """Matrix runs x d of corrected MI (corr = 'raw' for the uncorrected curve)."""
    cols = [f"{u}_d{d}_raw" for d in range(1, DMAX + 1)]
    raw = df[cols].to_numpy(dtype=float)
    if corr == "raw":
        return raw
    sh = df[[f"{u}_d{d}_sh_{corr}" for d in range(1, DMAX + 1)]].to_numpy(dtype=float)
    return raw - sh


def half_decay(curve: np.ndarray) -> float:
    """d at which the curve first falls to 50 % of its d = 1 value (linear interpolation); 10.5 = not within 10."""
    c1 = curve[0]
    if not np.isfinite(c1) or c1 <= 0:
        return float("nan")
    half = 0.5 * c1
    for i in range(1, len(curve)):
        if curve[i] <= half:
            prev = curve[i - 1]
            return float(i + (prev - half) / (prev - curve[i])) if prev != curve[i] else float(i + 1)
    return DMAX + 0.5


def curve_scalars(m: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame({"mi_d1": m[:, 0], "d_half": [half_decay(r) for r in m], "area_2_10": m[:, 1:].sum(axis=1),
                         "mi_d2": m[:, 1], "mi_d3": m[:, 2]})


def mi_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    long, summ = [], []
    corrs = ["raw"] + list(CORRECTIONS)
    for tag in TAGS:
        base = OUT / "mi" / tag
        point = load_runs(base / "real_point")
        if point.empty:
            continue
        sub = load_runs(base / "real_sub")
        nulls = {s: load_runs(base / s) for s in NULL_SOURCES["mi"]}
        for u in UNITS:
            for corr in corrs:
                p = corrected(point, u, corr)[0]
                sm = corrected(sub, u, corr) if not sub.empty else np.full((0, DMAX), np.nan)
                for d in range(1, DMAX + 1):
                    lo, hi = real_interval(p[d - 1], pd.Series(sm[:, d - 1]))
                    long.append({"dataset": tag, "U": u, "d": d, "source": "real", "correction": corr,
                                 "mi": p[d - 1], "lo": lo, "hi": hi, "n": len(sm)})
                ps = curve_scalars(p[None, :]).iloc[0]
                ss = curve_scalars(sm)
                row = {"dataset": tag, "U": u, "source": "real", "correction": corr, "n": len(sm)}
                for k in ps.index:
                    lo, hi = real_interval(ps[k], ss[k])
                    row.update({k: ps[k], f"{k}_lo": lo, f"{k}_hi": hi})
                lo_d = [real_interval(p[d], pd.Series(sm[:, d]))[0] for d in range(DMAX)]
                row["max_d_lo_pos"] = int(next((d for d in range(DMAX) if not lo_d[d] > 0), DMAX))
                row["d_half_not_reached_share"] = float((ss["d_half"] > DMAX).mean()) if len(ss) else np.nan
                summ.append(row)
                for s, df in nulls.items():
                    if df.empty:
                        continue
                    m = corrected(df, u, corr)
                    for d in range(1, DMAX + 1):
                        ns = null_summary(pd.Series(m[:, d - 1]))
                        long.append({"dataset": tag, "U": u, "d": d, "source": s, "correction": corr,
                                     "mi": ns["value"], "lo": ns["lo"], "hi": ns["hi"], "sd": ns["sd"], "n": len(m)})
                    cs = curve_scalars(m)
                    row = {"dataset": tag, "U": u, "source": s, "correction": corr, "n": len(m)}
                    for k in cs.columns:
                        ns = null_summary(cs[k])
                        row.update({k: ns["value"], f"{k}_lo": ns["lo"], f"{k}_hi": ns["hi"]})
                    mean_c = m.mean(axis=0)
                    row["max_d_lo_pos"] = int(next((d for d in range(DMAX)
                                                    if not np.quantile(m[:, d], Q[0]) > 0), DMAX))
                    row["d_half_of_mean_curve"] = half_decay(mean_c)
                    summ.append(row)
    long_df = pd.DataFrame(long)
    if len(long_df):
        long_df["flag_corr_lt_-0.01"] = (long_df["correction"] != "raw") & (long_df["mi"] < -0.01)
    return long_df, pd.DataFrame(summ)


# --------------------------------------------------------------------------- #
# metric 3
# --------------------------------------------------------------------------- #


def drift_table() -> pd.DataFrame:
    rows = []
    for tag in TAGS:
        base = OUT / "drift" / tag
        point = load_runs(base / "real_point")
        if point.empty:
            continue
        sub = load_runs(base / "real_sub")
        nulls = {s: load_runs(base / s) for s in NULL_SOURCES["drift"]}
        cols = [c for c in point.columns if c.split("__")[0] in DRIFT_AXES]
        for c in cols:
            axis, b, *rest = c.split("__")
            glyph = rest[0] if rest else b
            kind = "share" if rest else b
            bin_ = int(b) if rest else -1
            v = float(point[c].iloc[0])
            lo, hi = real_interval(v, sub[c]) if c in sub else (np.nan, np.nan)
            row = {"dataset": tag, "axis": axis, "bin": bin_, "glyph": glyph if rest else "", "kind": kind,
                   "real": v, "real_lo": lo, "real_hi": hi}
            for s, df in nulls.items():
                if not df.empty and c in df:
                    ns = null_summary(df[c])
                    row.update({f"{s}": ns["value"], f"{s}_lo": ns["lo"], f"{s}_hi": ns["hi"]})
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #

C = {"real": "#2a78d6", "autocopy": "#eb6834", "line": "#1baf7a", "posclass": "#8a8984"}
TXT2 = "#52514e"
LABEL = {"real": "reál", "autocopy": "autocopy 86", "line": "míchání v řádku", "posclass": "míchání v poz. třídě"}


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#b5b4ae")
    ax.tick_params(colors=TXT2, labelsize=7)
    ax.grid(axis="y", color="#e6e5e0", linewidth=0.6)
    ax.set_axisbelow(True)


def plot_decomposition(tab: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    terms = ["I_pos", "I_sec", "I_hand", "I_sec_g_pos", "I_pos_g_sec"]
    tlab = ["pozice", "sekce", "písař", "sekce | poz.", "poz. | sekce"]
    fig, axes = plt.subplots(len(TAGS), len(TARGETS), figsize=(12, 2.3 * len(TAGS)), sharey="col")
    for i, tag in enumerate(TAGS):
        for j, y in enumerate(TARGETS):
            ax = axes[i, j]
            _style(ax)
            sub = tab[(tab.dataset == tag) & (tab.Y == y)].set_index("term").reindex(terms)
            x = np.arange(len(terms))
            ax.bar(x - 0.2, sub["share"], 0.38, color=C["real"], label=LABEL["real"], edgecolor="white", linewidth=1)
            ax.errorbar(x - 0.2, sub["share"], yerr=[sub["share"] - sub["share_lo"], sub["share_hi"] - sub["share"]],
                        fmt="none", ecolor="#1d3f6e", elinewidth=0.8, capsize=2)
            if "autocopy_share" in sub:
                ac_lo = sub["autocopy_corr_lo"] / sub["autocopy_H"]
                ac_hi = sub["autocopy_corr_hi"] / sub["autocopy_H"]
                ax.bar(x + 0.2, sub["autocopy_share"], 0.38, color=C["autocopy"], label=LABEL["autocopy"],
                       edgecolor="white", linewidth=1)
                ax.errorbar(x + 0.2, sub["autocopy_share"], yerr=[sub["autocopy_share"] - ac_lo,
                                                                  ac_hi - sub["autocopy_share"]],
                            fmt="none", ecolor="#7a2f10", elinewidth=0.8, capsize=2)
            ax.set_xticks(x)
            ax.set_xticklabels(tlab, rotation=35, ha="right", fontsize=7)
            ax.axhline(0, color="#b5b4ae", linewidth=0.6)
            if i == 0:
                ax.set_title(y, fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{tag}\npodíl H(Y)", fontsize=8, color=TXT2)
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle("007: podíl H(Y) vysvětlený pozicí, sekcí a písařem (korigovaná MI / H(Y), 95% intervaly)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_mi(long: pd.DataFrame, path: Path, correction: str = "pc") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(UNITS), len(TAGS), figsize=(13, 2.4 * len(UNITS)), sharex=True)
    ulab = {"edge": "poslední→první glyf", "bpe": "BPE-64 jednotka", "affix": "afixová třída", "stem": "kmen ≥20"}
    for i, u in enumerate(UNITS):
        for j, tag in enumerate(TAGS):
            ax = axes[i, j]
            _style(ax)
            for s in ("posclass", "line", "autocopy", "real"):
                d = long[(long.dataset == tag) & (long.U == u) & (long.source == s) & (long.correction == correction)]
                if d.empty:
                    continue
                d = d.sort_values("d")
                ax.fill_between(d["d"], d["lo"], d["hi"], color=C[s], alpha=0.18, linewidth=0)
                ax.plot(d["d"], d["mi"], color=C[s], linewidth=2 if s == "real" else 1.4, marker="o", markersize=3,
                        label=LABEL[s])
            ax.axhline(0, color="#8a8984", linewidth=0.7)
            if i == 0:
                ax.set_title(tag, fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{ulab[u]}\nkorig. MI (bit)", fontsize=8, color=TXT2)
            if i == len(UNITS) - 1:
                ax.set_xlabel("vzdálenost d (tokeny, uvnitř odstavce)", fontsize=8, color=TXT2)
            ax.set_xticks(range(1, DMAX + 1))
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"007: MI(U_t; U_t+d) minus průměr míchání v poziční třídě odstavce ({correction}); "
                 "95% intervaly", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------- #


def summarize() -> None:
    dec = decomposition_table()
    dec.to_csv(OUT / "info_decomposition.csv", index=False, float_format="%.6g")
    long, summ = mi_tables()
    long.to_csv(OUT / "mi_distance.csv", index=False, float_format="%.6g")
    summ.to_csv(OUT / "mi_summary.csv", index=False, float_format="%.6g")
    drift = drift_table()
    drift.to_csv(OUT / "drift.csv", index=False, float_format="%.6g")
    # index of the per-run files (the rows themselves stay in the subfolders)
    idx = []
    for path in sorted(OUT.glob("*/*/*/runs.csv")):
        r = pd.read_csv(path)
        phase, tag, source = path.relative_to(OUT).parts[:3]
        idx.append({"phase": phase, "dataset": tag, "source": source, "path": str(path.relative_to(OUT)),
                    "n_runs": len(r), "seed_min": int(r["seed"].min()), "seed_max": int(r["seed"].max()),
                    "n_columns": r.shape[1], "elapsed_s_total": round(float(r["elapsed_s"].sum()), 1)})
    pd.DataFrame(idx).to_csv(OUT / "runs.csv", index=False)
    flags = {
        "decomp_I_gt_H": int(dec["flag_I_gt_H"].sum()) if len(dec) else 0,
        "decomp_corr_lt_-0.01": int(dec["flag_corr_lt_-0.01"].sum()) if len(dec) else 0,
        "mi_corr_lt_-0.01": int(long["flag_corr_lt_-0.01"].sum()) if len(long) else 0,
    }
    write_meta(OUT, "007_position_vs_content", {"n_sub": N_SUB, "sub_frac": SUB_FRAC, "n_null": N_NULL,
                                                "interval_scale": SCALE, "terms": [t[0] for t in TERMS],
                                                "units": list(UNITS), "corrections": list(CORRECTIONS)},
               flags=flags, files=[p.name for p in sorted(OUT.glob("*.csv"))])
    print(json.dumps(flags))
    try:
        if len(dec):
            plot_decomposition(dec, OUT / "decomposition.png")
        if len(long):
            plot_mi(long, OUT / "mi_distance.png", "pc")
    except ImportError:
        print("matplotlib not available: plots skipped")
