"""Experiment 014: tables and plots from ``results/014_place_rules_inventory`` and the reference runs of 009-011.

Writes ``runs.csv`` (all augmented-generator rows), ``augmented_profile.csv`` (value, interval, SD and z of
Voynich against every null source: plain ts2020 of 009-011 and ts2020 + R1 / R1+R2 / R1+R2+R3),
``augmented_lists.csv`` (per metric and variant: catches up / gets closer / still off / breaks / unchanged),
``calibration.csv``, ``position_pareto.png``, ``seam_pareto.png``, ``augmented.png``, ``stability.png`` and
``meta.json``.

Intervals as in 009: real text = point + 2 x the 2.5 / 97.5 % quantiles of the deviations of the 80 % page
subsamples; null sources = mean and quantiles over the 30 seeds.  z = (Voynich - null mean) / null SD; the
verdict of a metric for a variant compares |z| of the plain ts2020 (z0) with |z| of the variant (z1):
``dozene`` |z0| >= 3 and |z1| < 3; ``priblizi`` |z1| < |z0| / 2 but still >= 3; ``zhorsi`` |z1| > 1.5 |z0| (both
>= 3); ``stale_ne`` otherwise with both >= 3; ``pokazi`` |z0| < 3 and |z1| >= 3; ``beze_zmeny`` both < 3.
The combined verdict over split and join is ``dozene`` only when both are.
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
RES = ROOT / "results"
SCALE = 2.0
Q = (0.025, 0.975)
Z_OK = 3.0
M009 = list(R.M9.METRICS)
M010 = ["bnd_edge_in", "bnd_edge_line", "bnd_edge_ratio_line", "bnd_affix_ratio_line", "bnd_stem_ratio_line",
        "end_enrich_max", "end_mg_last", "end_mg_other", "end_len_diff", "head_rest", "vert_excess_above",
        "tpl_gain_pos", "tpl_gain_prev", "jac_samelen_effect"]
M011 = ["a__suffix__pos__H_given", "a__suffix__pos__explained", "a__prefix__pos__explained",
        "a__suffix__hand__explained", "a__suffix__sec__explained", "a__suffix__all-page__explained",
        "a__nb_next_prefix__mi", "a__nb_prev_suffix__mi", "a__cons_page__excess", "a__share_empty_stem"]
CALIB = ["calib__H2_glyph", "calib__mean_len", "calib__n_types", "calib__hapax_share"]
FAMILY_METRICS = {"009": M009, "010": M010, "011": M011, "calib": CALIB}
# Voynich-specific properties (SYNTHESIS section 1 and the line conventions of 009/010); everything else is
# "other" (calibration targets, vocabulary, paradigm, order, section, locality) and may break
SPECIFIC = {
    ("009", "prefix__I_pos"): "pozicni_vazby_afixu", ("009", "suffix__I_pos"): "pozicni_vazby_afixu",
    ("011", "a__prefix__pos__explained"): "pozicni_vazby_afixu", ("011", "a__suffix__pos__explained"): "pozicni_vazby_afixu",
    ("009", "edge_mi_excess"): "sev_slov", ("010", "bnd_edge_in"): "sev_slov",
    ("011", "a__nb_next_prefix__mi"): "sev_slov", ("011", "a__nb_prev_suffix__mi"): "sev_slov",
    ("010", "end_enrich_max"): "konec_radku", ("010", "end_mg_last"): "konec_radku", ("010", "end_len_diff"): "konec_radku",
    ("009", "jsd_last_mid"): "konec_radku", ("009", "last_top5_excess"): "konec_radku",
    ("009", "jsd_first_mid"): "zacatek_radku", ("009", "first_top5_excess"): "zacatek_radku",
    ("010", "bnd_edge_ratio_line"): "uzavreny_radek", ("010", "bnd_affix_ratio_line"): "uzavreny_radek",
    ("010", "bnd_stem_ratio_line"): "uzavreny_radek",
    ("009", "I_first_rest"): "zahlavi_radku", ("010", "head_rest"): "zahlavi_radku",
}
SOURCES_NEW = [f"ts2020+{v}" for v in R.VARIANTS]
ORDER = ["ts2020"] + SOURCES_NEW


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #


def stat_rows(g: pd.DataFrame, metrics: list[str], family: str, source: str) -> list[dict]:
    rows = []
    kind = g["kind"].fillna("null") if "kind" in g.columns else pd.Series("null", index=g.index)
    point = g[kind == "point"]
    for m in metrics:
        if m not in g.columns:
            continue
        if len(point):
            v = pd.to_numeric(point[m], errors="coerce").iloc[0]
            s = pd.to_numeric(g.loc[kind == "sub", m], errors="coerce").dropna()
            dev = s - s.mean()
            lo, hi = (v + SCALE * dev.quantile(Q[0]), v + SCALE * dev.quantile(Q[1])) if len(s) > 1 else (np.nan,) * 2
            sd = SCALE * s.std(ddof=1) if len(s) > 1 else np.nan
            kind_out = "real"
        else:
            s = pd.to_numeric(g[m], errors="coerce").dropna()
            v, lo, hi, sd = s.mean(), s.quantile(Q[0]), s.quantile(Q[1]), s.std(ddof=1)
            kind_out = "null"
        rows.append({"source": source, "family": family, "metric": m, "kind": kind_out, "value": float(v),
                     "lo": float(lo), "hi": float(hi), "sd": float(sd), "n": int(len(s))})
    return rows


def reference_runs(family: str) -> pd.DataFrame:
    path = RES / {"009": "009_register_profile", "010": "010_line_as_item", "011": "011_suffix_variation"}[family] / "runs.csv"
    r = pd.read_csv(path, low_memory=False)
    r["kind"] = r["kind"].fillna("null")
    keep = list(R.AUG_TAGS) + [f"ts2020:{t}" for t in R.AUG_TAGS]
    return r[r["source"].isin(keep)]


def new_runs() -> pd.DataFrame:
    frames = []
    base = OUT / "runs"
    if base.exists():
        for d in sorted(base.iterdir()):
            if (d / "runs.csv").exists():
                frames.append(load_runs(d))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def calibration_reference() -> pd.DataFrame:
    """Calibration metrics of the real text (point) and of the plain ts2020 (30 seeds), same code as ``augment``."""
    rows = []
    for tag in R.AUG_TAGS:
        ctx = R.ctx_of(tag)
        rows.append({"source": tag, "kind": "point", **R.calib_metrics(ctx.real)})
        for seed in range(R.N_NULL):
            rows.append({"source": f"ts2020:{tag}", "kind": "null", "seed": seed, **R.calib_metrics(ctx.corpus("ts2020", seed))})
    return pd.DataFrame(rows)


def build_profile(runs_new: pd.DataFrame, calib: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fam in ("009", "010", "011"):
        ref = reference_runs(fam)
        for src, g in ref.groupby("source", sort=False):
            rows += stat_rows(g, FAMILY_METRICS[fam], fam, src)
    for src, g in calib.groupby("source", sort=False):
        rows += stat_rows(g, CALIB, "calib", src)
    if len(runs_new):
        for src, g in runs_new.groupby("source", sort=False):
            for fam, metrics in FAMILY_METRICS.items():
                pref = "" if fam == "calib" else f"p{fam}__"
                sub = g[[c for c in g.columns if c.startswith(pref) or c == "kind"]].rename(
                    columns=lambda c: c[len(pref):] if pref and c.startswith(pref) else c)
                rows += stat_rows(sub, metrics, fam, src)
    prof = pd.DataFrame(rows)
    prof["tag"] = prof["source"].str.split(":").str[-1]
    prof["null"] = np.where(prof["kind"] == "real", "real", prof["source"].str.split(":").str[0])
    prof["role"] = [SPECIFIC.get((f, m), "other") for f, m in zip(prof["family"], prof["metric"])]
    # z of Voynich against every null source of the same tag
    P = prof.set_index(["source", "family", "metric"])
    z, vv, rel = [], [], []
    for _, r in prof.iterrows():
        if r["kind"] == "real":
            z.append(np.nan)
            vv.append(r["value"])
            rel.append(0.0)
            continue
        key = (r["tag"], r["family"], r["metric"])
        v = float(P.loc[key, "value"]) if key in P.index else np.nan
        vv.append(v)
        z.append((v - r["value"]) / r["sd"] if r["sd"] and r["sd"] > 0 else np.nan)
        rel.append((r["value"] - v) / abs(v) if v else np.nan)
    prof["voynich_value"] = vv
    prof["z"] = z
    prof["rel_diff"] = rel
    return prof


# --------------------------------------------------------------------------- #
# lists
# --------------------------------------------------------------------------- #


def verdict(z0: float, z1: float) -> str:
    if not (np.isfinite(z0) and np.isfinite(z1)):
        return "nelze"
    a0, a1 = abs(z0), abs(z1)
    if a0 >= Z_OK and a1 < Z_OK:
        return "dozene"
    if a0 >= Z_OK and a1 >= Z_OK:
        if a1 < a0 / 2:
            return "priblizi"
        if a1 > 1.5 * a0:
            return "zhorsi"
        return "stale_ne"
    if a0 < Z_OK and a1 >= Z_OK:
        return "pokazi"
    return "beze_zmeny"


def combine(v_split: str, v_join: str) -> str:
    if v_split == v_join:
        return v_split
    s = {v_split, v_join}
    if "pokazi" in s or "zhorsi" in s:
        return "pokazi" if "pokazi" in s else "zhorsi"
    if "dozene" in s:
        return "dozene_jen_" + ("split" if v_split == "dozene" else "join")
    if "priblizi" in s:
        return "priblizi"
    return "stale_ne" if "stale_ne" in s else "nelze"


def lists(prof: pd.DataFrame) -> pd.DataFrame:
    P = prof.set_index(["null", "tag", "family", "metric"])
    rows = []
    for fam, metrics in FAMILY_METRICS.items():
        for m in metrics:
            for var in R.VARIANTS:
                row = {"family": fam, "metric": m, "role": SPECIFIC.get((fam, m), "other"), "variant": var}
                vs = {}
                for tag in R.AUG_TAGS:
                    try:
                        ref = P.loc[("ts2020", tag, fam, m)]
                        new = P.loc[(f"ts2020+{var}", tag, fam, m)]
                    except KeyError:
                        vs[tag] = "nelze"
                        continue
                    row[f"voynich_{tag}"] = ref["voynich_value"]
                    row[f"ts2020_{tag}"] = ref["value"]
                    row[f"aug_{tag}"] = new["value"]
                    row[f"z_ts2020_{tag}"] = ref["z"]
                    row[f"z_aug_{tag}"] = new["z"]
                    vs[tag] = verdict(ref["z"], new["z"])
                    row[f"verdict_{tag}"] = vs[tag]
                row["verdict"] = combine(vs.get("ZL_split", "nelze"), vs.get("ZL_join", "nelze"))
                rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #


def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        return None


def pareto_plots() -> None:
    plt = _plt()
    if plt is None:
        return
    pc = pd.read_csv(OUT / "position_contributions.csv")
    sc = pd.read_csv(OUT / "seam_contributions.csv")
    styles = {"ZL_split": ("C0", "-"), "GC_split": ("C1", "-"), "ts2020:ZL_split:s0": ("C0", "--"),
              "ctl:linnaeus_items": ("C2", ":"), "ctl:gibberish_native": ("C3", ":")}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharey=True)
    for ax, target in zip(axes.ravel(), R.TARGETS):
        for src, (col, ls) in styles.items():
            sub = pc[(pc["source"] == src) & (pc["target"] == target) & (pc["items"] == "by_item")
                     & (pc["variant"].isin(["a", "-"]))].sort_values("rank")
            if len(sub):
                ax.plot(sub["rank"], sub["cum_share"], color=col, ls=ls, marker="o", ms=3,
                        label=f"{src.replace('ctl:', '')} (I={sub['contribution'].sum():.3f} b)")
        for lv in (0.5, 0.8, 0.95):
            ax.axhline(lv, color="grey", lw=0.5)
        ax.set_title(f"I({target}; position class): cumulative share by {target}")
        ax.set_xlabel("items (sorted)")
        ax.set_xlim(0.5, 20.5)
        ax.legend(fontsize=7)
    axes[0, 0].set_ylabel("cumulative share")
    axes[1, 0].set_ylabel("cumulative share")
    fig.tight_layout()
    fig.savefig(OUT / "position_pareto.png", dpi=110)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    for ax, (dec, items, title) in zip(axes, (("seam_glyph", "pair", "pairs (last glyph, first glyph)"),
                                              ("seam_glyph", "by_last", "by last glyph"),
                                              ("seam_affix", "pair", "pairs (suffix, next prefix)"))):
        for src, (col, ls) in styles.items():
            sub = sc[(sc["source"] == src) & (sc["decomposition"] == dec) & (sc["items"] == items)
                     & (sc["variant"].isin(["a", "-"]))]
            if not len(sub):
                continue
            c = np.sort(sub["contribution"].to_numpy())[::-1]
            tot = c.sum()
            ax.plot(np.arange(1, len(c) + 1), np.cumsum(c) / tot, color=col, ls=ls, marker="o", ms=2.5,
                    label=f"{src.replace('ctl:', '')} (I={tot:.3f} b)")
        for lv in (0.5, 0.8, 0.95):
            ax.axhline(lv, color="grey", lw=0.5)
        ax.set_xlim(0.5, 25.5)
        ax.set_ylim(0, 1.15)
        ax.set_title(f"seam: {title}")
        ax.set_xlabel("items (sorted)")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("cumulative share of I")
    fig.tight_layout()
    fig.savefig(OUT / "seam_pareto.png", dpi=110)
    plt.close(fig)


def augmented_plot(prof: pd.DataFrame) -> None:
    plt = _plt()
    if plt is None:
        return
    metrics = [(f, m) for f, ms in FAMILY_METRICS.items() for m in ms]
    rows = [(n, t) for t in R.AUG_TAGS for n in ORDER]
    P = prof.set_index(["null", "tag", "family", "metric"])["z"]
    Z = np.full((len(rows), len(metrics)), np.nan)
    for i, (n, t) in enumerate(rows):
        for j, (f, m) in enumerate(metrics):
            Z[i, j] = P.get((n, t, f, m), np.nan)
    fig, ax = plt.subplots(figsize=(18, 4.8))
    im = ax.imshow(np.clip(Z, -8, 8), cmap="RdBu_r", vmin=-8, vmax=8, aspect="auto")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{n} {t}" for n, t in rows], fontsize=8)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m for _, m in metrics], rotation=80, ha="right", fontsize=6)
    for j, (f, m) in enumerate(metrics):
        if SPECIFIC.get((f, m)):
            ax.get_xticklabels()[j].set_color("C3")
    x = 0
    for f, ms in FAMILY_METRICS.items():
        ax.axvline(x - 0.5, color="k", lw=0.8)
        ax.text(x + len(ms) / 2 - 0.5, -0.9, f, ha="center", fontsize=8)
        x += len(ms)
    ax.axhline(len(ORDER) - 0.5, color="k", lw=0.8)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if np.isfinite(Z[i, j]) and abs(Z[i, j]) < Z_OK:
                ax.text(j, i, "·", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.015, label="z = (Voynich - null) / SD(null), clipped at 8; · = |z| < 3")
    ax.set_title("ts2020 and ts2020 + place rules vs Voynich ZL (red = Voynich higher; red labels = Voynich-specific)",
                 pad=18)
    fig.tight_layout()
    fig.savefig(OUT / "augmented.png", dpi=110)
    plt.close(fig)


def stability_plot() -> None:
    plt = _plt()
    path = OUT / "stability.csv"
    if plt is None or not path.exists():
        return
    st = pd.read_csv(path)
    st = st[st["variant"].isin(["a", "-"])]
    tables = list(dict.fromkeys(st["table"]))
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, tag in zip(axes, ("ZL_split", "GC_split")):
        sub = st[st["source"] == tag]
        for k, (gname, col) in enumerate((("language", "C0"), ("hand", "C1"), ("section", "C2"))):
            g = sub[sub["grouping"] == gname]
            agg = g.groupby("table")["ratio"].agg(["median", "min", "max"]).reindex(tables)
            y = np.arange(len(tables)) + (k - 1) * 0.25
            ax.errorbar(agg["median"], y, xerr=[agg["median"] - agg["min"], agg["max"] - agg["median"]], fmt="o",
                        color=col, label=f"{gname} (median, range over pairs)", ms=4, capsize=2)
        ax.axvline(1, color="k", lw=0.8)
        ax.set_xscale("log")
        ax.set_yticks(range(len(tables)))
        ax.set_yticklabels(tables, fontsize=8)
        ax.set_xlabel("JSD between groups / JSD between halves of one group")
        ax.set_title(tag)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "stability.png", dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def summarize() -> None:
    runs = new_runs()
    if len(runs):
        runs.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    calib = calibration_reference()
    calib.to_csv(OUT / "calibration.csv", index=False, float_format="%.6g")
    prof = build_profile(runs, calib)
    prof.to_csv(OUT / "augmented_profile.csv", index=False, float_format="%.6g")
    lst = lists(prof)
    lst.to_csv(OUT / "augmented_lists.csv", index=False, float_format="%.5g")
    pareto_plots()
    augmented_plot(prof)
    stability_plot()
    counts = runs.groupby("source").size().to_dict() if len(runs) else {}
    write_meta(OUT, "014_place_rules_inventory", {
        "variants": {k: list(v) for k, v in R.VARIANTS.items()}, "tags": list(R.AUG_TAGS), "n_null": R.N_NULL,
        "ts2020_config": R.M9.TS2020_CONFIG, "n_rep_stability": R.N_REP, "min_hand": R.MIN_HAND,
        "min_section": R.MIN_SECTION, "min_count": R.MIN_COUNT, "z_ok": Z_OK, "runs": counts,
        "verdict_rule": "z0 = ts2020, z1 = variant; dozene |z0|>=3 & |z1|<3; pokazi |z0|<3 & |z1|>=3"})
    print("written:", ", ".join(p.name for p in sorted(OUT.glob("*.*"))))
    if len(lst):
        print(lst.groupby(["variant", "verdict"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    summarize()
