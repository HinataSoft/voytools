"""Experiment 017: tables, lists and plots from ``results/017_conditional_seam_rewrite/runs``.

Writes ``runs.csv`` (all generator rows), ``real_bootstrap.csv`` (point and page-bootstrap rows of the real
text), ``augmented_profile.csv`` (source x metric: mean, SD and 2.5 / 97.5 % over seeds, Voynich value,
bootstrap SD of the Voynich value, ``z_seed`` = (Voynich - mean) / SD over seeds, ``z_boot`` = (Voynich - mean)
/ bootstrap SD, ``z_strict`` = the larger |z| of the two), ``metric_lists.csv`` (metric x set: verdict against
the plain ts2020 as in 014, on ``z_strict``), ``augmented_lists.csv`` (property x set: reached / broken /
unchanged, split and join), ``minimal_unreachable.csv`` (the smallest combinations of properties that no set
reaches at once, plus the properties each set reaches), ``rewrite_cost.csv``, ``template_check.csv``
(the same lists for the 005 template), ``selfcheck.csv`` (tables learned on ts2020 applied to ts2020),
``augmented.png``, ``rewrite_cost.png``, ``meta.json``.

Verdict of a metric for a set (014): z0 = plain generator, z1 = set, both ``z_strict``; ``dozene`` |z0| >= 3
and |z1| < 3; ``priblizi`` |z1| < |z0| / 2 but >= 3; ``zhorsi`` |z1| > 1.5 |z0| (both >= 3); ``stale_ne``
otherwise with both >= 3; ``pokazi`` |z0| < 3 and |z1| >= 3; ``beze_zmeny`` both < 3.  A property (a group of
metrics) is *reached* by a set when every metric has |z_strict| < 3 in split **and** join; its verdict uses the
largest |z| of its metrics.  The combined verdict over split and join is ``dozene`` only when both are.
"""

from __future__ import annotations

import itertools
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

import run as R  # noqa: E402  (experiments/017_conditional_seam_rewrite/run.py)

OUT = R.OUT
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
PREFIX = {"009": "p009__", "010": "p010__", "011": "p011__", "calib": ""}
# properties = groups of metrics; the first block are the Voynich-specific properties of SYNTHESIS section 1
# (without A/B), then the word-cohesion and calibration targets (004) that a rewrite may break
PROPERTIES = {
    "pozicni_vazba_prefixu": [("009", "prefix__I_pos"), ("011", "a__prefix__pos__explained")],
    "pozicni_vazba_sufixu": [("009", "suffix__I_pos"), ("011", "a__suffix__pos__explained")],
    "zacatek_radku": [("009", "jsd_first_mid"), ("009", "first_top5_excess")],
    "konec_radku": [("009", "jsd_last_mid"), ("009", "last_top5_excess"), ("010", "end_mg_last"), ("010", "end_len_diff")],
    "konec_radku_17x": [("010", "end_enrich_max")],
    "sev_glyfovy": [("009", "edge_mi_excess"), ("010", "bnd_edge_in")],
    "sev_sufix_prefix": [("011", "a__nb_next_prefix__mi"), ("011", "a__nb_prev_suffix__mi")],
    "uzavreny_radek": [("010", "bnd_edge_ratio_line"), ("010", "bnd_affix_ratio_line"), ("010", "bnd_stem_ratio_line")],
    "zahlavi_radku": [("009", "I_first_rest"), ("010", "head_rest")],
    "vazba_prefix_sufix": [("009", "mi_prefix_suffix")],
    "sekcni_vazba_prefixu": [("009", "prefix__I_sec"), ("009", "prefix__I_sec_g_pos")],
    "sirka_paradigmatu": [("009", "mean_suffixes_per_stem"), ("009", "H_suffix_given_stem"), ("009", "frac_stems_multi_suffix")],
    "H_sufix_kmen_vysvetleno": [("011", "a__suffix__all-page__explained")],
    "cmi_excess": [("009", "cmi_excess")],
    "typy": [("calib", "calib__n_types")],
    "hapax": [("calib", "calib__hapax_share")],
    "H2_glyfu": [("calib", "calib__H2_glyph")],
    "delka": [("calib", "calib__mean_len")],
    "zipf": [("009", "zipf_slope")],
}
SPECIFIC = ("pozicni_vazba_prefixu", "pozicni_vazba_sufixu", "zacatek_radku", "konec_radku", "konec_radku_17x",
            "sev_glyfovy", "sev_sufix_prefix", "uzavreny_radek", "zahlavi_radku")
COHESION = ("vazba_prefix_sufix", "sekcni_vazba_prefixu", "sirka_paradigmatu", "H_sufix_kmen_vysvetleno", "typy", "hapax")
TARGET_PROPERTIES = SPECIFIC + COHESION          # candidates of the minimal unreachable combination
ROLE = {m: p for p, ms in PROPERTIES.items() for m in ms}
SETS = list(R.SETS)
COST_COLS = ["share_tokens_changed", "share_backoff", "share_rejected", "n_types_base", "n_types", "n_types_new",
             "n_tokens_new_types", "n_conditional_candidates", "n_rejected_total", "n_cand_total"]


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #


def all_runs() -> pd.DataFrame:
    frames = []
    base = OUT / "runs"
    if base.exists():
        for d in sorted(base.iterdir()):
            if (d / "runs.csv").exists():
                frames.append(load_runs(d))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def metric_col(fam: str, m: str) -> str:
    return PREFIX[fam] + m


def stat_rows(g: pd.DataFrame, source: str, real: pd.DataFrame | None = None) -> list[dict]:
    """Mean / SD / quantiles of every metric of one source (``real`` = point + bootstrap rows)."""
    rows = []
    tag = source.split(":")[-1]
    for fam, metrics in FAMILY_METRICS.items():
        for m in metrics:
            col = metric_col(fam, m)
            if col not in g.columns:
                continue
            if real is not None:
                point = pd.to_numeric(g.loc[g["kind"] == "point", col], errors="coerce")
                boot = pd.to_numeric(g.loc[g["kind"] == "boot", col], errors="coerce").dropna()
                v = float(point.iloc[0]) if len(point) else np.nan
                rows.append({"source": source, "tag": tag, "null": "real", "set": "real", "family": fam, "metric": m,
                             "kind": "real", "value": v, "sd": float(boot.std(ddof=1)) if len(boot) > 1 else np.nan,
                             "lo": float(v - (boot.quantile(0.5) - boot.quantile(Q[0]))) if len(boot) > 1 else np.nan,
                             "hi": float(v + (boot.quantile(Q[1]) - boot.quantile(0.5))) if len(boot) > 1 else np.nan,
                             "boot_mean": float(boot.mean()) if len(boot) else np.nan, "n": int(len(boot))})
            else:
                s = pd.to_numeric(g[col], errors="coerce").dropna()
                rows.append({"source": source, "tag": tag, "null": source.split(":")[0], "set": str(g["set"].iloc[0]),
                             "family": fam, "metric": m, "kind": "null", "value": float(s.mean()),
                             "sd": float(s.std(ddof=1)) if len(s) > 1 else np.nan, "lo": float(s.quantile(Q[0])),
                             "hi": float(s.quantile(Q[1])), "boot_mean": np.nan, "n": int(len(s))})
    return rows


def build_profile(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    real = runs[runs["generator"] == "real"]
    for tag, g in real.groupby("tag", sort=False):
        rows += stat_rows(g, tag, real=g)
    for src, g in runs[runs["generator"] != "real"].groupby("source", sort=False):
        rows += stat_rows(g, src)
    prof = pd.DataFrame(rows)
    prof["role"] = [ROLE.get((f, m), "other") for f, m in zip(prof["family"], prof["metric"])]
    P = prof[prof["kind"] == "real"].set_index(["tag", "family", "metric"])
    vv, sb, zs, zb = [], [], [], []
    for _, r in prof.iterrows():
        key = (r["tag"], r["family"], r["metric"])
        v = float(P.loc[key, "value"]) if key in P.index else np.nan
        s_boot = float(P.loc[key, "sd"]) if key in P.index else np.nan
        vv.append(v)
        sb.append(s_boot)
        if r["kind"] == "real":
            zs.append(np.nan)
            zb.append(np.nan)
        else:
            zs.append((v - r["value"]) / r["sd"] if r["sd"] and r["sd"] > 0 else np.nan)
            zb.append((v - r["value"]) / s_boot if s_boot and s_boot > 0 else np.nan)
    prof["voynich_value"] = vv
    prof["sd_boot"] = sb
    prof["z_seed"] = zs
    prof["z_boot"] = zb
    za, zbb = np.asarray(zs, dtype=float), np.asarray(zb, dtype=float)
    strict = np.where(np.isnan(zbb), za, np.where(np.isnan(za), zbb, np.where(np.abs(zbb) > np.abs(za), zbb, za)))
    prof["z_strict"] = strict
    prof["stricter"] = np.where(prof["kind"] == "real", "", np.where(np.abs(zbb) > np.abs(za), "boot", "seed"))
    prof["rel_diff"] = (prof["value"] - prof["voynich_value"]) / prof["voynich_value"].abs().replace(0, np.nan)
    return prof


# --------------------------------------------------------------------------- #
# verdicts
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


def metric_lists(prof: pd.DataFrame, generator: str, sets: list[str], tags: tuple) -> pd.DataFrame:
    P = prof[prof["kind"] == "null"].set_index(["null", "tag", "family", "metric"])
    base_null = generator
    rows = []
    for fam, metrics in FAMILY_METRICS.items():
        for m in metrics:
            for s in sets:
                if s == "base":
                    continue
                row = {"generator": generator, "family": fam, "metric": m, "role": ROLE.get((fam, m), "other"), "set": s}
                vs = {}
                for tag in tags:
                    try:
                        ref = P.loc[(base_null, tag, fam, m)]
                        new = P.loc[(f"{generator}+{s}", tag, fam, m)]
                    except KeyError:
                        vs[tag] = "nelze"
                        continue
                    row[f"voynich_{tag}"] = ref["voynich_value"]
                    row[f"base_{tag}"] = ref["value"]
                    row[f"aug_{tag}"] = new["value"]
                    row[f"z_base_{tag}"] = ref["z_strict"]
                    row[f"z_aug_{tag}"] = new["z_strict"]
                    row[f"z_aug_seed_{tag}"] = new["z_seed"]
                    row[f"z_aug_boot_{tag}"] = new["z_boot"]
                    vs[tag] = verdict(ref["z_strict"], new["z_strict"])
                    row[f"verdict_{tag}"] = vs[tag]
                row["verdict"] = combine(vs.get(tags[0], "nelze"), vs.get(tags[-1], "nelze")) if len(tags) > 1 \
                    else vs.get(tags[0], "nelze")
                rows.append(row)
    return pd.DataFrame(rows)


def property_lists(prof: pd.DataFrame, generator: str, sets: list[str], tags: tuple) -> pd.DataFrame:
    """Property x set: reached (all metrics |z_strict| < 3 in every tag), verdict against the plain generator."""
    P = prof[prof["kind"] == "null"].set_index(["null", "tag", "family", "metric"])
    rows = []
    for prop, metrics in PROPERTIES.items():
        for s in sets:
            null = generator if s == "base" else f"{generator}+{s}"
            row = {"generator": generator, "property": prop, "group": "specific" if prop in SPECIFIC else
                   ("cohesion" if prop in COHESION else "calibration"), "set": s, "n_metrics": len(metrics)}
            vs, reached = {}, True
            for tag in tags:
                z0, z1, worst = [], [], ""
                for fam, m in metrics:
                    try:
                        z1.append(float(P.loc[(null, tag, fam, m), "z_strict"]))
                        z0.append(float(P.loc[(generator, tag, fam, m), "z_strict"]))
                    except KeyError:
                        z1.append(np.nan)
                        z0.append(np.nan)
                z0a, z1a = np.asarray(z0), np.asarray(z1)
                ok = bool(np.all(np.abs(z1a) < Z_OK)) if np.isfinite(z1a).all() and len(z1a) else False
                reached &= ok
                i = int(np.nanargmax(np.abs(z1a))) if np.isfinite(z1a).any() else 0
                row[f"ok_{tag}"] = ok
                row[f"max_abs_z_{tag}"] = float(np.nanmax(np.abs(z1a))) if np.isfinite(z1a).any() else np.nan
                row[f"max_abs_z_base_{tag}"] = float(np.nanmax(np.abs(z0a))) if np.isfinite(z0a).any() else np.nan
                row[f"worst_metric_{tag}"] = metrics[i][1]
                row[f"n_ok_{tag}"] = int((np.abs(z1a) < Z_OK).sum())
                vs[tag] = "beze_zmeny" if s == "base" and ok else ("stale_ne" if s == "base" else
                                                                    verdict(row[f"max_abs_z_base_{tag}"], row[f"max_abs_z_{tag}"]))
                row[f"verdict_{tag}"] = vs[tag]
            row["reached"] = reached
            row["verdict"] = combine(vs.get(tags[0], "nelze"), vs.get(tags[-1], "nelze")) if len(tags) > 1 else vs.get(tags[0], "nelze")
            if s == "base":
                row["verdict"] = "dosazeno" if reached else "nedosazeno"
            else:
                row["verdict"] = {"dozene": "dosazeno", "beze_zmeny": "dosazeno" if reached else "beze_zmeny",
                                  "pokazi": "pokazeno"}.get(row["verdict"], row["verdict"])
            rows.append(row)
    return pd.DataFrame(rows)


def minimal_unreachable(plist: pd.DataFrame, sets: list[str], max_size: int = 4) -> pd.DataFrame:
    """Smallest combinations of target properties that no set reaches at once (every proper subset is
    reached by some set), plus one row per set with the properties it reaches."""
    reached = {s: set(plist.loc[(plist["set"] == s) & plist["reached"], "property"]) for s in sets}
    rows = [{"kind": "reached_by_set", "set": s, "size": len(reached[s] & set(TARGET_PROPERTIES)),
             "combination": " + ".join(p for p in TARGET_PROPERTIES if p in reached[s]),
             "not_reached": " + ".join(p for p in TARGET_PROPERTIES if p not in reached[s])} for s in sets]
    props = [p for p in TARGET_PROPERTIES if p in set(plist["property"])]

    def covered(comb) -> list[str]:
        return [s for s in sets if set(comb) <= reached[s]]

    found = []
    for k in range(1, max_size + 1):
        for comb in itertools.combinations(props, k):
            if covered(comb):
                continue
            if any(not covered(sub) for sub in itertools.combinations(comb, k - 1) if k > 1):
                continue  # a proper subset is already unreachable: not minimal
            found.append(comb)
        if found and k >= 2:
            break  # the smallest sizes are enough for the report
    for comb in found:
        subs = {" + ".join(sub): ",".join(covered(sub)) for sub in itertools.combinations(comb, len(comb) - 1)} if len(comb) > 1 else {}
        rows.append({"kind": "minimal_unreachable", "set": "", "size": len(comb), "combination": " + ".join(comb),
                     "not_reached": "", "subsets_reached_by": " | ".join(f"{a}: {b}" for a, b in subs.items())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# cost, template, self-check
# --------------------------------------------------------------------------- #


def rewrite_cost(runs: pd.DataFrame) -> pd.DataFrame:
    g = runs[runs["generator"] != "real"]
    rows = []
    for (src, s, tag, gen), d in g.groupby(["source", "set", "tag", "generator"], sort=False):
        row = {"source": src, "generator": gen, "set": s, "tag": tag, "n_seeds": len(d), "rules": d["rules"].iloc[0],
               "licence": bool(d["licence"].iloc[0])}
        for c in COST_COLS:
            if c in d.columns:
                v = pd.to_numeric(d[c], errors="coerce")
                row[f"{c}_mean"] = float(v.mean())
                row[f"{c}_sd"] = float(v.std(ddof=1)) if len(v) > 1 else np.nan
        for r in R.RULE_NAMES:
            for k in ("n_cand", "n_changed", "n_rejected"):
                c = f"{k}_{r}"
                if c in d.columns:
                    row[f"{c}_mean"] = float(pd.to_numeric(d[c], errors="coerce").mean())
        for lv in R.LEVELS:
            c = f"support_{lv.replace(':', '_')}"
            if c in d.columns:
                row[f"{c}_mean"] = float(pd.to_numeric(d[c], errors="coerce").mean())
        row["share_types_new_mean"] = row.get("n_types_new_mean", np.nan) / row["n_types_mean"] if row.get("n_types_mean") else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def selfcheck(prof: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    P = prof[prof["kind"] == "null"].set_index(["null", "tag", "family", "metric"])
    rows = []
    for tag in R.TAGS:
        sc = runs[runs["source"] == f"selfcheck+R12pA:{tag}"]
        if sc.empty:
            continue
        for fam, metrics in FAMILY_METRICS.items():
            for m in metrics:
                try:
                    ref = P.loc[("ts2020", tag, fam, m)]
                except KeyError:
                    continue
                col = metric_col(fam, m)
                s = pd.to_numeric(sc[col], errors="coerce").dropna()
                base = runs[(runs["source"] == f"ts2020:{tag}")].set_index("seed")
                paired = pd.to_numeric(sc.set_index("seed")[col], errors="coerce") - pd.to_numeric(base[col], errors="coerce").reindex(sc["seed"])
                paired = paired.dropna()
                rows.append({"tag": tag, "family": fam, "metric": m, "role": ROLE.get((fam, m), "other"),
                             "ts2020_mean": ref["value"], "ts2020_sd": ref["sd"], "selfcheck_mean": float(s.mean()),
                             "selfcheck_sd": float(s.std(ddof=1)) if len(s) > 1 else np.nan, "n": int(len(s)),
                             "z_mean": (float(s.mean()) - ref["value"]) / ref["sd"] if ref["sd"] > 0 else np.nan,
                             "paired_diff_mean": float(paired.mean()) if len(paired) else np.nan,
                             "paired_diff_sd": float(paired.std(ddof=1)) if len(paired) > 1 else np.nan,
                             "paired_t": float(paired.mean() / (paired.std(ddof=1) / np.sqrt(len(paired)))) if len(paired) > 1 and paired.std(ddof=1) > 0 else np.nan,
                             "voynich_value": ref["voynich_value"], "z_ts2020_strict": ref["z_strict"]})
        if len(sc):
            rows.append({"tag": tag, "family": "cost", "metric": "share_tokens_changed", "role": "cost",
                         "selfcheck_mean": float(pd.to_numeric(sc["share_tokens_changed"]).mean()), "n": int(len(sc))})
            rows.append({"tag": tag, "family": "cost", "metric": "share_backoff", "role": "cost",
                         "selfcheck_mean": float(pd.to_numeric(sc["share_backoff"]).mean()), "n": int(len(sc))})
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


def augmented_plot(prof: pd.DataFrame) -> None:
    plt = _plt()
    if plt is None:
        return
    metrics = [(f, m) for f, ms in FAMILY_METRICS.items() for m in ms]
    nulls = ["ts2020"] + [f"ts2020+{s}" for s in SETS if s != "base"]
    rows = [(n, t) for t in R.TAGS for n in nulls if ((prof["null"] == n) & (prof["tag"] == t)).any()]
    P = prof[prof["kind"] == "null"].set_index(["null", "tag", "family", "metric"])["z_strict"]
    Z = np.full((len(rows), len(metrics)), np.nan)
    for i, (n, t) in enumerate(rows):
        for j, (f, m) in enumerate(metrics):
            Z[i, j] = P.get((n, t, f, m), np.nan)
    fig, ax = plt.subplots(figsize=(18, 0.42 * len(rows) + 2.6))
    im = ax.imshow(np.clip(Z, -8, 8), cmap="RdBu_r", vmin=-8, vmax=8, aspect="auto")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{n} {t}" for n, t in rows], fontsize=8)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m for _, m in metrics], rotation=80, ha="right", fontsize=6)
    for j, (f, m) in enumerate(metrics):
        role = ROLE.get((f, m))
        if role in SPECIFIC:
            ax.get_xticklabels()[j].set_color("C3")
        elif role in COHESION:
            ax.get_xticklabels()[j].set_color("C2")
    x = 0
    for f, ms in FAMILY_METRICS.items():
        ax.axvline(x - 0.5, color="k", lw=0.8)
        ax.text(x + len(ms) / 2 - 0.5, -0.9, f, ha="center", fontsize=8)
        x += len(ms)
    ax.axhline(len(nulls) - 0.5, color="k", lw=0.8)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if np.isfinite(Z[i, j]) and abs(Z[i, j]) < Z_OK:
                ax.text(j, i, "·", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.015, label="z_strict = (Voynich - null) / SD (stricter of seeds / bootstrap), clipped at 8; · = |z| < 3")
    ax.set_title("ts2020 + conditional place rules vs Voynich ZL (red labels = Voynich-specific, green = word cohesion / vocabulary)", pad=18)
    fig.tight_layout()
    fig.savefig(OUT / "augmented.png", dpi=110)
    plt.close(fig)


def cost_plot(cost: pd.DataFrame, prof: pd.DataFrame) -> None:
    plt = _plt()
    if plt is None or cost.empty:
        return
    c = cost[cost["generator"] == "ts2020"]
    sets = [s for s in SETS if s in set(c["set"])]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    w = 0.38
    for k, tag in enumerate(R.TAGS):
        d = c[c["tag"] == tag].set_index("set").reindex(sets)
        x = np.arange(len(sets)) + (k - 0.5) * w
        axes[0].bar(x, d["share_tokens_changed_mean"], w, label=f"{tag}: tokens rewritten", color=f"C{k}")
        axes[0].bar(x, d["share_backoff_mean"].fillna(0) * d["n_conditional_candidates_mean"].fillna(0) / d["n_tokens_base_mean"] if "n_tokens_base_mean" in d else 0,
                    w, color="k", alpha=0.35)
        axes[1].bar(x, d["share_rejected_mean"].fillna(0), w, color=f"C{k}", label=tag)
        axes[2].bar(x, d["n_types_mean"], w, color=f"C{k}", label=f"{tag}: types")
        axes[2].bar(x, d["n_types_new_mean"], w, color="k", alpha=0.35, label="of which new vs plain ts2020" if k == 0 else None)
    real = prof[(prof["kind"] == "real") & (prof["metric"] == "calib__n_types")].set_index("tag")["value"]
    for k, tag in enumerate(R.TAGS):
        if tag in real.index:
            axes[2].axhline(real[tag], color=f"C{k}", ls="--", lw=1, label=f"Voynich {tag}")
    for ax, title in zip(axes, ("share of tokens rewritten", "share of rewrites rejected by the type licence",
                                "types (mean over seeds)")):
        ax.set_xticks(range(len(sets)))
        ax.set_xticklabels(sets, rotation=30, ha="right", fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "rewrite_cost.png", dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def summarize() -> None:
    runs = all_runs()
    if runs.empty:
        print("no runs")
        return
    gen = runs[runs["generator"] != "real"]
    gen.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    runs[runs["generator"] == "real"].to_csv(OUT / "real_bootstrap.csv", index=False, float_format="%.6g")
    prof = build_profile(runs)
    prof.to_csv(OUT / "augmented_profile.csv", index=False, float_format="%.6g")
    ts_sets = [s for s in SETS if (prof["null"] == ("ts2020" if s == "base" else f"ts2020+{s}")).any()]
    ml = metric_lists(prof, "ts2020", ts_sets, R.TAGS)
    ml.to_csv(OUT / "metric_lists.csv", index=False, float_format="%.5g")
    pl = property_lists(prof, "ts2020", ts_sets, R.TAGS)
    pl.to_csv(OUT / "augmented_lists.csv", index=False, float_format="%.5g")
    mu = minimal_unreachable(pl, ts_sets)
    mu.to_csv(OUT / "minimal_unreachable.csv", index=False)
    cost = rewrite_cost(runs)
    cost.to_csv(OUT / "rewrite_cost.csv", index=False, float_format="%.5g")
    tm_sets = [s for s in R.TEMPLATE_SETS if (prof["null"] == ("tmpl" if s == "base" else f"tmpl+{s}")).any()]
    if tm_sets:
        tl = property_lists(prof, "tmpl", tm_sets, R.TEMPLATE_TAGS)
        tm = metric_lists(prof, "tmpl", tm_sets, R.TEMPLATE_TAGS)
        tl["level"] = "property"
        tm["level"] = "metric"
        pd.concat([tl, tm], ignore_index=True).to_csv(OUT / "template_check.csv", index=False, float_format="%.5g")
        minimal_unreachable(tl, tm_sets).to_csv(OUT / "template_minimal_unreachable.csv", index=False)
    sc = selfcheck(prof, runs)
    sc.to_csv(OUT / "selfcheck.csv", index=False, float_format="%.5g")
    augmented_plot(prof)
    cost_plot(cost, prof)
    counts = runs.groupby("source").size().to_dict()
    write_meta(OUT, "017_conditional_seam_rewrite", {
        "sets": {k: {"rules": list(v[0]), "licence": v[1]} for k, v in R.SETS.items()}, "tags": list(R.TAGS),
        "template_tags": list(R.TEMPLATE_TAGS), "n_null": R.N_NULL, "n_boot": R.N_BOOT, "min_cell": R.MIN_CELL,
        "min_context": R.MIN_CONTEXT, "end_ratio": R.END_RATIO, "licence_min_count": R.LICENCE_MIN,
        "max_tries": R.MAX_TRIES, "ts2020_config": R.M9.TS2020_CONFIG, "z_ok": Z_OK, "runs": counts,
        "properties": {k: [f"{f}:{m}" for f, m in v] for k, v in PROPERTIES.items()},
        "verdict_rule": "z_strict = larger |z| of z_seed (SD over seeds) and z_boot (page-bootstrap SD of the real text); "
                        "property reached = all its metrics |z_strict| < 3 in split and join"})
    print("written:", ", ".join(p.name for p in sorted(OUT.glob("*.*"))))
    if len(pl):
        print(pl.pivot(index="property", columns="set", values="verdict").reindex(columns=ts_sets).to_string())
    if len(mu):
        print(mu.to_string())
    if len(sc):
        bad = sc[(sc["role"] != "cost") & (sc["z_mean"].abs() >= 1)]
        print(f"selfcheck: {len(bad)} metrics with |z_mean| >= 1 of {int((sc['role'] != 'cost').sum())}")


if __name__ == "__main__":
    summarize()
