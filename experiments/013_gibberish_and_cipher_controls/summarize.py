"""Experiment 013: tables and plot from ``results/013_gibberish_and_cipher_controls/runs/*`` and the results of 009-011.

Writes ``runs.csv`` (all 013 rows), ``profile.csv`` (009 profile, 010 line metrics and 011 decomposition of
the 013 sources plus the reference sources taken over from ``results/009..011``), ``distances.csv``,
``feature_lists.csv``, ``profile_heatmap.png`` and ``meta.json``.

Intervals as in 009: point + 2 x the 2.5 / 97.5 % quantiles of the deviations of the 80 % page subsamples.

Feature verdicts (``feature_lists.csv``): for every Voynich-specific property of ``results/SYNTHESIS.md``
section 1, f = (source - R) / (V - R), where V is the Voynich value and R the mean of the meaningful prose
controls, both at the source's size (gibberish: ``voynich_crop`` and ``caesar_crop``/``celsus_crop``;
Naibbe: ``ZL_split`` and Caesar/Celsus/Pliny of 009).  f >= 0.5: "ano" (the source has the property),
0.2 <= f < 0.5: "částečně", f < 0.2: "ne"; f is also given at the ends of the source's interval.
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
SCALE = 2.0  # sqrt(0.8 / 0.2)
Q = (0.025, 0.975)
M009 = R.M9.METRICS
M010 = ["bnd_edge_ratio_line", "bnd_affix_ratio_line", "end_enrich_max", "end_len_diff", "head_rest",
        "vert_excess_above", "tpl_gain_pos"]
M011 = ["a__suffix__pos__H_given", "a__suffix__pos__explained", "a__prefix__pos__explained",
        "a__suffix__hand__explained", "a__suffix__sec__explained", "a__suffix__all-page__explained",
        "a__nb_next_prefix__mi", "a__nb_prev_suffix__mi", "a__cons_page__excess"]
NEW = [f"ctl:{s}" for s in R.SOURCES]
PROSE_FULL = ["ctl:caesar_prose", "ctl:celsus_prose", "ctl:pliny_prose"]
PROSE_CROP = ["ctl:caesar_crop", "ctl:celsus_crop"]
REF_SOURCES = ["ZL_split", "GC_split", "ts2020:ZL_split", "tmpl:ZL_split", "ac86:ZL_split", *PROSE_FULL,
               "ctl:culpeper_prose", "ctl:linnaeus_items", "ctl:linnaeus_prose", "ctl:italian_prose"]

# Voynich-specific properties of results/SYNTHESIS.md section 1 -> (family, metric, note)
FEATURES = [
    ("pozicni_vazby_afixu", "009", "prefix__I_pos", ""),
    ("pozicni_vazby_afixu", "009", "suffix__I_pos", ""),
    ("sekcni_vazby_afixu", "009", "prefix__I_sec_g_pos", "blábol: sekce = pisatel"),
    ("sekcni_vazby_afixu", "009", "suffix__I_sec_g_pos", "blábol: sekce = pisatel"),
    ("registr_AB_uvnitr_pisare", "011", "a__suffix__hand__explained",
     "A/B nelze přímo; analog: podíl H(sufix|kmen) vysvětlený písařem/pisatelem"),
    ("sirka_paradigmatu", "009", "mean_suffixes_per_stem", ""),
    ("sirka_paradigmatu", "009", "H_suffix_given_stem", ""),
    ("sirka_paradigmatu", "009", "frac_stems_multi_suffix", ""),
    ("hapax", "009", "hapax_share", ""),
    ("nizke_poradi_tokenu", "009", "token_mi_excess_share", "Voynich pod textem"),
    ("nizke_poradi_tokenu", "009", "mi_stem_d1", "Voynich pod textem"),
    ("okrajova_MI", "009", "edge_mi_excess", "Voynich uvnitř rozpětí textů, nad generátory"),
    ("konec_radku", "010", "end_enrich_max", "blábol: vlastní řádky (gibberish_native)"),
    ("konec_radku", "009", "jsd_last_mid", "blábol: vlastní řádky jen v 010"),
    ("uzavreny_radek", "010", "bnd_edge_ratio_line", "Voynich ~0; blábol: vlastní řádky"),
    ("uzavreny_radek", "010", "bnd_affix_ratio_line", "Voynich ~0; blábol: vlastní řádky"),
    ("pokles_sdileni_se_vzdalenosti", "009", "loc_log_dist", "na hraně i u Voynichu"),
    ("paradigma_neni_pozicni", "011", "a__suffix__all-page__explained", "011: Voynich 7–8 %"),
    ("sufix_prefix_dalsiho_slova", "011", "a__nb_next_prefix__mi", "011: Voynich 0,038 bitu"),
    ("BPE_zlab", "-", "-", "nepočítáno (není v profilech 009–011)"),
]


def stat_table(runs: pd.DataFrame, metrics: list[str], family: str, origin: str) -> pd.DataFrame:
    rows = []
    for src, g in runs.groupby("source", sort=False):
        point = g[g["kind"] == "point"]
        for m in metrics:
            if m not in g.columns:
                continue
            if len(point):
                v = pd.to_numeric(point[m], errors="coerce").iloc[0]
                s = pd.to_numeric(g.loc[g["kind"] == "sub", m], errors="coerce").dropna()
                dev = s - s.mean()
                lo, hi = (v + SCALE * dev.quantile(Q[0]), v + SCALE * dev.quantile(Q[1])) if len(s) > 1 else (np.nan,) * 2
                sd = SCALE * s.std(ddof=1) if len(s) > 1 else np.nan
            else:
                s = pd.to_numeric(g[m], errors="coerce").dropna()
                v, lo, hi, sd = s.mean(), s.quantile(Q[0]), s.quantile(Q[1]), s.std(ddof=1)
            rows.append({"source": src, "family": family, "metric": m, "value": float(v), "lo": float(lo),
                         "hi": float(hi), "sd": float(sd), "n": len(s), "origin": origin})
    return pd.DataFrame(rows)


def load_new(kind: str) -> pd.DataFrame:
    base = OUT / "runs" / kind
    frames = [load_runs(d) for d in sorted(base.iterdir()) if (d / "runs.csv").exists()] if base.exists() else []
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def reference(kind: str) -> pd.DataFrame:
    r = pd.read_csv(RES / {"010": "010_line_as_item", "011": "011_suffix_variation"}[kind] / "runs.csv",
                    low_memory=False)
    r["kind"] = r["kind"].fillna("null")
    return r[r["source"].isin(REF_SOURCES)]


def build_profile() -> tuple[pd.DataFrame, pd.DataFrame]:
    new = {k: load_new(k) for k in ("p009", "p010", "p011")}
    allruns = pd.concat([d.assign(family=k[1:]) for k, d in new.items() if len(d)], ignore_index=True)
    parts = [stat_table(new["p009"], M009, "009", "013"), stat_table(new["p010"], M010, "010", "013"),
             stat_table(new["p011"], M011, "011", "013")]
    p9 = pd.read_csv(RES / "009_register_profile" / "profile.csv")
    p9 = p9[p9["source"].isin(REF_SOURCES)][["source", "metric", "value", "lo", "hi", "sd", "n"]]
    parts.append(p9.assign(family="009", origin="009"))
    parts.append(stat_table(reference("010"), M010, "010", "010"))
    parts.append(stat_table(reference("011"), M011, "011", "011"))
    return pd.concat(parts, ignore_index=True), allruns


def distances(prof: pd.DataFrame) -> pd.DataFrame:
    """RMS of standardized 009-profile differences (method 009; SD over all sources here, shuffles excluded)."""
    W = prof[prof["family"] == "009"].pivot_table(index="source", columns="metric", values="value")
    W = W[[m for m in M009 if m in W.columns]]
    scale = W.std(ddof=1).replace(0, np.nan)
    pairs = [("ctl:voynich_crop", s) for s in ("ctl:gibberish_pooled", "ctl:gibberish_native", "ctl:caesar_crop",
                                                "ctl:celsus_crop")]
    pairs += [("ZL_split", s) for s in ("ctl:gibberish_pooled", "ctl:naibbe_caesar", "ctl:naibbe_pliny",
                                        "ctl:naibbe_plaintext", *PROSE_FULL, "ts2020:ZL_split", "tmpl:ZL_split",
                                        "ctl:linnaeus_items")]
    pairs += [("GC_split", s) for s in ("ctl:naibbe_caesar", "ctl:naibbe_pliny", "ctl:gibberish_pooled")]
    pairs += [("ctl:naibbe_pliny", "ctl:pliny_prose"), ("ctl:naibbe_caesar", "ctl:caesar_prose"),
              ("ctl:naibbe_plaintext", "ctl:pliny_prose"), ("ctl:naibbe_pliny", "ctl:naibbe_plaintext"),
              ("ctl:naibbe_pliny", "ts2020:ZL_split"), ("ctl:naibbe_caesar", "ctl:naibbe_pliny")]
    pairs += [("ctl:gibberish_pooled", s) for s in ("ctl:caesar_crop", "ctl:celsus_crop", "ctl:naibbe_pliny")]
    rows = []
    for a, b in pairs:
        if a not in W.index or b not in W.index:
            continue
        z = (W.loc[b] - W.loc[a]) / scale
        rows.append({"a": a, "b": b, "metrics": "all", "distance": float(np.sqrt(np.nanmean(z ** 2))),
                     "n_metrics": int(z.notna().sum())})
        for g, ms in R.M9.GROUPS.items():
            zz = z[[m for m in ms if m in z.index]]
            rows.append({"a": a, "b": b, "metrics": f"only_{g}", "distance": float(np.sqrt(np.nanmean(zz ** 2))),
                         "n_metrics": int(zz.notna().sum())})
    return pd.DataFrame(rows)


def feature_lists(prof: pd.DataFrame) -> pd.DataFrame:
    P = prof.drop_duplicates(["source", "family", "metric"]).set_index(["source", "family", "metric"])

    def get(src, fam, m):
        try:
            return P.loc[(src, fam, m)]
        except KeyError:
            return None

    targets = {  # source -> (Voynich at that size, prose references at that size)
        "ctl:gibberish_pooled": ("ctl:voynich_crop", PROSE_CROP),
        "ctl:gibberish_native": ("ctl:voynich_crop", PROSE_CROP),
        "ctl:naibbe_caesar": ("ZL_split", PROSE_FULL),
        "ctl:naibbe_pliny": ("ZL_split", PROSE_FULL),
        "ctl:naibbe_plaintext": ("ZL_split", PROSE_FULL),
    }
    rows = []
    for feat, fam, m, note in FEATURES:
        for src, (vsrc, prose) in targets.items():
            base = {"feature": feat, "family": fam, "metric": m, "source": src, "voynich_ref": vsrc, "note": note}
            if fam == "-":
                rows.append({**base, "verdict": "nepočítáno"})
                continue
            if src == "ctl:gibberish_pooled" and fam == "010":
                rows.append({**base, "verdict": "viz gibberish_native (řádky v layoutu ZL jsou mechanické)"})
                continue
            s, v = get(src, fam, m), get(vsrc, fam, m)
            refs = [get(p, fam, m) for p in prose]
            refs = [r["value"] for r in refs if r is not None and np.isfinite(r["value"])]
            if s is None or v is None or not refs or not np.isfinite(s["value"]):
                rows.append({**base, "verdict": "nelze (chybí hodnota)"})
                continue
            r0 = float(np.mean(refs))
            den = v["value"] - r0
            f = (s["value"] - r0) / den if den else np.nan
            flo = (s["lo"] - r0) / den if den else np.nan
            fhi = (s["hi"] - r0) / den if den else np.nan
            lo_f, hi_f = (min(flo, fhi), max(flo, fhi)) if np.isfinite(flo) and np.isfinite(fhi) else (np.nan, np.nan)
            verdict = "ano" if f >= 0.5 else ("částečně" if f >= 0.2 else "ne")
            rows.append({**base, "value": s["value"], "lo": s["lo"], "hi": s["hi"], "voynich_value": v["value"],
                         "prose_mean": r0, "f": f, "f_lo": lo_f, "f_hi": hi_f, "verdict": verdict})
    return pd.DataFrame(rows)


def heatmap(prof: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    order = ["ZL_split", "ctl:voynich_crop", "ctl:gibberish_pooled", "ctl:gibberish_native", "ctl:naibbe_caesar",
             "ctl:naibbe_pliny", "ctl:naibbe_plaintext", "ctl:caesar_prose", "ctl:celsus_prose", "ctl:pliny_prose",
             "ctl:caesar_crop", "ctl:celsus_crop", "ctl:linnaeus_items", "ts2020:ZL_split", "tmpl:ZL_split"]
    W = prof[prof["family"] == "009"].pivot_table(index="source", columns="metric", values="value")
    order = [s for s in order if s in W.index]
    W = W.loc[order, [m for m in M009 if m in W.columns]]
    Z = (W - W.loc["ZL_split"]) / W.std(ddof=1).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(16, 6.5))
    im = ax.imshow(Z.clip(-3, 3).to_numpy(dtype=float), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([s.replace("ctl:", "") for s in order], fontsize=8)
    ax.set_xticks(range(W.shape[1]))
    ax.set_xticklabels(W.columns, rotation=75, ha="right", fontsize=7)
    x = 0
    for g, ms in R.M9.GROUPS.items():
        ax.axvline(x - 0.5, color="k", lw=0.8)
        ax.text(x + len(ms) / 2 - 0.5, -0.9, g, ha="center", fontsize=8)
        x += len(ms)
    fig.colorbar(im, ax=ax, fraction=0.02, label="(source - ZL split) / SD over sources, clipped at 3")
    ax.set_title("Register profile (009 metrics): gibberish and Naibbe vs Voynich ZL split (red = source higher)",
                 pad=18)
    fig.tight_layout()
    fig.savefig(OUT / "profile_heatmap.png", dpi=110)
    plt.close(fig)


def summarize() -> None:
    prof, runs = build_profile()
    runs.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    prof.to_csv(OUT / "profile.csv", index=False, float_format="%.6g")
    distances(prof).to_csv(OUT / "distances.csv", index=False, float_format="%.4g")
    feature_lists(prof).to_csv(OUT / "feature_lists.csv", index=False, float_format="%.5g")
    heatmap(prof)
    counts = runs.groupby(["family", "source", "kind"]).size().reset_index(name="n")
    write_meta(OUT, "013_gibberish_and_cipher_controls", {
        "sources": R.SOURCES, "reference_sources": REF_SOURCES, "writer_seed": R.WRITER_SEED,
        "feature_rule": "f = (source - prose mean) / (Voynich - prose mean); ano >= 0.5, castecne >= 0.2",
        "runs": counts.to_dict(orient="records")})
    print("written:", ", ".join(p.name for p in sorted(OUT.glob("*.*"))))


if __name__ == "__main__":
    summarize()
