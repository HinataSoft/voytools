"""Experiment 008: the canonical copying null (ts2020, 004 config 21) and the 005 template in 006/007,
the "positional stem" under three parsers, and the A/B suffix shift conditioned on position.

Steps 1 and 3 run the existing scripts of 006 and 007 with the sources added there (their outputs for
the earlier sources do not change)::

    python experiments/007_position_vs_content/run.py decomp --file ZL3b-n.txt --space-policy split \\
        --source selfcitation --out-root results/008_ts2020_recheck/007      # also: mi; --source template
    python experiments/006_topic_at_distance/run.py pages --combo ZL_split \\
        --sources selfcitation_21,template_005 --out-root results/008_ts2020_recheck/006

This script does the rest (from the repository root; every command is resumable)::

    python experiments/008_ts2020_recheck/run.py parser --key ZL      # step 2 (point + 100 subsamples)
    python experiments/008_ts2020_recheck/run.py ab --tag ZL_split    # step 4 (point, 200 subsamples, 200 perms)
    python experiments/008_ts2020_recheck/run.py profile --key ZL     # step 5: template on rows 1-3 of SYNTHESIS
    python experiments/008_ts2020_recheck/run.py check                # step 3 of the criteria: real values vs 005-007
    python experiments/008_ts2020_recheck/run.py summarize            # all tables of results/008_ts2020_recheck/

Seeds: ts2020 = 2100 + rep (the 004 phase-B samples of config 21), template = 100000 + rep with
``pool_seed = rep`` (the 005 ablation 'best' samples), real page subsamples as in 007 (seeds 0..99,
80 % of the pages without replacement) and 005 (seeds 1..200).  Decisions: ``notes.md``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import load_runs, run_seeds  # noqa: E402
from voytools._common import glyph_lists  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402
from voytools.morphology import apply_bpe, bpe_dependence_curve, decompose_table  # noqa: E402
from voytools.stats import edge_glyph_mi, jsd, positional_profile, vocabulary_summary  # noqa: E402

OUT = ROOT / "results" / "008_ts2020_recheck"
R004 = ROOT / "results" / "004_selfcitation_port"
R005 = ROOT / "results" / "005_structured_template_null"
R006 = ROOT / "results" / "006_topic_at_distance"
R007 = ROOT / "results" / "007_position_vs_content"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


E7 = _load("run", ROOT / "experiments" / "007_position_vs_content" / "run.py")  # 'run': summarize imports it
S7 = _load("summarize007", ROOT / "experiments" / "007_position_vs_content" / "summarize.py")
E6 = _load("exp006", ROOT / "experiments" / "006_topic_at_distance" / "run.py")
E5 = _load("exp005", ROOT / "experiments" / "005_structured_template_null" / "run.py")

FILES = E7.FILES
TAGS = [f"{k}_{p}" for k in FILES for p in E7.POLICIES]
N_NULL = 30
NULLS7 = {"autocopy": R007, "selfcitation": OUT / "007", "template": OUT / "007"}
NULLS6 = {"autocopy_86": R006, "perm_global": R006, "perm_section": R006, "perm_hand": R006,
          "selfcitation_21": OUT / "006", "template_005": OUT / "006"}
Z_LIMIT = 3.0

# --------------------------------------------------------------------------- #
# step 2: parser variants
# --------------------------------------------------------------------------- #

# (b) positional glyphs added to the learned affix lists.  Eva: gallows p f t k and the bench gallows
# cph cfh cth ckh (one glyph each under preset 'basic'; the spec's "cp, cf, ct, ck"), line-final m g am om.
# v101 equivalents from token alignment of ZL/GC loci with equal token counts (cmd ``check``,
# parser_affix_map.csv): p = g/j, f = f/u, t = k, k = h, cph = G/J, cfh = F/U, cth = K, ckh = H; m = p/P, g = *.
POS_AFFIXES = {
    "ZL": {"prefixes": ["p", "f", "t", "k", "cph", "cfh", "cth", "ckh"], "suffixes": ["m", "g", "am", "om"]},
    "GC": {"prefixes": ["g", "j", "f", "u", "k", "h", "G", "J", "F", "U", "K", "H"],
           "suffixes": ["p", "P", "*", "ap", "op"]},
}
EVA_POS = {"pre": ["p", "f", "t", "k", "cph", "cfh", "cth", "ckh"], "suf": ["m", "g", "am", "om"]}
VARIANTS = ("default", "positional", "bpe_mid", "bpe_mid3")  # bpe_mid3: only tokens with >= 3 BPE units
PARSER_TERMS = E7.TERMS + (("I_place", "place", None),)
N_SUB_PARSER = 100


def stem_codes(ctx, df: pd.DataFrame, variant: str) -> tuple[tuple, dict, np.ndarray | None]:
    """Stem class (top 200 + <other>) of ``variant``, descriptive numbers and the token mask (bpe_mid3)."""
    mask = None
    if variant.startswith("bpe_mid"):
        units = apply_bpe(glyph_lists(df, ctx.gpreset), ctx.merges)
        stem = np.asarray(["".join(u[1:-1]) for u in units], dtype=object)
        if variant == "bpe_mid3":
            mask = np.asarray([len(u) >= 3 for u in units])
            stem = stem[mask]
    else:
        pre, suf = list(ctx.prefixes), list(ctx.suffixes)
        if variant == "positional":
            add = POS_AFFIXES[ctx.tag.split("_")[0]]
            pre += [a for a in add["prefixes"] if a not in pre]
            suf += [a for a in add["suffixes"] if a not in suf]
        stem = decompose_table(df, prefixes=pre, suffixes=suf)["stem"].to_numpy(dtype=object)
    sc = pd.Series(stem).value_counts()
    top = set(sc.index[:E7.TOP_STEMS])
    info = {"n_stem_types": int(len(sc)), "share_empty_stem": float((stem == "").mean()),
            "share_top200": float(np.isin(stem, list(top)).mean()),
            "n_tokens": int(len(stem))}
    return E7.codes([s if s in top else "<other>" for s in stem]), info, mask


def place_codes(df: pd.DataFrame) -> tuple[np.ndarray, int]:
    f, last = df["is_first"].to_numpy(dtype=bool), df["is_last"].to_numpy(dtype=bool)
    return E7.codes(np.where(f & last, "single", np.where(f, "first", np.where(last, "last", "mid"))))


def parser_terms(ys: dict, X: dict, k_bias: int, seed: int) -> dict:
    """``E7.decomposition`` for several stem variants (same estimator, shuffle-corrected plug-in MI).

    ``ys``: variant -> (stem codes, token mask or None); the mask restricts X to the variant's tokens."""
    rng = np.random.default_rng(seed)
    out = {}
    for v, (y, mask) in ys.items():
        Xv = X if mask is None else {k: x[mask] for k, x in X.items()}
        zero = np.zeros(len(y), dtype=np.int64)
        sh = {z: E7.Shuffler(zero if z is None else Xv[z]) for z in {t[2] for t in PARSER_TERMS}}
        out[f"{v}__H"] = E7.H(y)
        for name, xn, zn in PARSER_TERMS:
            z = None if zn is None else Xv[zn]
            raw = E7.cmi(y, Xv[xn], z)
            bias = float(np.mean([E7.cmi(y[sh[zn].src(rng)], Xv[xn], z) for _ in range(k_bias)]))
            out[f"{v}__{name}__corr"] = raw - bias
    return out


def cmd_parser(args) -> None:
    ctx = E7.Context(FILES[args.key], "split")
    enc = ctx.encode(ctx.real)
    ys, info, masks = {}, {}, {}
    for v in VARIANTS:
        ys[v], info[v], masks[v] = stem_codes(ctx, ctx.real, v)
    X = {k: v[0] for k, v in enc["X"].items()}
    X["place"] = place_codes(ctx.real)[0]
    base = OUT / "parser" / ctx.tag
    base.mkdir(parents=True, exist_ok=True)
    if not (base / "point.csv").exists():
        pt = parser_terms({v: (c[0], masks[v]) for v, c in ys.items()}, X, E7.K_BIAS["point"], seed=0)
        pd.DataFrame([{"tag": ctx.tag, "source": "real", **pt,
                       **{f"{v}__{k}": x for v, d in info.items() for k, x in d.items()}}]).to_csv(
            base / "point.csv", index=False)

    def one(seed: int) -> dict:
        m = E7.page_subsample(enc, seed)
        Xs = {k: v[m] for k, v in X.items()}
        sub = {}
        for v, c in ys.items():  # stem codes of a masked variant are indexed by its own tokens
            mv = masks[v]
            sub[v] = (c[0][m], None) if mv is None else (c[0][m[mv]], mv[m])
        return {"tag": ctx.tag, "source": "real_sub", **parser_terms(sub, Xs, E7.K_BIAS["sub"], seed)}

    run_seeds("008_parser", list(range(N_SUB_PARSER)), one, base / "real_sub",
              params={"tag": ctx.tag, "variants": VARIANTS, "positional_affixes": POS_AFFIXES[args.key],
                      "prefixes": ctx.prefixes, "suffixes": ctx.suffixes, "bpe_merges": E7.BPE_MERGES,
                      "k_bias": E7.K_BIAS, "top_stems": E7.TOP_STEMS}, progress=E7.progress)


# --------------------------------------------------------------------------- #
# step 4: A/B suffix shift conditioned on position
# --------------------------------------------------------------------------- #

AB = (("AB_all", "lang", "A", "B"), ("H3_A_vs_B", "lh", "A3", "B3"))
MIN_CLASS = 20      # a position class enters the conditional JSD only if both groups have >= 20 tokens there
MIN_STD = 5         # direct standardisation: classes with >= 5 tokens in each group
N_AB_SUB, N_AB_PERM = 200, 200


def ab_values(vals: np.ndarray, key: np.ndarray, a: str, b: str, pc: np.ndarray, mid: np.ndarray,
              place: np.ndarray, support: list[str]) -> dict:
    ma, mb = key == a, key == b
    out = {"all": E5._jsd_of(vals, ma, mb, support), "mid": E5._jsd_of(vals, ma & mid, mb & mid, support)}
    for pl in ("first", "last"):
        out[pl] = E5._jsd_of(vals, ma & (place == pl), mb & (place == pl), support)
    for name, sel in (("cond", np.ones(len(vals), dtype=bool)), ("mid_cond", mid)):
        num = den = 0.0
        tot = int(((ma | mb) & sel).sum())
        for c in np.unique(pc[sel]):
            cm = sel & (pc == c)
            na, nb = int((ma & cm).sum()), int((mb & cm).sum())
            if na >= MIN_CLASS and nb >= MIN_CLASS:
                num += (na + nb) * E5._jsd_of(vals, ma & cm, mb & cm, support)
                den += na + nb
        out[name] = num / den if den else float("nan")
        out[f"{name}_coverage"] = den / tot if tot else float("nan")
    # direct standardisation: both suffix distributions reweighted to the pooled position-class mix of
    # the classes present in both groups, then one JSD (no small cells averaged separately)
    sidx = {v: i for i, v in enumerate(support)}
    vi = np.asarray([sidx[v] for v in vals])
    for name, sel in (("std", np.ones(len(vals), dtype=bool)), ("mid_std", mid)):
        pa, pb = np.zeros(len(support)), np.zeros(len(support))
        den = 0.0
        tot = int(((ma | mb) & sel).sum())
        for c in np.unique(pc[sel]):
            cm = sel & (pc == c)
            ia, ib = vi[ma & cm], vi[mb & cm]
            if len(ia) >= MIN_STD and len(ib) >= MIN_STD:
                w = len(ia) + len(ib)
                pa += w * np.bincount(ia, minlength=len(support)) / len(ia)
                pb += w * np.bincount(ib, minlength=len(support)) / len(ib)
                den += w
        out[name] = jsd(pa, pb) if den else float("nan")
        out[f"{name}_coverage"] = den / tot if tot else float("nan")
    return out


def ab_frame(dfd: pd.DataFrame) -> dict:
    pc = E7.position_classes(dfd)
    f, last = dfd["is_first"].to_numpy(dtype=bool), dfd["is_last"].to_numpy(dtype=bool)
    place = np.where(f & last, "single", np.where(f, "first", np.where(last, "last", "mid")))
    return {"vals": dfd["suffix"].astype(str).to_numpy(), "lang": dfd["language"].fillna("?").astype(str).to_numpy(),
            "lh": E5.lang_hand(dfd), "hand": dfd["hand"].fillna("?").astype(str).to_numpy(),
            "page": dfd["folio"].astype(str).to_numpy(), "pc": pc, "mid": ~(f | last), "place": place}


def ab_all(fr: dict, support: list[str], lang: np.ndarray | None = None) -> dict:
    lang = fr["lang"] if lang is None else lang
    lh = np.char.add(lang.astype(str), fr["hand"].astype(str))
    out = {}
    for name, kind, a, b in AB:
        for k, v in ab_values(fr["vals"], lang if kind == "lang" else lh, a, b, fr["pc"], fr["mid"], fr["place"],
                              support).items():
            out[f"{name}__{k}"] = v
    return out


def ab_perm(fr: dict, support: list[str], within_hand: bool, rng: np.random.Generator) -> dict:
    """Language labels permuted among pages (inside each scribe, or among all pages)."""
    pages, idx = np.unique(fr["page"], return_inverse=True)
    first = np.unique(idx, return_index=True)[1]
    plang, phand = fr["lang"][first].copy(), fr["hand"][first]
    if within_hand:
        for h in np.unique(phand):
            m = np.flatnonzero(phand == h)
            plang[m] = plang[m][rng.permutation(len(m))]
    else:
        plang = plang[rng.permutation(len(plang))]
    return ab_all(fr, support, plang[idx])


def cmd_ab(args) -> None:
    for tag in args.tag.split(","):
        ctx = E5.context(tag)  # identical real text and affix lists as 005 step 4
        fr = ab_frame(ctx.real_d)
        support = sorted(set(fr["vals"].tolist()))
        base = OUT / "ab" / tag
        base.mkdir(parents=True, exist_ok=True)
        if not (base / "perm.csv").exists():
            rng = np.random.default_rng(8008)
            rows = []
            for r in range(N_AB_PERM):
                for wh in (True, False):
                    rows.append({"rep": r, "null": "within_hand" if wh else "all_pages",
                                 **ab_perm(fr, support, wh, rng)})
            pd.DataFrame(rows).to_csv(base / "perm.csv", index=False, float_format="%.6g")

        def one(seed: int) -> dict:
            if seed == 0:
                return {"tag": tag, "source": "real", **ab_all(fr, support)}
            d = decompose_table(E5.page_subsample(ctx.real, seed), prefixes=ctx.prefixes, suffixes=ctx.suffixes)
            return {"tag": tag, "source": "real_sub", **ab_all(ab_frame(d), support)}

        run_seeds("008_ab", list(range(N_AB_SUB + 1)), one, base / "real",
                  params={"tag": tag, "min_class": MIN_CLASS, "n_perm": N_AB_PERM, "note": "seed 0 = full text"},
                  progress=E7.progress)


# --------------------------------------------------------------------------- #
# step 5: template on the SYNTHESIS rows the 005 runs do not cover
# --------------------------------------------------------------------------- #


def lines_of(df: pd.DataFrame) -> list[list[str]]:
    return df.groupby("locus_id", sort=False)["token"].apply(list).tolist()


def profile_metrics(df: pd.DataFrame, seed: int) -> dict:
    """Rows 1-3 and 5 of SYNTHESIS section 1, with the calls of 004 ``stats_for``."""
    curve = bpe_dependence_curve(lines_of(df), merge_checkpoints=(0, 32, 64, 128)).set_index("merges")["gap"]
    vs = vocabulary_summary(df).iloc[0]
    return {"edge_mi_excess": float(edge_glyph_mi(df, shuffles=20, seed=seed).iloc[0]["excess_bits"]),
            "pos_jsd_first_mid": float(positional_profile(df, which="first").jsd.loc["line_first", "mid"]),
            "bpe_gap64": float(curve.get(64, np.nan)), "hapax_share": float(vs["hapax_share"]),
            "n_types": float(vs["n_types"])}


def cmd_profile(args) -> None:
    ctx = E7.Context(FILES[args.key], "split")
    base = OUT / "profile" / ctx.tag

    def one(seed: int) -> dict:
        if seed < 0:
            return {"tag": ctx.tag, "source": "real", **profile_metrics(ctx.real, 0)}
        return {"tag": ctx.tag, "source": "template", **profile_metrics(ctx.corpus("template", seed), seed)}

    run_seeds("008_profile", [-1] + E7.null_seeds("template"), one, base,
              params={"tag": ctx.tag, "template_005": E7.template_config(args.key)}, progress=E7.progress)


# --------------------------------------------------------------------------- #
# reproduction check
# --------------------------------------------------------------------------- #


def cmd_check(args) -> None:
    """Real values recomputed with the unchanged code paths vs the files of 005/006/007; affix map v101."""
    res = {}
    for tag in TAGS:
        key, pol = tag.split("_")
        ctx = E7.Context(FILES[key], pol)
        enc = ctx.encode(ctx.real)
        pt = E7.decomposition(enc, E7.K_BIAS["point"], seed=0)
        old = pd.read_csv(R007 / "decomp" / tag / "point.csv")
        res[f"007_decomp_{tag}"] = float(max(abs(pt[f"{r.Y}__{r.term}__corr"] - r.corr) for r in old.itertuples()))
        if pol == "split":
            ms = E7.mi_stats(enc, E7.K_SHUF["point"], 0)
            om = load_runs(R007 / "mi" / tag / "real_point").iloc[0]
            res[f"007_mi_{tag}"] = float(max(abs(v - float(om[k])) for k, v in ms.items() if k in om.index))
        c6 = E6.Context(tag)
        r6 = E6.make_runner(c6, "real", c6.similarity(c6.real_d), None)(0)
        o6 = load_runs(R006 / "parts" / tag / "real").iloc[0]
        res[f"006_pages_{tag}"] = float(max(abs(v - float(o6[k])) for k, v in r6.items()
                                            if k in o6.index and isinstance(v, (int, float)) and np.isfinite(v)))
    # v101 equivalents of the Eva positional glyphs (aligned loci with equal token counts)
    z = E7.Context(FILES["ZL"], "split").real
    g = E7.Context(FILES["GC"], "split").real
    zl = z.groupby("locus_id", sort=False)["token"].apply(list)
    gl = g.groupby("locus_id", sort=False)["token"].apply(list)
    from voytools.glyphs import collapse_eva
    rows = []
    for loc in zl.index.intersection(gl.index):
        a, b = zl[loc], gl[loc]
        if len(a) != len(b):
            continue
        for x, y in zip(a, b):
            gx = collapse_eva(x, "basic")
            if gx[0] in EVA_POS["pre"]:
                rows.append(("prefix", gx[0], y[:1]))
            if gx[-1] in ("m", "g"):
                rows.append(("suffix", gx[-1], y[-1:]))
            if "".join(gx[-2:]) in ("am", "om"):
                rows.append(("suffix", "".join(gx[-2:]), y[-2:]))
    m = pd.DataFrame(rows, columns=["side", "eva", "v101"]).value_counts().rename("n").reset_index()
    m["share"] = m["n"] / m.groupby(["side", "eva"])["n"].transform("sum")
    m = m.sort_values(["side", "eva", "n"], ascending=[True, True, False])
    OUT.mkdir(parents=True, exist_ok=True)
    m[m["share"] >= 0.02].round(4).to_csv(OUT / "parser_affix_map.csv", index=False)
    (OUT / "check_reproduction.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #


def zrow(real: float, null: pd.Series) -> dict:
    s = null.dropna()
    sd = float(s.std(ddof=1)) if len(s) > 1 else float("nan")
    mean = float(s.mean()) if len(s) else float("nan")
    z = (real - mean) / sd if sd and np.isfinite(sd) and sd > 0 else float("nan")
    return {"null_mean": mean, "null_lo": float(s.quantile(0.025)) if len(s) else np.nan,
            "null_hi": float(s.quantile(0.975)) if len(s) else np.nan, "null_sd": sd, "n_seeds": int(len(s)),
            "z": z, "pct_of_real": 100 * mean / real if real else float("nan")}


def decomp_table() -> pd.DataFrame:
    rows = []
    for tag in TAGS:
        point = pd.read_csv(R007 / "decomp" / tag / "point.csv")
        sub = load_runs(R007 / "decomp" / tag / "real_sub")
        nulls = {s: load_runs(p / "decomp" / tag / s) for s, p in NULLS7.items()}
        for r in point.itertuples():
            col = f"{r.Y}__{r.term}__corr"
            lo, hi = S7.real_interval(r.corr, sub[col])
            for s, df in nulls.items():
                base = {"dataset": tag, "Y": r.Y, "term": r.term, "H": r.H, "real": r.corr, "real_lo": lo,
                        "real_hi": hi, "null": s}
                if df.empty:
                    rows.append({**base, "n_seeds": 0, "note": "no null (ts2020 generates Eva only)"})
                    continue
                rows.append({**base, **zrow(r.corr, df[col]), "null_H": float(df[f"{r.Y}__H"].mean()), "note": ""})
    return pd.DataFrame(rows)


def mi_table() -> pd.DataFrame:
    rows = []
    for tag in TAGS:
        point = load_runs(R007 / "mi" / tag / "real_point")
        sub = load_runs(R007 / "mi" / tag / "real_sub")
        nulls = {s: load_runs(p / "mi" / tag / s) for s, p in NULLS7.items()}
        for u in E7.UNITS:
            for corr in ("pc", "line"):
                p = S7.corrected(point, u, corr)[0]
                sm = S7.corrected(sub, u, corr)
                ps, ss = S7.curve_scalars(p[None, :]).iloc[0], S7.curve_scalars(sm)
                stats = {f"mi_d{d}": (p[d - 1], pd.Series(sm[:, d - 1])) for d in range(1, E7.DMAX + 1)}
                stats["area_2_10"] = (ps["area_2_10"], ss["area_2_10"])
                for stat, (pv, sv) in stats.items():
                    lo, hi = S7.real_interval(float(pv), sv)
                    for s, df in nulls.items():
                        base = {"dataset": tag, "U": u, "correction": corr, "stat": stat, "real": float(pv),
                                "real_lo": lo, "real_hi": hi, "null": s}
                        if df.empty:
                            rows.append({**base, "n_seeds": 0, "note": "no null (ts2020 generates Eva only)"})
                            continue
                        m = S7.corrected(df, u, corr)
                        nv = S7.curve_scalars(m)["area_2_10"] if stat == "area_2_10" else \
                            pd.Series(m[:, int(stat[4:]) - 1])
                        rows.append({**base, **zrow(float(pv), nv), "note": ""})
    return pd.DataFrame(rows)


PAGE_TERMS = ["s_log_dist", "s_same_section", "s_same_language", "s_same_hand", "contrast", "contrast_adj",
              "h_contrast", "b_log_dist", "b_same_section"]


def pages_table() -> pd.DataFrame:
    rows = []
    for combo in E6.COMBOS:
        real = load_runs(R006 / "parts" / combo / "real").iloc[0]
        boot = load_runs(R006 / "parts" / combo / "real_boot")
        nulls = {s: load_runs(p / "parts" / combo / s) for s, p in NULLS6.items()}
        for m in E6.METRICS:
            for k in E6.KS:
                pre = f"{m}_k{k}_"
                curve_cols = [f"c_{sm}_{b}" for sm in ("same", "diff") for b in E6.BIN_NAMES]
                for t in PAGE_TERMS + curve_cols:
                    c = pre + t
                    lo, hi = boot[c].quantile(0.025), boot[c].quantile(0.975)
                    kind = "curve" if t.startswith("c_") else ("contrast" if "contrast" in t else "ols")
                    for s, df in nulls.items():
                        base = {"combo": combo, "metric": m, "K": k, "table": kind, "term": t, "real": float(real[c]),
                                "real_lo": float(lo), "real_hi": float(hi), "null": s}
                        if df.empty:
                            rows.append({**base, "n_seeds": 0, "note": "no null (ts2020 generates Eva only)"})
                            continue
                        v = df[c]
                        zr = zrow(float(real[c]), v)
                        zr["p_null_le_real"] = float((1 + (v <= real[c]).sum()) / (1 + v.notna().sum()))
                        zr["p_null_ge_real"] = float((1 + (v >= real[c]).sum()) / (1 + v.notna().sum()))
                        rows.append({**base, **zr, "note": ""})
    return pd.DataFrame(rows)


def parser_table() -> pd.DataFrame:
    rows = []
    for key in FILES:
        tag = f"{key}_split"
        point = pd.read_csv(OUT / "parser" / tag / "point.csv").iloc[0]
        sub = load_runs(OUT / "parser" / tag / "real_sub")
        for v in VARIANTS:
            for name, _, _ in PARSER_TERMS:
                col = f"{v}__{name}__corr"
                lo, hi = S7.real_interval(float(point[col]), sub[col])
                rows.append({"dataset": tag, "variant": v, "term": name, "corr": float(point[col]), "lo": lo, "hi": hi,
                             "H": float(point[f"{v}__H"]), "share": float(point[col]) / float(point[f"{v}__H"]),
                             "n_sub": int(sub[col].notna().sum()), "n_stem_types": int(point[f"{v}__n_stem_types"]),
                             "share_empty_stem": float(point[f"{v}__share_empty_stem"]),
                             "share_top200": float(point[f"{v}__share_top200"]),
                             "n_tokens": int(point[f"{v}__n_tokens"])})
            # paired differences (same subsample): position minus section, position minus section | hand
            for a, b in (("I_pos", "I_sec"), ("I_pos", "I_sec_g_hand"), ("I_pos_g_sec", "I_sec_g_pos"),
                         ("I_place", "I_sec")):
                d0 = float(point[f"{v}__{a}__corr"] - point[f"{v}__{b}__corr"])
                lo, hi = S7.real_interval(d0, sub[f"{v}__{a}__corr"] - sub[f"{v}__{b}__corr"])
                rows.append({"dataset": tag, "variant": v, "term": f"{a}-{b}", "corr": d0, "lo": lo, "hi": hi,
                             "n_sub": len(sub)})
    return pd.DataFrame(rows)


def ab_table() -> pd.DataFrame:
    rows = []
    for tag in TAGS:
        runs = load_runs(OUT / "ab" / tag / "real")
        real, sub = runs[runs["source"] == "real"].iloc[0], runs[runs["source"] == "real_sub"]
        perm = pd.read_csv(OUT / "ab" / tag / "perm.csv")
        for name, _, _, _ in AB:
            for stat in ("all", "cond", "std", "mid", "mid_cond", "mid_std", "first", "last"):
                c = f"{name}__{stat}"
                sd_real = 2 * float(sub[c].std(ddof=1))  # delete-d jackknife, as in 005
                row = {"dataset": tag, "contrast": name, "stat": stat, "jsd": float(real[c]),
                       "lo": float(real[c] - 1.96 * sd_real), "hi": float(real[c] + 1.96 * sd_real),
                       "coverage": float(real.get(f"{c}_coverage", np.nan)) if f"{c}_coverage" in real else np.nan}
                for nl, pm in perm.groupby("null"):
                    v = pm[c].dropna()
                    ex = row["jsd"] - float(v.mean())
                    hw = 1.96 * np.sqrt(sd_real ** 2 + float(v.std(ddof=1)) ** 2)
                    row.update({f"perm_{nl}_mean": float(v.mean()), f"perm_{nl}_q975": float(v.quantile(0.975)),
                                f"excess_{nl}": ex, f"excess_{nl}_lo": ex - hw, f"excess_{nl}_hi": ex + hw,
                                f"p_{nl}": float((1 + (v >= row["jsd"]).sum()) / (1 + len(v)))})
                rows.append(row)
    out = pd.DataFrame(rows)
    # ratio conditional / unconditional (raw JSD and excess over the within-scribe permutation)
    for (tag, name), g in out.groupby(["dataset", "contrast"]):
        g = g.set_index("stat")
        for stat, ref in (("cond", "all"), ("std", "all"), ("mid_cond", "mid"), ("mid_std", "mid"), ("mid", "all"),
                          ("mid_std", "all")):
            for col in ("jsd", "excess_within_hand", "excess_all_pages"):
                out.loc[(out["dataset"] == tag) & (out["contrast"] == name) & (out["stat"] == stat),
                        f"ratio_{col}_vs_{ref}"] = g.loc[stat, col] / g.loc[ref, col]
    return out


# ---- step 5: synthesis rows

SYN_ROWS = [  # (row, property, metric, where the real value and the null values come from)
    (1, "Okrajová MI (poslední→první glyf)", "edge_mi_excess", "profile"),
    (2, "Poziční JSD (první glyf řádku vs. střed)", "pos_jsd_first_mid", "profile"),
    (3, "BPE žlab @64", "bpe_gap64", "profile"),
    (4, "Sufixů na kmen (≥5 tok.)", "mean_suffixes_per_stem", "paradigm"),
    (5, "Hapax (podíl typů)", "hapax_share", "profile"),
    (6, "jsd_suffix_AB (frekvence sufixů A vs B)", "jsd_suffix_AB", "paradigm"),
    (7, "Sekční informace kmenů nezávislá na písaři (007)", "stem__I_sec_g_hand", "decomp"),
    (7, "Sekční informace kmenů nezávislá na písaři (007)", "stem__I_sec_g_pos", "decomp"),
    (7, "Sekční informace kmenů nezávislá na písaři (007)", "stem__I_sec", "decomp"),
    (8, "Strmý pokles sdílení kmenů se vzdáleností (006; cos, K=30)", "cos_k30_s_log_dist", "pages"),
    (8, "Strmý pokles sdílení kmenů se vzdáleností (006; cos, K=30)", "cos_k30_contrast_adj", "pages"),
    (9, "Žádná reziduální MI za d = 2 (007; kmen)", "stem_area_2_10_pc", "mi"),
    (9, "Žádná reziduální MI za d = 2 (007; kmen)", "stem_area_2_10_line", "mi"),
    (9, "Žádná reziduální MI za d = 2 (007; kmen)", "stem_mi_d3_pc", "mi"),
    (9, "Žádná reziduální MI za d = 2 (007; kmen)", "stem_mi_d1_pc", "mi"),
]


def synthesis_table(dec: pd.DataFrame, mi: pd.DataFrame, pages: pd.DataFrame) -> pd.DataFrame:
    rows = []
    r004 = {pol: load_runs(R004 / "phaseB" / f"ZL_{pol}") for pol in ("split", "join")}
    real004 = {pol: pd.read_csv(R004 / "real" / f"ZL_{pol}" / "point.csv").iloc[0] for pol in ("split", "join")}
    abl = {k: load_dir_005(R005 / "ablation" / f"{k}_split" / "best") for k in FILES}
    real005 = {k: E5.real_point(f"{k}_split") for k in FILES}
    prof = {k: load_runs(OUT / "profile" / f"{k}_split") for k in FILES}
    for row, prop, metric, kind in SYN_ROWS:
        def add(dataset, null, real, vals, source, note=""):
            if vals is None:
                rows.append({"row": row, "property": prop, "metric": metric, "dataset": dataset, "null": null,
                             "real": real, "n_seeds": 0, "source_file": source, "note": note})
                return
            zr = zrow(float(real), pd.Series(vals, dtype=float))
            rows.append({"row": row, "property": prop, "metric": metric, "dataset": dataset, "null": null,
                         "real": float(real), **zr, "reproduces": bool(abs(zr["z"]) < Z_LIMIT),
                         "source_file": source, "note": note})

        if kind in ("profile", "paradigm"):
            for pol in ("split", "join"):
                cid = r004[pol][r004[pol]["config_id"] == 21]
                add(f"ZL_{pol}", "ts2020", real004[pol][metric], cid[metric],
                    f"results/004_selfcitation_port/phaseB/ZL_{pol}/runs.csv (config 21, 50 seeds)")
            add("GC_split", "ts2020", np.nan, None, "", "no null: ts2020 generates Eva only")
            for k in FILES:
                if kind == "paradigm" or metric == "hapax_share":
                    add(f"{k}_split", "tmpl", real005[k][metric], abl[k][metric],
                        f"results/005_structured_template_null/ablation/{k}_split/best (50 seeds)")
                else:
                    p = prof[k]
                    add(f"{k}_split", "tmpl", p[p["source"] == "real"].iloc[0][metric],
                        p[p["source"] == "template"][metric],
                        f"results/008_ts2020_recheck/profile/{k}_split (30 seeds)")
        elif kind == "decomp":
            y, term = metric.split("__")
            for _, r in dec[(dec["Y"] == y) & (dec["term"] == term) & dec["null"].isin(["selfcitation", "template"])].iterrows():
                null = "ts2020" if r["null"] == "selfcitation" else "tmpl"
                df = load_runs(OUT / "007" / "decomp" / r["dataset"] / r["null"])
                add(r["dataset"], null, r["real"], df[f"{metric}__corr"] if not df.empty else None,
                    f"results/008_ts2020_recheck/007/decomp/{r['dataset']}/{r['null']}",
                    "" if not df.empty else "no null: ts2020 generates Eva only")
        elif kind == "pages":
            for _, r in pages[(pages["metric"] == "cos") & (pages["K"] == 30) & (pages["term"] == metric[8:])
                              & pages["null"].isin(["selfcitation_21", "template_005"])].iterrows():
                null = "ts2020" if r["null"] == "selfcitation_21" else "tmpl"
                df = load_runs(OUT / "006" / "parts" / r["combo"] / r["null"])
                add(r["combo"], null, r["real"], df[metric] if not df.empty else None,
                    f"results/008_ts2020_recheck/006/parts/{r['combo']}/{r['null']}",
                    "" if not df.empty else "no null: ts2020 generates Eva only")
        elif kind == "mi":
            u, rest = metric.split("_", 1)
            stat, corr = rest.rsplit("_", 1)
            for tag in TAGS:
                point = load_runs(R007 / "mi" / tag / "real_point")
                real = float(S7.curve_scalars(S7.corrected(point, u, corr)).iloc[0][stat])
                for s, null in (("selfcitation", "ts2020"), ("template", "tmpl")):
                    df = load_runs(OUT / "007" / "mi" / tag / s)
                    vals = S7.curve_scalars(S7.corrected(df, u, corr))[stat] if not df.empty else None
                    add(tag, null, real, vals, f"results/008_ts2020_recheck/007/mi/{tag}/{s}",
                        "" if vals is not None else "no null: ts2020 generates Eva only")
    out = pd.DataFrame(rows)
    # verdict per row and null: reproduces in every set that has the null (ts2020: ZL split and join)
    key = ["row", "metric", "null"]
    v = out[out["n_seeds"] > 0].groupby(key)["reproduces"].agg(["all", "any", "size"])
    out["verdict"] = [("ano" if v.loc[k, "all"] else ("částečně" if v.loc[k, "any"] else "ne"))
                      if k in v.index else "" for k in zip(out["row"], out["metric"], out["null"])]
    return out


def load_dir_005(path: Path) -> pd.DataFrame:
    return pd.concat([load_runs(p) for p in sorted(path.glob("part*"))], ignore_index=True)


def cmd_summarize(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dec = decomp_table()
    dec.round(6).to_csv(OUT / "decomp_ts2020.csv", index=False)
    mi = mi_table()
    mi.round(6).to_csv(OUT / "mi_distance_ts2020.csv", index=False)
    pages = pages_table()
    pages.round(6).to_csv(OUT / "pages_ts2020.csv", index=False)
    parser_table().round(6).to_csv(OUT / "parser_variants.csv", index=False)
    ab_table().round(6).to_csv(OUT / "ab_by_position.csv", index=False)
    syn = synthesis_table(dec, mi, pages)
    syn.round(6).to_csv(OUT / "synthesis_update.csv", index=False)
    n = {"decomp": len(dec), "mi": len(mi), "pages": len(pages), "synthesis": len(syn)}
    write_meta(OUT, "008_ts2020_recheck", {"z_limit": Z_LIMIT, "n_null": N_NULL, "nulls_007": list(NULLS7),
                                           "nulls_006": list(NULLS6)}, n_rows=n)
    print(syn[["row", "metric", "dataset", "null", "real", "null_mean", "z", "reproduces", "verdict"]].to_string())


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("parser")
    p.add_argument("--key", default="ZL", choices=list(FILES))
    p = sub.add_parser("ab")
    p.add_argument("--tag", default=",".join(TAGS))
    p = sub.add_parser("profile")
    p.add_argument("--key", default="ZL", choices=list(FILES))
    sub.add_parser("check")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    {"parser": cmd_parser, "ab": cmd_ab, "profile": cmd_profile, "check": cmd_check,
     "summarize": cmd_summarize}[args.cmd](args)


if __name__ == "__main__":
    main()
