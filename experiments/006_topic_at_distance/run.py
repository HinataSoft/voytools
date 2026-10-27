"""Experiment 006: topic at a distance - stem sharing between pages, labels versus running text.

Usage (from the repository root; see ``spec.md``)::

    # page pairs: real, page bootstrap and all null sources of one combination (resumable)
    python experiments/006_topic_at_distance/run.py pages --combo ZL_split
    python experiments/006_topic_at_distance/run.py pages --combo GC_join --sources autocopy_86 --seeds 0:25
    # labels (hits of label stems in the text, permutation of labels within sections)
    python experiments/006_topic_at_distance/run.py labels --combo ZL_split
    # ok-/ot- tokens in running text (position, neighbours, length) vs qo-, ch-/sh-, d-
    python experiments/006_topic_at_distance/run.py okot --combo ZL_split
    # merge, curve.csv, ols.csv, contrast.csv, curve.png, labels.png
    python experiments/006_topic_at_distance/run.py summarize
    # outside the spec: pairs.csv, v101 initial classes, sensitivity to the T / C / A pages
    python experiments/006_topic_at_distance/run.py pairs
    python experiments/006_topic_at_distance/run.py initials
    python experiments/006_topic_at_distance/run.py sections

Combinations: {ZL3b, GC2a} x {split, join}.  Sources of the page part (``seed`` = repetition):

* ``real`` - the real text, seed 0;
* ``real_boot`` - bootstrap over pages (pairs of resampled pages), seeds 0..199;
* ``perm_global`` / ``perm_section`` / ``perm_hand`` - page contents permuted among page positions
  globally / inside the section (``illustration``) / inside the scribe (``hand``); metadata
  (position, section, language, hand) stay with the position; seeds 0..99 each;
* ``autocopy_86`` - auto-copy generator, config 86 of 003 (``edit_model='trigram'``), seeds 0..49.

``voytools.selfcitation`` (experiment 004) is not merged, so the self-citation null is skipped.

Stems: ``decompose_table`` with affixes learned on the real P text of the combination; every sample is
decomposed with the same lists.  Function stems: the K most frequent stems of the sample itself
(K = 0, 30, 100; for the permutation nulls identical to the real list).  Similarities: cosine of
``log(1+f)`` stem vectors and Jaccard of stems with >= 2 occurrences, pages with >= 40 clean P tokens.
Distance = |rank of page i - rank of page j| in file order over all pages of the file.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_runs, load_tables, run_seeds, seed_range  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.nulls import autocopy_like, shuffle_within_lines  # noqa: E402
from voytools.stats import (  # noqa: E402
    DISTANCE_BIN_EDGES,
    distance_bins,
    jsd,
    ols_fit,
    page_similarity,
    page_term_counts,
)
from voytools.wordgrammar import learn_word_grammar  # noqa: E402

DATA = ROOT / "data"
OUT = ROOT / "results" / "006_topic_at_distance"
BEST003 = ROOT / "results" / "003_autocopy_strict" / "best_configs.csv"

FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
POLICIES = ("split", "join")
COMBOS = [f"{k}_{p}" for k in FILES for p in POLICIES]
MIN_TOKENS = 40
KS = (0, 30, 100)
MAIN_K = 30
METRICS = ("cos", "jac")
FAR, NEAR = 16, 3  # contrast: same section dist >= FAR minus other section dist <= NEAR
SOURCES = {"real": range(1), "real_boot": range(200), "perm_global": range(100),
           "perm_section": range(100), "perm_hand": range(100), "autocopy_86": range(50)}
GEN_PARAMS = ("edit_model", "new_token_source", "p_new", "window", "p_copy_exact", "p_two_edits", "recency")
LABEL_TYPES = ("Lz", "Lf", "Lc", "Ls", "Ln")
LABEL_PERMS = 200
LABEL_BOOT = 200
OKOT_SHUFFLES = 20
OKOT_BOOT = 200

# Token-initial classes.  Eva prefixes; v101 (GC) prefixes derived by aligning ZL and GC loci with the
# same number of tokens (cmd ``initials`` -> initial_map.csv): Eva ok- = v101 oh-/Ah-, ot- = ok-/Ak-,
# qo- = 4o-, ch- = 1-, sh- = 2-/3-/5-, d- = 8-/7-.
CLASS_PREFIXES = {
    "ZL": {"ok": ("ok",), "ot": ("ot",), "qo": ("qo",), "chsh": ("ch", "sh"), "d": ("d",)},
    "GC": {"ok": ("oh", "Ah"), "ot": ("ok", "Ak"), "qo": ("4o",), "chsh": ("1", "2", "3", "5"), "d": ("8", "7")},
}
BASE_CLASSES = ("ok", "ot", "qo", "chsh", "d", "other")
NEIGH = ("okot", "qo", "chsh", "d", "other", "edge")


def parse_seeds(spec: str) -> list[int]:
    """``'0:20'`` -> seeds 0..19; ``'0:100/4:1'`` -> chunk 1 of 4 of 0..99; ``'3,5,8'`` -> those seeds."""
    if "/" in spec:
        rng, part = spec.split("/")
        chunks, k = part.split(":")
        a, b = rng.split(":")
        return seed_range(int(a), int(b), int(chunks))[int(k)]
    if ":" in spec:
        a, b = spec.split(":")
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(",") if s.strip()]


def autocopy_params() -> dict:
    """Generator keyword arguments of config 86 (row 1 of 003's best_configs.csv)."""
    best = pd.read_csv(BEST003)
    r = best[best["config_id"].astype(str) == "86"].iloc[0]
    p = {k: r[k] for k in GEN_PARAMS}
    for k in ("p_new", "p_copy_exact", "p_two_edits", "recency"):
        p[k] = float(p[k])
    p["window"] = int(p["window"])
    p["edit_type_probs"] = (0.6, 0.2, 0.2)
    return p


def classify(tokens, key: str) -> np.ndarray:
    """Token-initial class (``BASE_CLASSES``) of every token, alphabet ``key`` ('ZL' or 'GC')."""
    pref = CLASS_PREFIXES[key]
    out = np.full(len(tokens), "other", dtype=object)
    for i, t in enumerate(tokens):
        for c, ps in pref.items():
            if t.startswith(ps):
                out[i] = c
                break
    return out


def col(name: str) -> str:
    return name.replace("-", "_").replace("+", "p")


BIN_NAMES = [col(b) for b in distance_bins(np.asarray(DISTANCE_BIN_EDGES))]


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #


class Context:
    """Real text of one combination, page metadata, pair arrays and affix lists."""

    def __init__(self, combo: str, with_grammar: bool = False, drop_sections: tuple = ()):
        self.combo = combo
        self.key, self.policy = combo.split("_")
        self.file = FILES[self.key]
        tables = load_tables(DATA / self.file, space_policy=self.policy)
        self.all_tokens = tables.tokens
        self.real = filter_tokens(tables.tokens, generic_type="P", clean=True)
        d = decompose_table(self.real)
        self.real_d = d
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]
        order = pd.unique(tables.tokens["folio"].to_numpy())
        rank = {p: i for i, p in enumerate(order)}
        size = self.real.groupby("folio", sort=False).size()
        sec_of = self.real.groupby("folio", sort=False)["illustration"].first().astype(str)
        self.pages = [p for p in pd.unique(self.real["folio"].to_numpy())
                      if size[p] >= MIN_TOKENS and sec_of[p] not in drop_sections]
        meta = self.real.groupby("folio", sort=False)[["illustration", "language", "hand"]].first().loc[self.pages]
        self.meta = meta.assign(rank=[rank[p] for p in self.pages])
        n = len(self.pages)
        self.n = n
        self.size = size.loc[self.pages].to_numpy(dtype=float)
        self.rank = self.meta["rank"].to_numpy()
        self.sec = self.meta["illustration"].astype(str).to_numpy()
        self.lang = self.meta["language"].astype(str).to_numpy()
        self.hand = self.meta["hand"].astype(str).to_numpy()
        self.iu, self.ju = np.triu_indices(n, 1)
        self.grammar = None
        if with_grammar:
            self.grammar = learn_word_grammar(self.real["token"].tolist(), self.real.attrs.get("glyph_preset"))

    def columns(self) -> dict:
        return {"combo": self.combo, "file": self.file, "space_policy": self.policy}

    def decompose(self, df: pd.DataFrame) -> pd.DataFrame:
        return decompose_table(df, prefixes=self.prefixes, suffixes=self.suffixes)

    def top_stems(self, dfd: pd.DataFrame, k: int) -> list[str]:
        vc = dfd.loc[dfd["stem"] != "", "stem"].value_counts()
        return vc.index[:k].tolist()

    def similarity(self, dfd: pd.DataFrame) -> dict:
        """``{(metric, K): n x n similarity matrix}`` over ``self.pages``."""
        out = {}
        top = self.top_stems(dfd, max(KS))
        for k in KS:
            counts = page_term_counts(dfd, "stem", exclude=top[:k]).reindex(self.pages, fill_value=0)
            out[("cos", k)] = page_similarity(counts.to_numpy(), "cosine")
            out[("jac", k)] = page_similarity(counts.to_numpy(), "jaccard", min_count=2)
        return out


def pair_design(ctx: Context, a: np.ndarray, b: np.ndarray) -> dict:
    """Covariates of page pairs given as position indices ``a``, ``b`` (a != b)."""
    dist = np.abs(ctx.rank[a] - ctx.rank[b])
    return {
        "dist": dist,
        "log_dist": np.log(dist),
        "same_section": (ctx.sec[a] == ctx.sec[b]).astype(float),
        "same_language": (ctx.lang[a] == ctx.lang[b]).astype(float),
        "same_hand": (ctx.hand[a] == ctx.hand[b]).astype(float),
        "bin": np.asarray([col(x) for x in distance_bins(dist)], dtype=object),
    }


def pair_metrics(share: np.ndarray, des: dict, size_a: np.ndarray, size_b: np.ndarray) -> dict:
    """OLS coefficients, contrast and distance curve for one vector of pair similarities.

    ``size_a``/``size_b`` are the token counts of the two pages (they travel with the content); the
    ``s_*`` model adds ``log n_i + log n_j`` and ``|log n_i - log n_j|`` to the spec's model and
    ``contrast_adj`` is the contrast of the residuals of ``share ~ size terms`` (outside the spec:
    cosine and Jaccard grow with page size, and page size differs strongly between sections).
    """
    out = {"mean": float(share.mean())}
    X = np.column_stack([des["log_dist"], des["same_section"], des["same_language"], des["same_hand"]])
    r = ols_fit(share, X, ["log_dist", "same_section", "same_language", "same_hand"])
    for t in ("intercept", "log_dist", "same_section", "same_language", "same_hand", "r2"):
        out[f"b_{t}"] = float(r[t])
    ss, dist = des["same_section"] == 1, des["dist"]
    far, near = ss & (dist >= FAR), ~ss & (dist <= NEAR)
    out["contrast"] = float(share[far].mean() - share[near].mean()) if far.any() and near.any() else np.nan
    # within-scribe control: only pairs written by the same hand
    h = des["same_hand"] == 1
    rh = ols_fit(share[h], np.column_stack([des["log_dist"][h], des["same_section"][h]]), ["log_dist", "same_section"])
    out["h_log_dist"], out["h_same_section"] = float(rh["log_dist"]), float(rh["same_section"])
    fh, nh = far & h, near & h
    out["h_contrast"] = float(share[fh].mean() - share[nh].mean()) if fh.any() and nh.any() else np.nan
    la, lb = np.log(size_a), np.log(size_b)
    size_terms = np.column_stack([la + lb, np.abs(la - lb)])
    rs = ols_fit(share, np.column_stack([X, size_terms]),
                 ["log_dist", "same_section", "same_language", "same_hand", "size_sum", "size_diff"])
    for t in ("log_dist", "same_section", "same_language", "same_hand", "size_sum", "size_diff", "r2"):
        out[f"s_{t}"] = float(rs[t])
    r0 = ols_fit(share, size_terms, ["size_sum", "size_diff"])
    resid = share - (r0["intercept"] + size_terms @ np.asarray([r0["size_sum"], r0["size_diff"]]))
    out["contrast_adj"] = float(resid[far].mean() - resid[near].mean()) if far.any() and near.any() else np.nan
    out["h_contrast_adj"] = float(resid[fh].mean() - resid[nh].mean()) if fh.any() and nh.any() else np.nan
    for same, mask in (("same", ss), ("diff", ~ss)):
        for b in BIN_NAMES:
            m = mask & (des["bin"] == b)
            out[f"c_{same}_{b}"] = float(share[m].mean()) if m.any() else np.nan
    return out


def all_metrics(ctx: Context, sims: dict, a: np.ndarray, b: np.ndarray, sa: np.ndarray, sb: np.ndarray) -> dict:
    """Metrics for every (metric, K); covariates from positions ``a, b``, similarities from contents ``sa, sb``."""
    des = pair_design(ctx, a, b)
    out = {}
    for (m, k), S in sims.items():
        for name, v in pair_metrics(S[sa, sb], des, ctx.size[sa], ctx.size[sb]).items():
            out[f"{m}_k{k}_{name}"] = round(v, 6) if np.isfinite(v) else v
    s = sims[("cos", 0)][sa, sb]
    out["max_cos_k0"] = round(float(s.max()), 6)
    out["n_cos_k0_gt095"] = int((s > 0.95).sum())
    out["n_pairs"] = len(a)
    return out


def permutation(ctx: Context, groups: np.ndarray | None, seed: int) -> np.ndarray:
    """Content permutation: position p gets the content of page ``perm[p]`` (inside ``groups``)."""
    rng = np.random.default_rng(seed)
    n = ctx.n
    if groups is None:
        return rng.permutation(n)
    perm = np.arange(n)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        perm[idx] = idx[rng.permutation(len(idx))]
    return perm


def make_runner(ctx: Context, source: str, real_sims: dict, ac_params: dict | None):
    iu, ju = ctx.iu, ctx.ju

    def one(seed: int) -> dict:
        base = {**ctx.columns(), "source": source}
        if source == "real":
            return {**base, **all_metrics(ctx, real_sims, iu, ju, iu, ju)}
        if source == "real_boot":
            rng = np.random.default_rng(seed)
            idx = rng.integers(ctx.n, size=ctx.n)
            i2, j2 = np.triu_indices(ctx.n, 1)
            a, b = idx[i2], idx[j2]
            keep = a != b
            a, b = a[keep], b[keep]
            return {**base, **all_metrics(ctx, real_sims, a, b, a, b)}
        if source.startswith("perm_"):
            groups = {"perm_global": None, "perm_section": ctx.sec, "perm_hand": ctx.hand}[source]
            p = permutation(ctx, groups, seed)
            return {**base, **all_metrics(ctx, real_sims, iu, ju, p[iu], p[ju])}
        if source == "autocopy_86":
            sample = autocopy_like(ctx.real, allowed_from=ctx.real["token"].tolist(), word_grammar=ctx.grammar,
                                   seed=10_000 + seed, **ac_params)
            sims = ctx.similarity(ctx.decompose(sample))
            return {**base, **all_metrics(ctx, sims, iu, ju, iu, ju)}
        raise ValueError(source)

    return one


def progress(i: int, n: int, s: int) -> None:
    if i % 25 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_pages(args) -> None:
    sources = list(SOURCES) if args.sources == "all" else args.sources.split(",")
    ctx = Context(args.combo, with_grammar="autocopy_86" in sources)
    real_sims = ctx.similarity(ctx.real_d)
    ac = autocopy_params()
    parts = OUT / "parts" / args.combo
    # list of excluded function stems (real text) and the page table, written once per combination
    top = ctx.top_stems(ctx.real_d, max(KS))
    counts = ctx.real_d.loc[ctx.real_d["stem"] != "", "stem"].value_counts()
    parts.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"rank": range(1, len(top) + 1), "stem": top, "count": [int(counts[s]) for s in top]}).to_csv(
        parts / "function_stems.csv", index=False)
    ctx.meta.reset_index().to_csv(parts / "pages.csv", index=False)
    for source in sources:
        seeds = parse_seeds(args.seeds) if args.seeds else list(SOURCES[source])
        print(f"{args.combo} {source}: {len(seeds)} seeds", flush=True)
        params = {"combo": args.combo, "source": source, "min_tokens": MIN_TOKENS, "Ks": KS,
                  "prefixes": ctx.prefixes, "suffixes": ctx.suffixes, "n_pages": ctx.n}
        if source == "autocopy_86":
            params["autocopy"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in ac.items()}
            params["seed_offset"] = 10_000
        run_seeds(f"006_{args.combo}_{source}", seeds, make_runner(ctx, source, real_sims, ac),
                  parts / source, params=params, progress=progress)


SECTION_DROPS = {"none": (), "T": ("T",), "T_C_A": ("T", "C", "A")}
SENS_TERMS = ["contrast", "contrast_adj", "h_contrast", "s_log_dist", "s_same_section", "s_same_language",
              "s_same_hand"]


def cmd_sections(args) -> None:
    """sensitivity_sections.csv: key numbers without the text-only (T) or also C/A pages (outside the spec).

    Pages whose illustration code is ``T`` lie inside other sections (e.g. f76r in the bio quires), so
    they inflate "other section, near" pairs.  Real text with page bootstrap, the three permutation
    nulls and the auto-copy null (same samples as ``pages``: seeds 10000+); cosine and Jaccard, all K.
    """
    rows = []
    srcs = ("real", "real_boot", "perm_global", "perm_section", "perm_hand", "autocopy_86")
    ac = autocopy_params()
    combos = COMBOS if args.combo == "all" else [args.combo]
    for combo in combos:
        for name, drop in SECTION_DROPS.items():
            ctx = Context(combo, with_grammar=True, drop_sections=drop)
            sims = ctx.similarity(ctx.real_d)
            res = {src: [make_runner(ctx, src, sims, ac)(sd) for sd in SOURCES[src]] for src in srcs}
            for m in METRICS:
                for k in KS:
                    for t in SENS_TERMS:
                        c = f"{m}_k{k}_{t}"
                        real = res["real"][0][c]
                        boot = pd.Series([r[c] for r in res["real_boot"]], dtype=float)
                        row = {**ctx.columns(), "dropped": name, "n_pages": ctx.n, "metric": m, "K": k, "term": t,
                               "real": real, "lo": boot.quantile(0.025), "hi": boot.quantile(0.975)}
                        for src in srcs[2:]:
                            v = pd.Series([r[c] for r in res[src]], dtype=float)
                            row[f"{src}_mean"] = v.mean()
                            row[f"{src}_p_le"] = (1 + (v <= real).sum()) / (1 + v.notna().sum())
                            row[f"{src}_p_ge"] = (1 + (v >= real).sum()) / (1 + v.notna().sum())
                        rows.append(row)
            print(f"sections {combo} {name}: {ctx.n} pages", flush=True)
    out = pd.DataFrame(rows).round(5)
    for combo, sub in out.groupby("combo", sort=False):
        sub.to_csv(OUT / "parts" / combo / "sensitivity_sections.csv", index=False)
    merge_sections()


def merge_sections() -> None:
    parts = [OUT / "parts" / c / "sensitivity_sections.csv" for c in COMBOS]
    frames = [pd.read_csv(p) for p in parts if p.exists()]
    if frames:
        pd.concat(frames, ignore_index=True).to_csv(OUT / "sensitivity_sections.csv", index=False)


def cmd_pairs(args) -> None:
    """pairs.csv: every real page pair (split policy, both alphabets) with covariates and similarities."""
    frames = []
    for key in FILES:
        ctx = Context(f"{key}_split")
        sims = ctx.similarity(ctx.real_d)
        des = pair_design(ctx, ctx.iu, ctx.ju)
        df = pd.DataFrame({
            "file": ctx.file, "space_policy": "split",
            "page_i": np.asarray(ctx.pages)[ctx.iu], "page_j": np.asarray(ctx.pages)[ctx.ju],
            "section_i": ctx.sec[ctx.iu], "section_j": ctx.sec[ctx.ju],
            "hand_i": ctx.hand[ctx.iu], "hand_j": ctx.hand[ctx.ju],
            "dist": des["dist"], "same_section": des["same_section"].astype(int),
            "same_language": des["same_language"].astype(int), "same_hand": des["same_hand"].astype(int),
        })
        for k in KS:
            df[f"cos_k{k}"] = sims[("cos", k)][ctx.iu, ctx.ju].round(4)
        df[f"jac_k{MAIN_K}"] = sims[("jac", MAIN_K)][ctx.iu, ctx.ju].round(4)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "pairs.csv", index=False)
    sus = out[out["cos_k0"] > 0.95]
    print(f"pairs.csv: {len(out)} rows; cos_k0 > 0.95: {len(sus)}")


# --------------------------------------------------------------------------- #
# labels
# --------------------------------------------------------------------------- #


def label_items(ctx: Context) -> pd.DataFrame:
    """Clean label tokens (L loci) with page, section, locus type, initial class and stem."""
    lab = filter_tokens(ctx.all_tokens, generic_type="L", clean=True)
    lab = lab[lab["token"].str.len() > 0].reset_index(drop=True)
    d = decompose_table(lab, prefixes=ctx.prefixes, suffixes=ctx.suffixes)
    cls = classify(d["token"].tolist(), ctx.key)
    return pd.DataFrame({
        "folio": d["folio"].to_numpy(), "section": d["illustration"].astype(str).to_numpy(),
        "locus_type": d["locus_type"].astype(str).to_numpy(), "token": d["token"].to_numpy(),
        "stem": d["stem"].to_numpy(), "okot": np.isin(cls, ("ok", "ot")),
    })


def label_hits(page_code: np.ndarray, stem_code: np.ndarray, sec_of_page: np.ndarray, presence: np.ndarray) -> dict:
    """Hit indicators of label items: stem in the text of the same page / of another page of the section / elsewhere."""
    same_page = presence[page_code, stem_code] > 0
    # pages of the same section carrying the stem, minus the page itself
    sec_codes, sec_inv = np.unique(sec_of_page, return_inverse=True)
    sec_presence = np.zeros((len(sec_codes), presence.shape[1]))
    np.add.at(sec_presence, sec_inv, presence > 0)
    tot = (presence > 0).sum(axis=0)
    in_sec = sec_presence[sec_inv[page_code], stem_code] - same_page
    other_sec = tot[stem_code] - sec_presence[sec_inv[page_code], stem_code]
    return {"same_page": same_page, "same_section": in_sec > 0, "elsewhere": other_sec > 0,
            "any_other_page": (in_sec + other_sec) > 0}


def cmd_labels(args) -> None:
    rows = []
    for combo in ([args.combo] if args.combo != "all" else COMBOS):
        ctx = Context(combo)
        items = label_items(ctx)
        top = ctx.top_stems(ctx.real_d, max(KS))
        for target, types in (("P", ("P",)), ("PCR", ("P", "C", "R"))):
            text = filter_tokens(ctx.all_tokens, generic_type=list(types), clean=True)
            text = ctx.decompose(text)
            pages = sorted(set(text["folio"]) | set(items["folio"]))
            pidx = {p: i for i, p in enumerate(pages)}
            sec_map = pd.concat([text[["folio", "illustration"]], ctx.all_tokens[["folio", "illustration"]]]) \
                .drop_duplicates("folio").set_index("folio")["illustration"].astype(str)
            sec_of_page = np.asarray([sec_map[p] for p in pages])
            stems = sorted(set(text["stem"]) | set(items["stem"]))
            sidx = {s: i for i, s in enumerate(stems)}
            presence = np.zeros((len(pages), len(stems)), dtype=np.int32)
            np.add.at(presence, ([pidx[p] for p in text["folio"]], [sidx[s] for s in text["stem"]]), 1)
            for k in KS:
                ex = set(top[:k])
                it = items[(items["stem"] != "") & ~items["stem"].isin(ex) & items["locus_type"].isin(LABEL_TYPES)]
                it = it.reset_index(drop=True)
                pc = np.asarray([pidx[p] for p in it["folio"]])
                sc = np.asarray([sidx[s] for s in it["stem"]])
                sec = it["section"].to_numpy()
                groups = {"all": np.ones(len(it), bool), "okot": it["okot"].to_numpy(), "other": ~it["okot"].to_numpy()}
                for t in LABEL_TYPES:
                    groups[t] = it["locus_type"].to_numpy() == t
                obs = label_hits(pc, sc, sec_of_page, presence)
                # null: permute the page assignment of label items among items of the same section
                rng = np.random.default_rng(1000 + k)
                null = {h: {g: [] for g in groups} for h in obs}
                for _ in range(LABEL_PERMS):
                    perm = np.arange(len(it))
                    for s in np.unique(sec):
                        idx = np.flatnonzero(sec == s)
                        perm[idx] = idx[rng.permutation(len(idx))]
                    hp = label_hits(pc[perm], sc, sec_of_page, presence)
                    for h in hp:
                        for g, m in groups.items():
                            null[h][g].append(hp[h][m].mean() if m.any() else np.nan)
                # page-cluster bootstrap of the observed rates
                upages = np.unique(pc)
                by_page = {p: np.flatnonzero(pc == p) for p in upages}
                boot = {h: {g: [] for g in groups} for h in obs}
                rngb = np.random.default_rng(2000 + k)
                for _ in range(LABEL_BOOT):
                    pick = np.concatenate([by_page[p] for p in rngb.choice(upages, len(upages))])
                    for h in obs:
                        for g, m in groups.items():
                            mm = m[pick]
                            boot[h][g].append(obs[h][pick][mm].mean() if mm.any() else np.nan)
                for h in obs:
                    for g, m in groups.items():
                        nv = np.asarray(null[h][g], float)
                        o = float(obs[h][m].mean()) if m.any() else np.nan
                        bv = np.asarray(boot[h][g], float)
                        rows.append({
                            **ctx.columns(), "text": target, "K": k, "group": g, "hit": h, "n_items": int(m.sum()),
                            "n_pages": int(len(np.unique(pc[m]))), "observed": o,
                            "obs_lo": np.nanquantile(bv, 0.025), "obs_hi": np.nanquantile(bv, 0.975),
                            "null_mean": np.nanmean(nv), "null_lo": np.nanquantile(nv, 0.025),
                            "null_hi": np.nanquantile(nv, 0.975),
                            "ratio": o / np.nanmean(nv) if np.nanmean(nv) > 0 else np.nan,
                            "p_upper": (1 + np.sum(nv >= o)) / (1 + len(nv)),
                        })
        print(f"labels {combo}: {len(items)} label tokens", flush=True)
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "labels.csv"
    if args.combo != "all" and path.exists():
        old = pd.read_csv(path)
        out = pd.concat([old[old["combo"] != args.combo], out], ignore_index=True)
    out.round(5).to_csv(path, index=False)
    write_meta(OUT / "meta_labels", "006_labels", {"label_types": LABEL_TYPES, "perms": LABEL_PERMS,
                                                   "boot": LABEL_BOOT, "Ks": KS})


# --------------------------------------------------------------------------- #
# ok-/ot- profile
# --------------------------------------------------------------------------- #

FEATS = ["n", "n_elig", "sum_relpos", "n_first", "n_last", "sum_len", "n_parfirst"] + \
        [f"h{i}" for i in range(5)] + [f"prev_{c}" for c in NEIGH] + [f"next_{c}" for c in NEIGH]


def profile_counts(df: pd.DataFrame, key: str, pages: list[str]) -> np.ndarray:
    """Array page x class x FEATS of counts/sums for the profile of token-initial classes."""
    cls = classify(df["token"].tolist(), key)
    neigh = np.where(np.isin(cls, ("ok", "ot")), "okot", cls)
    lid = df["locus_id"].to_numpy()
    same_prev = np.r_[False, lid[1:] == lid[:-1]]
    same_next = np.r_[lid[:-1] == lid[1:], False]
    prev = np.where(same_prev, np.concatenate([["edge"], neigh[:-1]]), "edge")
    nxt = np.where(same_next, np.concatenate([neigh[1:], ["edge"]]), "edge")
    elig = df["n_tokens"].to_numpy() >= 2
    rel = df["rel_pos"].to_numpy(float)
    hb = np.minimum((rel * 5).astype(int), 4)
    pidx = {p: i for i, p in enumerate(pages)}
    pc = np.asarray([pidx[p] for p in df["folio"]])
    cc = np.asarray([BASE_CLASSES.index(c) for c in cls])
    arr = np.zeros((len(pages), len(BASE_CLASSES), len(FEATS)))
    f = {n: i for i, n in enumerate(FEATS)}

    def add(feat, w, mask=None):
        m = np.ones(len(df), bool) if mask is None else mask
        np.add.at(arr, (pc[m], cc[m], np.full(m.sum(), f[feat])), np.asarray(w, float)[m] if np.ndim(w) else w)

    add("n", 1.0)
    add("n_elig", 1.0, elig)
    add("sum_relpos", rel, elig)
    add("n_first", 1.0, elig & df["is_first"].to_numpy(bool))
    add("n_last", 1.0, elig & df["is_last"].to_numpy(bool))
    add("sum_len", df["len_glyphs"].to_numpy(float))
    add("n_parfirst", 1.0, df["is_par_first_line"].to_numpy(bool))
    for i in range(5):
        add(f"h{i}", 1.0, elig & (hb == i))
    for c in NEIGH:
        add(f"prev_{c}", 1.0, prev == c)
        add(f"next_{c}", 1.0, nxt == c)
    return arr


def profile_stats(tot: np.ndarray) -> dict:
    """Derived statistics per class (incl. merged 'okot') from summed counts (class x FEATS)."""
    f = {n: i for i, n in enumerate(FEATS)}
    classes = {c: tot[i] for i, c in enumerate(BASE_CLASSES)}
    classes["okot"] = classes["ok"] + classes["ot"]
    hist = [f"h{i}" for i in range(5)]
    out = {}
    for c, v in classes.items():
        n, ne = v[f["n"]], v[f["n_elig"]]
        out[(c, "n")] = n
        out[(c, "share_tokens")] = n / tot[:, f["n"]].sum()
        out[(c, "mean_relpos")] = v[f["sum_relpos"]] / ne if ne else np.nan
        out[(c, "share_first")] = v[f["n_first"]] / ne if ne else np.nan
        out[(c, "share_last")] = v[f["n_last"]] / ne if ne else np.nan
        out[(c, "mean_len")] = v[f["sum_len"]] / n if n else np.nan
        out[(c, "share_parfirst")] = v[f["n_parfirst"]] / n if n else np.nan
        h = np.asarray([v[f[x]] for x in hist])
        for i in range(5):
            out[(c, f"relpos_bin{i}")] = h[i] / h.sum() if h.sum() else np.nan
        for side in ("prev", "next"):
            p = np.asarray([v[f[f"{side}_{x}"]] for x in NEIGH])
            for x, val in zip(NEIGH, p):
                out[(c, f"{side}_{x}")] = val / p.sum() if p.sum() else np.nan
    # divergences of okot from the comparison classes
    for c in ("qo", "chsh", "d", "other"):
        a, b = classes["okot"], classes[c]
        out[(c, "jsd_relpos_vs_okot")] = jsd([a[f[x]] for x in hist], [b[f[x]] for x in hist])
        for side in ("prev", "next"):
            out[(c, f"jsd_{side}_vs_okot")] = jsd([a[f[f"{side}_{x}"]] for x in NEIGH],
                                                  [b[f[f"{side}_{x}"]] for x in NEIGH])
        for s in ("mean_relpos", "share_first", "share_last", "mean_len", "share_parfirst"):
            out[(c, f"diff_{s}_okot_minus")] = out[("okot", s)] - out[(c, s)]
    return out


def cmd_okot(args) -> None:
    rows = []
    for combo in ([args.combo] if args.combo != "all" else COMBOS):
        ctx = Context(combo)
        pages = pd.unique(ctx.real["folio"].to_numpy()).tolist()
        arr = profile_counts(ctx.real, ctx.key, pages)
        real = profile_stats(arr.sum(axis=0))
        rng = np.random.default_rng(0)
        boot = [profile_stats(arr[rng.integers(len(pages), size=len(pages))].sum(axis=0)) for _ in range(OKOT_BOOT)]
        shuf = [profile_stats(profile_counts(shuffle_within_lines(ctx.real, seed=s), ctx.key, pages).sum(axis=0))
                for s in range(OKOT_SHUFFLES)]
        for key, v in real.items():
            bv = np.asarray([b[key] for b in boot], float)
            sv = np.asarray([s[key] for s in shuf], float)
            rows.append({**ctx.columns(), "cls": key[0], "stat": key[1], "real": v,
                         "real_lo": np.nanquantile(bv, 0.025), "real_hi": np.nanquantile(bv, 0.975),
                         "shuffle_mean": np.nanmean(sv), "shuffle_lo": np.nanquantile(sv, 0.025),
                         "shuffle_hi": np.nanquantile(sv, 0.975)})
        print(f"okot {combo} done", flush=True)
    out = pd.DataFrame(rows)
    path = OUT / "okot_profile.csv"
    OUT.mkdir(parents=True, exist_ok=True)
    if args.combo != "all" and path.exists():
        old = pd.read_csv(path)
        out = pd.concat([old[old["combo"] != args.combo], out], ignore_index=True)
    out.round(5).to_csv(path, index=False)
    write_meta(OUT / "meta_okot", "006_okot", {"boot": OKOT_BOOT, "shuffles": OKOT_SHUFFLES,
                                               "class_prefixes": CLASS_PREFIXES})


def cmd_initials(args) -> None:
    """initial_map.csv: ZL class versus GC initial for aligned loci (justifies CLASS_PREFIXES['GC'])."""
    z = load_tables(DATA / FILES["ZL"], space_policy="split").tokens
    g = load_tables(DATA / FILES["GC"], space_policy="split").tokens
    zl = z.groupby("locus_id", sort=False)["token"].apply(list)
    gl = g.groupby("locus_id", sort=False)["token"].apply(list)
    pairs = []
    for lid in zl.index.intersection(gl.index):
        a, b = zl[lid], gl[lid]
        if len(a) == len(b):
            pairs.extend(zip(a, b))
    zc = classify([a for a, _ in pairs], "ZL")
    gc = classify([b for _, b in pairs], "GC")
    tab = pd.crosstab(pd.Series(zc, name="ZL_class"), pd.Series(gc, name="GC_class"))
    tab.to_csv(OUT / "initial_map.csv")
    print(tab)


# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #


def interval(v: pd.Series) -> tuple[float, float, float]:
    v = v.dropna()
    return float(v.mean()), float(v.quantile(0.025)), float(v.quantile(0.975))


def cmd_summarize(args) -> None:
    frames = []
    for combo in COMBOS:
        for source in SOURCES:
            r = load_runs(OUT / "parts" / combo / source)
            if len(r):
                frames.append(r)
    runs = pd.concat(frames, ignore_index=True)
    runs.to_csv(OUT / "runs.csv", index=False)
    write_meta(OUT, "006_topic_at_distance", {"sources": {k: len(v) for k, v in SOURCES.items()}, "Ks": KS,
                                              "min_tokens": MIN_TOKENS, "far": FAR, "near": NEAR},
               n_rows=len(runs))
    curve, ols, contrast = [], [], []
    terms = ["b_intercept", "b_log_dist", "b_same_section", "b_same_language", "b_same_hand", "b_r2",
             "h_log_dist", "h_same_section", "s_log_dist", "s_same_section", "s_same_language", "s_same_hand",
             "s_size_sum", "s_size_diff", "s_r2"]
    for combo, sub in runs.groupby("combo", sort=False):
        real = sub[sub["source"] == "real"].iloc[0]
        boot = sub[sub["source"] == "real_boot"]
        for m in METRICS:
            for k in KS:
                pre = f"{m}_k{k}_"
                base = {"combo": combo, "file": real["file"], "space_policy": real["space_policy"], "metric": m, "K": k}
                for source, s in sub.groupby("source", sort=False):
                    if source == "real_boot":
                        continue
                    is_real = source == "real"
                    for same in ("same", "diff"):
                        for b in BIN_NAMES:
                            c = f"{pre}c_{same}_{b}"
                            if is_real:
                                est, lo, hi = real[c], *interval(boot[c])[1:]
                            else:
                                est, lo, hi = interval(s[c])
                            curve.append({**base, "source": source, "section": same, "bin": b,
                                          "estimate": est, "lo": lo, "hi": hi})
                    for t in terms:
                        c = pre + t
                        if is_real:
                            est, lo, hi, p, q = real[c], *interval(boot[c])[1:], np.nan, np.nan
                        else:
                            est, lo, hi = interval(s[c])
                            p = (1 + (s[c] >= real[c]).sum()) / (1 + s[c].notna().sum())
                            q = (1 + (s[c] <= real[c]).sum()) / (1 + s[c].notna().sum())
                        ols.append({**base, "source": source, "model": {"b": "spec", "h": "same_hand", "s": "spec+size"}[t[0]],
                                    "term": t.split("_", 1)[1], "estimate": est, "lo": lo, "hi": hi,
                                    "p_null_ge_real": p, "p_null_le_real": q})
                    for t in ("contrast", "h_contrast", "contrast_adj", "h_contrast_adj"):
                        c = pre + t
                        if is_real:
                            est, lo, hi, p, q = real[c], *interval(boot[c])[1:], np.nan, np.nan
                        else:
                            est, lo, hi = interval(s[c])
                            p = (1 + (s[c] >= real[c]).sum()) / (1 + s[c].notna().sum())
                            q = (1 + (s[c] <= real[c]).sum()) / (1 + s[c].notna().sum())
                        contrast.append({**base, "source": source, "pairs": "same_hand" if t.startswith("h_") else "all",
                                         "size_adjusted": t.endswith("_adj"),
                                         "estimate": est, "lo": lo, "hi": hi, "p_null_ge_real": p, "p_null_le_real": q,
                                         "n_reps": 1 if is_real else int(s[c].notna().sum())})
    pd.DataFrame(curve).round(5).to_csv(OUT / "curve.csv", index=False)
    pd.DataFrame(ols).round(5).to_csv(OUT / "ols.csv", index=False)
    pd.DataFrame(contrast).round(5).to_csv(OUT / "contrast.csv", index=False)
    stems = []
    for combo in COMBOS:
        p = OUT / "parts" / combo / "function_stems.csv"
        if p.exists():
            stems.append(pd.read_csv(p).assign(combo=combo))
    if stems:
        pd.concat(stems, ignore_index=True).to_csv(OUT / "function_stems.csv", index=False)
    merge_sections()
    plot_curve(pd.DataFrame(curve))
    if (OUT / "labels.csv").exists():
        plot_labels(pd.read_csv(OUT / "labels.csv"))
    print(f"runs.csv {len(runs)} rows; curve/ols/contrast written")


def plot_curve(curve: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = list(distance_bins(np.asarray(DISTANCE_BIN_EDGES)))
    x = np.arange(len(BIN_NAMES))
    styles = {"real": ("#1f3b73", "-", "real"), "perm_section": ("#2a9d8f", "--", "perm. in section"),
              "perm_global": ("#8a8a8a", ":", "perm. global"), "autocopy_86": ("#d1495b", "-.", "auto-copy 86")}
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=True)
    for ax, combo in zip(axes.flat, COMBOS):
        sub = curve[(curve["combo"] == combo) & (curve["metric"] == "cos") & (curve["K"] == MAIN_K)]
        for source, (color, ls, name) in styles.items():
            for same, marker in (("same", "o"), ("diff", "s")):
                s = sub[(sub["source"] == source) & (sub["section"] == same)].set_index("bin").reindex(BIN_NAMES)
                if s["estimate"].isna().all():
                    continue
                lab = f"{name}, {'same' if same == 'same' else 'other'} section"
                ax.plot(x, s["estimate"], ls, color=color, marker=marker, ms=4, lw=1.5,
                        alpha=1.0 if same == "same" else 0.6, label=lab)
                if source == "real":
                    ax.fill_between(x, s["lo"], s["hi"], color=color, alpha=0.12, lw=0)
        ax.set_title(combo.replace("_", " "), fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.grid(axis="y", alpha=0.3)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("page distance (file order)")
    for ax in axes[:, 0]:
        ax.set_ylabel(f"cosine of log(1+f) stem vectors, K={MAIN_K}")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT / "curve.png", dpi=110)
    plt.close(fig)


def plot_labels(lab: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = ["all", "okot", "other"] + list(LABEL_TYPES)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=True)
    for ax, combo in zip(axes.flat, COMBOS):
        s = lab[(lab["combo"] == combo) & (lab["text"] == "P") & (lab["K"] == MAIN_K) & (lab["hit"] == "same_page")]
        s = s.set_index("group").reindex(groups)
        x = np.arange(len(groups))
        ax.bar(x - 0.2, s["observed"], 0.38, color="#1f3b73", label="observed (page bootstrap 95 %)")
        ax.errorbar(x - 0.2, s["observed"], yerr=[s["observed"] - s["obs_lo"], s["obs_hi"] - s["observed"]],
                    fmt="none", ecolor="#333333", lw=1)
        ax.bar(x + 0.2, s["null_mean"], 0.38, color="#bbbbbb", label="labels permuted in section (95 %)")
        ax.errorbar(x + 0.2, s["null_mean"], yerr=[s["null_mean"] - s["null_lo"], s["null_hi"] - s["null_mean"]],
                    fmt="none", ecolor="#666666", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{g}\nn={int(n) if n == n else 0}" for g, n in zip(groups, s["n_items"])], fontsize=7)
        ax.set_title(combo.replace("_", " "), fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel(f"label stem in P text of same page, K={MAIN_K}")
    h, lab_ = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab_, loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT / "labels.png", dpi=110)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pages")
    p.add_argument("--combo", required=True, choices=COMBOS)
    p.add_argument("--sources", default="all")
    p.add_argument("--seeds", default=None)
    for name in ("labels", "okot"):
        q = sub.add_parser(name)
        q.add_argument("--combo", default="all", choices=COMBOS + ["all"])
    sub.add_parser("pairs")
    q = sub.add_parser("sections")
    q.add_argument("--combo", default="all", choices=COMBOS + ["all"])
    sub.add_parser("initials")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    {"pages": cmd_pages, "pairs": cmd_pairs, "labels": cmd_labels, "okot": cmd_okot,
     "initials": cmd_initials, "summarize": cmd_summarize,
     "sections": cmd_sections}[args.cmd](args)


if __name__ == "__main__":
    main()
