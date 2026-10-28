"""Experiment 012: Currier A and B as registers - change of writing convention or change of content?

Measures (``voytools.registers``), all at equal size (whole pages drawn to exactly N tokens,
N = half the smaller language of the scope) so that vocabulary measures compare like with like:

1. stem sharing (Jaccard, coverage) A vs B, compared with two halves of the same language (random
   pages, 200 draws; contiguous halves by page order), a permutation of the language label among
   the pages of the same scribe, and the ts2020 / template nulls (A/B = layout labels);
2. within-stem affix transfer matrix A -> B (suffix, prefix): row entropy, diagonal mass, best
   one-to-one relabelling (Hungarian), JSD before / after relabelling; same nulls;
3. gradient: leave-one-page-out projection of every page on the A-B suffix axis, Spearman drift
   inside A and inside B along the page order, jump at the A/B boundaries;
4. contributors to the JSD A/B (suffix, prefix, stem, glyph bigram) and the known Currier markers;
5. sections inside B (and A) vs A/B inside the herbal (JSD of stems and suffixes at equal size).

Scopes: ``all`` (paragraph text), ``hand3`` (scribe 3, the only scribe writing both), ``herbal``
(illustration H).  Decompositions (a) learned affixes, (b) + positional glyphs (011 / 008).

Usage (from the repository root)::

    python experiments/012_ab_register/run.py measure --set ZL_split                 # point + 200 draws
    python experiments/012_ab_register/run.py measure --set ZL_split --source ts2020  # 30 seeds
    python experiments/012_ab_register/run.py extras --set ZL_split                  # mapping, gradient, ...
    python experiments/012_ab_register/run.py summarize
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import registers as G  # noqa: E402
from voytools import run_seeds  # noqa: E402
from voytools._common import glyph_lists  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.nulls import as_table  # noqa: E402
from voytools._common import line_ids  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M11 = _load("exp011", "011_suffix_variation")
M9 = M11.M9

OUT = ROOT / "results" / "012_ab_register"
CACHE = ROOT / "scratch" / "012_cache"
SETS = ("ZL_split", "GC_split", "ZL_join", "GC_join")
NULL_SETS = {"ts2020": ("ZL_split", "ZL_join"), "tmpl": ("ZL_split", "GC_split")}
SCOPES = ("all", "hand3", "herbal")
N_DRAW = 200
N_NULL = 30
MIN_SHARE = 3     # stems with >= 3 tokens (sharing)
MIN_MAP = 5       # stems with >= 5 tokens in A and in B (mapping)
TARGETS = ("suffix", "prefix")
MIN_PAGE = 20     # gradient: pages with >= 20 tokens
N_PERM = 200
SEC_N = 1500      # sections: tokens per side
SEC_REPS = 100
SHARE_F = ("jaccard", "cover_xy", "cover_yx", "n_x", "n_y")
MAP_F = ("row_entropy", "diag_mass", "matched_mass", "jsd_identity", "jsd_matched", "n_moved", "n_stems")


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


class Data:
    """One corpus with page metadata, decompositions (a)/(b) and scope page pools."""

    def __init__(self, df: pd.DataFrame, ctx, rank: dict):
        self.df = df.reset_index(drop=True)
        self.page = self.df["folio"].to_numpy(dtype=object)
        meta = self.df.groupby("folio", sort=False).agg(lang=("language", "first"), hand=("hand", "first"),
                                                        sec=("illustration", "first"), n=("token", "size"))
        meta["order"] = [rank[p] for p in meta.index]
        self.meta = meta.sort_values("order")
        self.dec = {}
        for var, (pre, suf) in M11.variants(ctx).items():
            d = decompose_table(self.df, prefixes=pre, suffixes=suf)
            self.dec[var] = {c: d[c].to_numpy(dtype=object) for c in ("prefix", "stem", "suffix")}

    def pool(self, scope: str, lang: str, meta: pd.DataFrame | None = None) -> list:
        m = self.meta if meta is None else meta
        sel = m["lang"] == lang
        if scope == "hand3":
            sel &= m["hand"].astype(str) == "3"
        elif scope == "herbal":
            sel &= m["sec"] == "H"
        return m.index[sel].tolist()

    def size(self, pages) -> int:
        return int(self.meta.loc[list(pages), "n"].sum())


def scope_n(data: Data, scope: str) -> int:
    return min(data.size(data.pool(scope, "A")), data.size(data.pool(scope, "B"))) // 2


def compare(data: Data, mx: np.ndarray, my: np.ndarray, prefix: str) -> dict:
    """Sharing and mapping measures between the token sets ``mx`` and ``my`` for both variants."""
    out = {}
    for var, d in data.dec.items():
        st = d["stem"]
        r = G.stem_sharing(st[mx], st[my], MIN_SHARE)
        for f in SHARE_F:
            out[f"{prefix}__{var}__share__{f}"] = r[f]
        for t in TARGETS:
            af = d[t]
            M, n = G.transfer_matrix(st[mx], af[mx], st[my], af[my], MIN_MAP, labels=sorted(set(af)))
            fit = G.mapping_fit(M, pd.Series(af[mx]).value_counts(), pd.Series(af[my]).value_counts())
            fit["n_stems"] = n
            for f in MAP_F:
                out[f"{prefix}__{var}__{t}__{f}"] = fit[f]
    return out


def permuted_meta(meta: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Language labels permuted among the pages of the same scribe."""
    m = meta.copy()
    for _, idx in m.groupby("hand").groups.items():
        idx = list(idx)
        m.loc[idx, "lang"] = m.loc[idx, "lang"].to_numpy()[rng.permutation(len(idx))]
    return m


def draw_measures(data: Data, seed: int, with_perm: bool = True) -> dict:
    """One equal-size draw per scope: A vs B, two halves of A, two halves of B, permuted labels."""
    rng = np.random.default_rng(seed)
    out = {}
    pm = permuted_meta(data.meta, rng) if with_perm else None
    for scope in SCOPES:
        n = scope_n(data, scope)
        out[f"{scope}__N"] = n
        pa, pb = data.pool(scope, "A"), data.pool(scope, "B")
        ma, mb = G.draw_tokens(data.page, pa, n, rng), G.draw_tokens(data.page, pb, n, rng)
        out.update(compare(data, ma, mb, f"{scope}__AB"))
        for lab, pool, first in (("AA", pa, ma), ("BB", pb, mb)):
            second = G.draw_tokens(data.page, pool, n, rng, exclude=np.unique(data.page[first]))
            out.update(compare(data, first, second, f"{scope}__{lab}"))
        if with_perm:
            qa, qb = data.pool(scope, "A", pm), data.pool(scope, "B", pm)
            if qa and qb:
                k = min(n, data.size(qa), data.size(qb))
                out.update(compare(data, G.draw_tokens(data.page, qa, k, rng), G.draw_tokens(data.page, qb, k, rng),
                                   f"{scope}__PERM"))
    return out


def point_measures(data: Data) -> dict:
    """Full A vs full B (unequal sizes) and contiguous halves (page order) at N tokens."""
    out = {}
    for scope in SCOPES:
        n = scope_n(data, scope)
        pa, pb = data.pool(scope, "A"), data.pool(scope, "B")
        out.update(compare(data, np.isin(data.page, pa), np.isin(data.page, pb), f"{scope}__ABfull"))
        for lab, pool in (("AAcontig", pa), ("BBcontig", pb)):
            first, second = G.split_pages(pool, data.meta.loc[pool, "n"].tolist())
            m1 = _first_n(data.page, first, n)
            m2 = _first_n(data.page, second, n)
            out.update(compare(data, m1, m2, f"{scope}__{lab}"))
    return out


def _first_n(page: np.ndarray, pages: list, n: int) -> np.ndarray:
    """Tokens of ``pages`` in page order, cut to ``n``."""
    mask = np.zeros(len(page), dtype=bool)
    got = 0
    for p in pages:
        idx = np.flatnonzero(page == p)[: max(n - got, 0)]
        mask[idx] = True
        got += len(idx)
    return mask


def gradient(data: Data, var: str = "a", target: str = "suffix") -> tuple[pd.DataFrame, dict]:
    """Page scores on the A-B axis (leave-one-page-out), drift inside A and B, boundary jumps."""
    d = data.dec[var]
    sc = G.page_scores(data.page, d[target], data.df["language"].astype(str).to_numpy(), min_tokens=MIN_PAGE)
    sc = sc.set_index("page").join(data.meta[["hand", "sec", "order"]]).reset_index().sort_values("order")
    stats = {}
    for lab in ("A", "B"):
        s = sc[sc["label"] == lab]
        r = G.order_drift(s["score"], s["order"], N_PERM, seed=1)
        stats.update({f"drift_{lab}_rho": r["rho"], f"drift_{lab}_p": r["p"], f"n_pages_{lab}": r["n_pages"]})
    j = G.boundary_jumps(sc["score"].to_numpy(), sc["label"].to_numpy())
    stats.update({f"jump_{k}": v for k, v in j.items()})
    a, b = sc[sc["label"] == "A"]["score"], sc[sc["label"] == "B"]["score"]
    stats.update({"mean_A": float(a.mean()), "mean_B": float(b.mean()), "sd_A": float(a.std()), "sd_B": float(b.std()),
                  "separation": float((b.mean() - a.mean()) / np.sqrt((a.var() + b.var()) / 2)),
                  "misclassified": int((a > 0).sum() + (b < 0).sum())})
    return sc, stats


# --------------------------------------------------------------------------- #
# sources
# --------------------------------------------------------------------------- #


def cached_corpus(ctx, kind: str, tag: str, seed: int) -> pd.DataFrame:
    """Null corpus of 009's generators, cached as token lines under scratch/ (gitignored)."""
    path = CACHE / tag / f"{kind}_{seed}.txt"
    real = ctx.real
    if path.exists():
        lines = [ln.split(" ") if ln else [] for ln in path.read_text(encoding="utf-8").split("\n")]
        return as_table(lines, real, glyph_preset=real.attrs.get("glyph_preset"))
    x = ctx.corpus(kind, seed)
    lines = x.groupby(line_ids(x), sort=True)["token"].apply(list).tolist()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(" ".join(ln) for ln in lines), encoding="utf-8")
    return x


def cmd_measure(args) -> None:
    tag, source = args.set, args.source
    ctx = M9.Ctx(tag) if source == "real" else M9.Ctx(f"real:{tag}")
    rank = ctx.rank
    if source == "real":
        data = Data(ctx.real, ctx, rank)

        def one(seed: int) -> dict:
            if seed < 0:
                return {"set": tag, "source": source, "kind": "point", **point_measures(data)}
            return {"set": tag, "source": source, "kind": "draw", **draw_measures(data, seed)}
        default = list(range(N_DRAW))
        point_dir = OUT / "runs" / f"point__{tag}"  # the point row has other columns than the draws
        run_seeds("012_ab_register", [-1], one, point_dir, params={"set": tag, "source": source, "kind": "point"})
    else:
        def one(seed: int) -> dict:
            data = Data(cached_corpus(ctx, source, tag, seed), ctx, rank)
            _, gst = gradient(data)
            return {"set": tag, "source": source, "kind": "null", **draw_measures(data, seed, with_perm=False),
                    **{f"grad__{k}": v for k, v in gst.items()}}
        default = list(range(N_NULL))
    seeds = M9.parse_seeds(args.seeds, default)
    params = {"set": tag, "source": source, "variants": {k: {"prefixes": v[0], "suffixes": v[1]}
                                                         for k, v in M11.variants(ctx).items()},
              "scopes": SCOPES, "n_draw": N_DRAW, "n_null": N_NULL, "min_share": MIN_SHARE, "min_map": MIN_MAP,
              "ts2020_config": M9.TS2020_CONFIG, "tmpl_seed0": M9.TMPL_SEED0}
    run_seeds("012_ab_register", seeds, one, OUT / "runs" / f"{source}__{tag}", params=params, progress=M11.progress)


# --------------------------------------------------------------------------- #
# extras: mapping matrices, transfers, gradient, contributors, sections
# --------------------------------------------------------------------------- #


def mapping_matrices(data: Data, var: str, target: str, n_rep: int, seed: int) -> pd.DataFrame:
    """Mean normalised transfer matrices over equal-size draws: A->B and halves (A|A, B|B pooled)."""
    rng = np.random.default_rng(seed)
    d = data.dec[var]
    labels = sorted(set(d[target]))
    acc = {"AB": 0.0, "halves": 0.0}
    n = scope_n(data, "all")
    pa, pb = data.pool("all", "A"), data.pool("all", "B")
    for _ in range(n_rep):
        ma, mb = G.draw_tokens(data.page, pa, n, rng), G.draw_tokens(data.page, pb, n, rng)
        M, _ = G.transfer_matrix(d["stem"][ma], d[target][ma], d["stem"][mb], d[target][mb], MIN_MAP, labels)
        acc["AB"] = acc["AB"] + M.to_numpy() / max(M.to_numpy().sum(), 1e-12) / n_rep
        for pool, first in ((pa, ma), (pb, mb)):
            second = G.draw_tokens(data.page, pool, n, rng, exclude=np.unique(data.page[first]))
            H, _ = G.transfer_matrix(d["stem"][first], d[target][first], d["stem"][second], d[target][second],
                                     MIN_MAP, labels)
            acc["halves"] = acc["halves"] + H.to_numpy() / max(H.to_numpy().sum(), 1e-12) / (2 * n_rep)
    rows = []
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            rows.append({"target": target, "variant": var, "affix_A": a or "<none>", "affix_B": b or "<none>",
                         "mass_AB": acc["AB"][i, j], "mass_halves": acc["halves"][i, j],
                         "excess": acc["AB"][i, j] - acc["halves"][i, j],
                         "row_share_AB": acc["AB"][i, j] / acc["AB"][i].sum() if acc["AB"][i].sum() else np.nan})
    return pd.DataFrame(rows)


def contributors(data: Data, scope: str, gpreset) -> pd.DataFrame:
    rows = []
    pa, pb = data.pool(scope, "A"), data.pool(scope, "B")
    ma, mb = np.isin(data.page, pa), np.isin(data.page, pb)
    glyphs = glyph_lists(data.df, gpreset)
    big = np.asarray([" ".join(f"{x}{y}" for x, y in zip(("<",) + g, g + (">",))) for g in glyphs], dtype=object)
    units = {"suffix": data.dec["a"]["suffix"], "prefix": data.dec["a"]["prefix"], "stem": data.dec["a"]["stem"],
             "token": data.df["token"].to_numpy(dtype=object)}
    for unit, vals in list(units.items()) + [("glyph_bigram", None)]:
        if vals is None:
            ca = pd.Series([b for s in big[ma] for b in s.split(" ")]).value_counts()
            cb = pd.Series([b for s in big[mb] for b in s.split(" ")]).value_counts()
        else:
            ca, cb = pd.Series(vals[ma]).value_counts(), pd.Series(vals[mb]).value_counts()
        idx = ca.index.union(cb.index)
        ca, cb = ca.reindex(idx, fill_value=0), cb.reindex(idx, fill_value=0)
        c = G.jsd_contributions(ca.to_numpy(), cb.to_numpy())
        tot = c.sum()
        order = np.argsort(-c)
        for rnk, k in enumerate(order[:10]):
            rows.append({"scope": scope, "unit": unit, "rank": rnk + 1, "value": idx[k] or "<none>",
                         "p_A": ca.iloc[k] / ca.sum(), "p_B": cb.iloc[k] / cb.sum(), "contribution": c[k],
                         "share_of_jsd": c[k] / tot, "jsd_total": tot})
    return pd.DataFrame(rows)


def sections_within(data: Data, seed: int) -> pd.DataFrame:
    """JSD of stems and suffixes between groups at SEC_N tokens each (token subsamples), and the noise floor."""
    rng = np.random.default_rng(seed)
    m = data.meta
    groups = {}
    for lang in ("A", "B"):
        for sec in sorted(m.loc[m["lang"] == lang, "sec"].unique()):
            pages = m.index[(m["lang"] == lang) & (m["sec"] == sec)]
            idx = np.flatnonzero(np.isin(data.page, pages))
            if len(idx) >= SEC_N:
                groups[f"{sec}_{lang}"] = idx
    names = sorted(groups)
    rows = []
    for unit in ("stem", "suffix"):
        v = data.dec["a"][unit]
        for i, g1 in enumerate(names):
            for g2 in names[i:]:
                vals = []
                for _ in range(SEC_REPS):
                    if g1 == g2:  # noise floor: two disjoint token samples of the same group
                        if len(groups[g1]) < 2 * SEC_N:
                            break
                        perm = rng.permutation(groups[g1])
                        x, y = perm[:SEC_N], perm[SEC_N: 2 * SEC_N]
                    else:
                        x = rng.choice(groups[g1], SEC_N, replace=False)
                        y = rng.choice(groups[g2], SEC_N, replace=False)
                    cx, cy = pd.Series(v[x]).value_counts(), pd.Series(v[y]).value_counts()
                    idx = cx.index.union(cy.index)
                    vals.append(G.jsd(cx.reindex(idx, fill_value=0).to_numpy(), cy.reindex(idx, fill_value=0).to_numpy()))
                if not vals:
                    continue
                kind = ("same_group" if g1 == g2 else
                        "A_vs_B_herbal" if {g1, g2} == {"H_A", "H_B"} else
                        "within_B" if g1.endswith("_B") and g2.endswith("_B") else
                        "within_A" if g1.endswith("_A") and g2.endswith("_A") else "A_vs_B_other")
                rows.append({"unit": unit, "group_1": g1, "group_2": g2, "kind": kind, "jsd": float(np.mean(vals)),
                             "lo": float(np.quantile(vals, 0.025)), "hi": float(np.quantile(vals, 0.975)),
                             "n_tokens_each": SEC_N, "reps": len(vals)})
    return pd.DataFrame(rows)


def cmd_extras(args) -> None:
    tag = args.set
    ctx = M9.Ctx(tag)
    data = Data(ctx.real, ctx, ctx.rank)
    key = tag.split("_")[0]
    out = OUT / "extras" / tag
    out.mkdir(parents=True, exist_ok=True)
    mats = pd.concat([mapping_matrices(data, var, t, N_DRAW // 4, seed=7) for var in data.dec for t in TARGETS],
                     ignore_index=True)
    mats.round(6).to_csv(out / "mapping.csv", index=False)
    sc, st = gradient(data)
    sc.round(5).to_csv(out / "page_scores.csv", index=False)
    rows = [{"set": tag, "var": "a", "target": "suffix", **st}]
    for var in data.dec:
        for t in TARGETS:
            if (var, t) != ("a", "suffix"):
                rows.append({"set": tag, "var": var, "target": t, **gradient(data, var, t)[1]})
    pd.DataFrame(rows).round(6).to_csv(out / "gradient.csv", index=False)
    pd.concat([contributors(data, s, ctx.gpreset) for s in SCOPES], ignore_index=True).round(6).to_csv(
        out / "contributors.csv", index=False)
    sections_within(data, seed=11).round(6).to_csv(out / "sections_within.csv", index=False)
    print(f"extras {tag}: done ({key})", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("measure")
    p.add_argument("--set", required=True, choices=SETS)
    p.add_argument("--source", default="real", choices=("real", "ts2020", "tmpl"))
    p.add_argument("--seeds", default=None)
    p = sub.add_parser("extras")
    p.add_argument("--set", required=True, choices=SETS)
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from summarize import summarize  # noqa: E402  (experiments/012_ab_register/summarize.py)

        summarize()
        return
    {"measure": cmd_measure, "extras": cmd_extras}[args.cmd](args)


if __name__ == "__main__":
    main()
