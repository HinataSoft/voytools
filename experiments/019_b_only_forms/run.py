"""Experiment 019: what is the fifth of B tokens that A does not have - names, other morphology, or a second vocabulary?

Whole-token forms of B (A symmetrically) are classed against the vocabulary of the other language
(``voytools.bforms``): shared (>= 1 token in A), B-only (>= 2 in B, 0 in A) split into d-variants
(a counterpart in A reached by one Currier operation: ins ``d`` / ``e`` / ``ee`` / ``q``, ``ch`` <-> ``sh``;
v101 ``8`` / ``7`` / ``c`` / ``C`` / ``4``, ``1`` <-> ``2``), other counterpart within one edit, no counterpart
within one edit; and hapaxes of B.  Stems of decomposition (a); a new stem has no token in A.

1. composition of the fifth (token shares of the classes, new stems vs new forms of shared stems);
2. profile of every class against shared forms matched by token count (paradigm width, dispersion
   over pages, section entropy, labels, line edges, co-occurrence in one line, glyph-bigram
   log-likelihood under a model of the shared forms, seam) - in A x B, in halves of B (noise) and A x B
   in the other direction;
3. the ``-ey`` / ``-edy`` choice: P(``-edy`` | glyph before, prefix, line position, last glyph of the
   previous token) in B and in A, entropy with and without the context, stems with both forms (and on
   one page), against ``-y`` / ``-dy`` and ``-ol`` / ``-or`` as the yardstick of ordinary variation;
4. Latin references (Pliny other topic, Caesar vs Celsus other register) and generators (ts2020,
   template) through the same code;
5. lists of the most frequent new stems of B and of A.

Usage (from the repository root)::

    python experiments/019_b_only_forms/run.py measure --set ZL_split                  # 200 draws
    python experiments/019_b_only_forms/run.py measure --set ZL_split --source ts2020  # 30 seeds
    python experiments/019_b_only_forms/run.py measure --set pliny                     # control
    python experiments/019_b_only_forms/run.py edy --set ZL_split                      # -ey / -edy rule
    python experiments/019_b_only_forms/run.py lists --set ZL_split                    # new stems, classes
    python experiments/019_b_only_forms/run.py summarize
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import bforms as F  # noqa: E402
from voytools import registers as G  # noqa: E402
from voytools import run_seeds  # noqa: E402
from voytools._common import line_ids  # noqa: E402
from voytools.glyphs import collapse_eva  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.nulls import as_table  # noqa: E402
from voytools.table import filter_tokens, load_tables  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M16 = _load("exp016", "016_ab_stem_variants")
M12 = M16.M12
M11 = M16.M11
M9 = M16.M9

OUT = ROOT / "results" / "019_b_only_forms"
CACHE = ROOT / "scratch" / "019_cache"
SETS = ("ZL_split", "GC_split", "ZL_join", "GC_join")
CTL_SETS = ("pliny", "caesar_celsus")
NULL_SETS = {"ts2020": ("ZL_split", "ZL_join"), "tmpl": ("ZL_split", "GC_split")}
N_DRAW = 200
N_NULL = 30
MIN_COUNT = 2
MAX_DIST = 2
N_SHUF = 20          # co-occurrence shuffles per profile
N_PERM = 200         # context permutations (-edy rule)
MIN_CTX = 5          # tokens per context for the determinism shares
PLINY_A = M16.PLINY_A
CMPS = (("AB", "y_only"), ("AB", "x_only"), ("BB", "x_only"), ("AA", "x_only"))
PROFILE_CLASSES = ("d_variant", "other_le1", "no_counterpart", "hapax")
COMP_F = ("n_tokens", "n_forms", "share_only", "share_shared", "share_d_variant", "share_other_le1",
          "share_no_counterpart", "share_hapax", "only_d_variant", "only_other_le1", "only_no_counterpart",
          "only_new_stem", "only_new_form", "d_variant_new_stem", "other_le1_new_stem", "no_counterpart_new_stem",
          "hapax_new_stem", "n_forms_d_variant", "n_forms_other_le1", "n_forms_no_counterpart", "n_forms_hapax",
          "n_forms_shared")
PROFILE_F = ("n_forms", "n_tokens", "paradigm_forms", "H_suffix_given_stem", "pages_per_token", "gries_dp",
             "top_page_share", "section_entropy", "label_share", "label_any", "first_share", "last_share",
             "par_first_share", "cooc_ratio", "bigram_ll", "seam_jsd")
# suffix pairs: (name, s1, s2, exclude base ending) per scheme
PAIRS = {
    "extended": [("ey_edy", "ey", "edy", None), ("y_dy", "y", "dy", "e"), ("ol_or", "ol", "or", None)],
    "v101": [("c9_c89", "c9", "c89", None), ("c9_c79", "c9", "c79", None), ("9_89", "9", "89", "c"),
             ("oe_oy", "oe", "oy", None)],
    None: [("us_um", "us", "um", None), ("is_es", "is", "es", None)],
}
CONTEXTS = ("prev_unit", "prefix", "line_pos", "prev_token_last", "prev_unit+prefix")


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


def _units(strings, scheme) -> np.ndarray:
    cache: dict = {}
    out = np.empty(len(strings), dtype=object)
    for i, s in enumerate(strings):
        t = cache.get(s)
        if t is None:
            t = cache[s] = (tuple(collapse_eva(s, "extended")) if scheme == "extended" else tuple(s)) if s else ()
        out[i] = t
    return out


def scheme_of(df: pd.DataFrame, ctx) -> str | None:
    if ctx.gpreset == "extended":
        return "extended"
    alph = str(df.attrs.get("alphabet", "")).lower()
    return "v101" if "v101" in alph else None


class Data:
    """One corpus: whole-token forms (unit tuples), stems / suffixes of (a), per-token metadata, label counts."""

    def __init__(self, df: pd.DataFrame, ctx, rank: dict, labels: pd.DataFrame | None = None):
        self.df = df.reset_index(drop=True)
        self.page = self.df["folio"].to_numpy(dtype=object)
        meta = self.df.groupby("folio", sort=False).agg(lang=("language", "first"), hand=("hand", "first"),
                                                        sec=("illustration", "first"), n=("token", "size"))
        meta["order"] = [rank.get(p, i) for i, p in enumerate(meta.index)]
        self.meta = meta.sort_values("order")
        self.scheme = scheme_of(self.df, ctx)
        pre, suf = M11.variants(ctx)["a"]
        self.prefixes, self.suffixes = pre, suf
        d = decompose_table(self.df, prefixes=pre, suffixes=suf)
        self.stem = d["stem"].to_numpy(dtype=object)
        self.suffix = d["suffix"].to_numpy(dtype=object)
        self.prefix = d["prefix"].to_numpy(dtype=object)
        tokens = self.df["token"].astype(str).tolist()
        self.tokens = np.asarray(tokens, dtype=object)
        self.form = _units(tokens, self.scheme)
        self.code, uniq = pd.factorize(pd.Series(self.form, dtype=object))
        self.uniq = list(uniq)
        self.code_of = {f: i for i, f in enumerate(self.uniq)}
        self.stem_of = {}
        for f, s in zip(self.form, self.stem):
            if f not in self.stem_of:
                self.stem_of[f] = s
        self.line = line_ids(self.df)
        self.section = self.df["illustration"].astype(str).to_numpy(dtype=object)
        self.is_first = self.df["is_first"].to_numpy(dtype=bool)
        self.is_last = self.df["is_last"].to_numpy(dtype=bool)
        self.is_par_first = self.df["is_par_first_line"].fillna(False).to_numpy(dtype=bool)
        nxt = np.empty(len(self.form), dtype=object)
        prev_last = np.empty(len(self.form), dtype=object)
        for i in range(len(self.form)):
            nxt[i] = self.form[i + 1][0] if i + 1 < len(self.form) and self.line[i + 1] == self.line[i] and self.form[i + 1] else ""
            prev_last[i] = self.form[i - 1][-1] if i > 0 and self.line[i - 1] == self.line[i] and self.form[i - 1] else "#"
        self.next_first = nxt
        self.prev_token_last = prev_last
        self.vocab = sorted({u for f in self.uniq for u in f})
        self.label_counts = None
        if labels is not None:
            lang_of = dict(zip(self.meta.index, self.meta["lang"]))
            lab = labels[labels["folio"].map(lang_of).notna()]
            lf = _units(lab["token"].astype(str).tolist(), self.scheme)
            ll = lab["folio"].map(lang_of).to_numpy(dtype=object)
            self.label_counts = {}
            for lang in ("A", "B"):
                cl = Counter(lf[ll == lang].tolist())
                cp = Counter(self.form[self.df["language"].to_numpy() == lang].tolist())
                self.label_counts[lang] = {f: (cl.get(f, 0), cp.get(f, 0)) for f in set(cl) | set(cp)}

    def pool(self, lang: str) -> list:
        return self.meta.index[self.meta["lang"] == lang].tolist()

    def size(self, pages) -> int:
        return int(self.meta.loc[list(pages), "n"].sum())

    def counts(self, mask: np.ndarray) -> dict:
        c = Counter(self.form[mask].tolist())
        c.pop((), None)
        return dict(c)

    def tok(self, mask: np.ndarray) -> dict:
        return {"form": self.form[mask], "form_code": self.code[mask], "form_uniq": self.uniq,
                "stem": self.stem[mask], "suffix": self.suffix[mask], "page": self.page[mask],
                "section": self.section[mask], "line": self.line[mask], "is_first": self.is_first[mask],
                "is_last": self.is_last[mask], "is_par_first": self.is_par_first[mask],
                "next_first": self.next_first[mask]}


def scope_n(data: Data) -> int:
    return min(data.size(data.pool("A")), data.size(data.pool("B"))) // 2


# --------------------------------------------------------------------------- #
# draw measures
# --------------------------------------------------------------------------- #


def one_side(data: Data, mx: np.ndarray, my: np.ndarray, lang_x: str | None, rng: np.random.Generator,
             profiles: bool = True) -> tuple[dict, pd.DataFrame]:
    """Classes of the X sample against the Y sample, their composition and (optionally) class profiles."""
    cx, cy = data.counts(mx), data.counts(my)
    cl = F.classify_forms(cx, cy, data.stem_of, data.scheme, MIN_COUNT, MAX_DIST)
    out = {f"comp__{k}": v for k, v in F.composition(cl).items() if k in COMP_F}
    if not profiles:
        return out, cl
    tok = data.tok(mx)
    labels = data.label_counts.get(lang_x) if (data.label_counts and lang_x) else None
    sh = cl["cls"] == "shared"
    shared = dict(zip(cl.loc[sh, "form"], cl.loc[sh, "n_x"]))
    seam_ref = np.isin(tok["form_code"], np.asarray([data.code_of[f] for f in shared], dtype=np.int64))
    for cls in PROFILE_CLASSES:
        forms, matched = F.match_by_count(cl, cls, rng)
        if not forms:
            continue
        held = {f: c for f, c in shared.items() if f not in set(matched)}
        model = F.bigram_model(held, data.vocab)
        p = F.form_profile(forms, tok, labels, model, rng, N_SHUF, seam_ref)
        q = F.form_profile(matched, tok, labels, model, rng, N_SHUF, seam_ref)
        out.update({f"{cls}__{k}": v for k, v in p.items() if k in PROFILE_F})
        out.update({f"{cls}_m__{k}": v for k, v in q.items() if k in PROFILE_F})
    return out, cl


def draw_measures(data: Data, seed: int, profiles: bool = True) -> dict:
    """One equal-size draw: A vs B both directions, two halves of B, two halves of A."""
    rng = np.random.default_rng(seed)
    n = scope_n(data)
    out = {"N": n}
    pa, pb = data.pool("A"), data.pool("B")
    ma, mb = G.draw_tokens(data.page, pa, n, rng), G.draw_tokens(data.page, pb, n, rng)
    r, _ = one_side(data, mb, ma, "B", rng, profiles)
    out.update({f"AB__y_only__{k}": v for k, v in r.items()})
    r, _ = one_side(data, ma, mb, "A", rng, profiles)
    out.update({f"AB__x_only__{k}": v for k, v in r.items()})
    for lab, pool, first, lang in (("BB", pb, mb, "B"), ("AA", pa, ma, "A")):
        second = G.draw_tokens(data.page, pool, n, rng, exclude=np.unique(data.page[first]))
        r, _ = one_side(data, first, second, lang, rng, profiles)
        out.update({f"{lab}__x_only__{k}": v for k, v in r.items()})
    return out


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


def control_frame(name: str) -> tuple[pd.DataFrame, object]:
    """Latin controls as two 'languages' (as 016): Pliny books 16-21 vs 23-30, Caesar vs Celsus."""
    if name == "pliny":
        df = M9.control_table("pliny_prose").reset_index(drop=True)
        df["language"] = np.where(df["illustration"].astype(str).isin(PLINY_A), "A", "B")
    elif name == "caesar_celsus":
        parts = []
        for lab, ctl in (("A", "caesar_prose"), ("B", "celsus_prose")):
            d = M9.control_table(ctl).reset_index(drop=True)
            d["language"] = lab
            d["folio"] = lab + ":" + d["folio"].astype(str)
            parts.append(d)
        df = pd.concat(parts, ignore_index=True)
        df.attrs.update(parts[0].attrs)
    else:
        raise ValueError(name)
    df["hand"] = "ctl"
    dec = decompose_table(df)
    ctx = SimpleNamespace(source=f"ctl:{name}", base_tag=name, gpreset=None, real=df,
                          prefixes=dec.attrs["prefixes"], suffixes=dec.attrs["suffixes"])
    return df, ctx


def control_data(name: str) -> Data:
    df, ctx = control_frame(name)
    rank = {p: i for i, p in enumerate(pd.unique(df["folio"].to_numpy()))}
    return Data(df, ctx, rank)


def labels_table(tag: str) -> pd.DataFrame:
    key, policy = tag.split("_")
    t = load_tables(ROOT / "data" / M9.FILES[key], space_policy=policy)
    return filter_tokens(t.tokens, generic_type="L", clean=True).reset_index(drop=True)


def real_data(tag: str) -> tuple:
    ctx = M9.Ctx(tag)
    return Data(ctx.real, ctx, ctx.rank, labels=labels_table(tag)), ctx


def cmd_measure(args) -> None:
    tag, source = args.set, args.source
    if tag in CTL_SETS:
        data = control_data(tag)

        def one(seed: int) -> dict:
            return {"set": tag, "source": "control", "kind": "draw", **draw_measures(data, seed)}
        default = list(range(N_DRAW))
        extra = {"pliny_A": PLINY_A}
    elif source == "real":
        data, ctx = real_data(tag)

        def one(seed: int) -> dict:
            return {"set": tag, "source": source, "kind": "draw", **draw_measures(data, seed)}
        default = list(range(N_DRAW))
        extra = {"prefixes": list(data.prefixes), "suffixes": list(data.suffixes)}
    else:
        ctx = M9.Ctx(f"real:{tag}")

        def one(seed: int) -> dict:
            data = Data(cached_corpus(ctx, source, tag, seed), ctx, ctx.rank)
            return {"set": tag, "source": source, "kind": "null", **draw_measures(data, seed)}
        default = list(range(N_NULL))
        extra = {"ts2020_config": M9.TS2020_CONFIG, "tmpl_seed0": M9.TMPL_SEED0}
    seeds = M9.parse_seeds(args.seeds, default)
    params = {"set": tag, "source": source, "n_draw": N_DRAW, "n_null": N_NULL, "min_count": MIN_COUNT,
              "max_dist": MAX_DIST, "n_shuf": N_SHUF, "classes": list(F.CLASSES), "d_ops": {k: {"units": sorted(v["units"]),
              "subs": [sorted(s) for s in v["subs"]]} for k, v in F.D_OPS.items()}, **extra}
    run_seeds("019_b_only_forms", seeds, one, OUT / "runs" / f"{source}__{tag}", params=params, progress=M11.progress)


# --------------------------------------------------------------------------- #
# the -ey / -edy rule (full text)
# --------------------------------------------------------------------------- #


def edy_rule(data: Data, tag: str, source: str = "real", n_perm: int = N_PERM) -> tuple[pd.DataFrame, pd.DataFrame]:
    """P(second form | context) for every suffix pair, language and context; summary and context tables."""
    rows, ctx_rows = [], []
    lang = data.df["language"].to_numpy(dtype=object)
    line_pos = np.where(data.is_first, "first", np.where(data.is_last, "last", "mid")).astype(object)
    pairs = PAIRS.get(data.scheme, PAIRS[None])
    for name, s1, s2, excl in pairs:
        sp = F.suffix_pairs(data.tokens, s1, s2, excl)
        if sp.empty:
            continue
        idx = sp.index.to_numpy()
        base_units = _units(sp["base"].tolist(), data.scheme)
        prev_unit = np.asarray([b[-1] if b else "#" for b in base_units], dtype=object)
        ctx_vals = {"prev_unit": prev_unit, "prefix": data.prefix[idx].astype(object), "line_pos": line_pos[idx],
                    "prev_token_last": data.prev_token_last[idx],
                    "prev_unit+prefix": np.asarray([f"{a}|{b}" for a, b in zip(prev_unit, data.prefix[idx])], dtype=object)}
        for lg in ("A", "B"):
            sel = lang[idx] == lg
            ch = sp["choice"].to_numpy()[sel]
            base = {"set": tag, "source": source, "pair": name, "language": lg, "n": int(sel.sum()),
                    "n_s1": int((ch == 0).sum()), "n_s2": int((ch == 1).sum())}
            bf = F.both_forms(sp["base"].to_numpy()[sel], ch, data.page[idx][sel], MIN_COUNT)
            for cname in CONTEXTS:
                r = F.pair_choice(ch, ctx_vals[cname][sel], MIN_CTX, n_perm, seed=1)
                rows.append({**base, "context": cname, **{k: v for k, v in r.items() if k != "n"}, **bf})
                if cname in ("prev_unit", "prefix", "line_pos", "prev_token_last") and sel.any():
                    d = pd.DataFrame({"c": ctx_vals[cname][sel], "ch": ch})
                    g = d.groupby("c")["ch"].agg(["size", "mean"])
                    for c, (n, m) in g.iterrows():
                        ctx_rows.append({"set": tag, "source": source, "pair": name, "language": lg, "context": cname,
                                         "value": c, "n": int(n), "p_s2": float(m)})
    return pd.DataFrame(rows), pd.DataFrame(ctx_rows)


def cmd_edy(args) -> None:
    tag, source = args.set, args.source
    out = OUT / "extras" / f"{source}__{tag}"
    out.mkdir(parents=True, exist_ok=True)
    if tag in CTL_SETS:
        data = control_data(tag)
        rule, ctxs = edy_rule(data, tag, "control")
    elif source == "real":
        data, _ = real_data(tag)
        rule, ctxs = edy_rule(data, tag)
    else:
        ctx = M9.Ctx(f"real:{tag}")
        parts, cparts = [], []
        for seed in range(N_NULL):
            data = Data(cached_corpus(ctx, source, tag, seed), ctx, ctx.rank)
            r, c = edy_rule(data, tag, source, n_perm=20)
            r.insert(0, "seed", seed)
            parts.append(r)
            if seed == 0:
                cparts.append(c)
            print(f"  edy {source} seed {seed}", flush=True)
        rule, ctxs = pd.concat(parts, ignore_index=True), pd.concat(cparts, ignore_index=True)
    rule.round(6).to_csv(out / "edy_rule.csv", index=False)
    ctxs.round(6).to_csv(out / "edy_contexts.csv", index=False)
    print(f"edy {source} {tag}: done", flush=True)


# --------------------------------------------------------------------------- #
# lists: classes and new stems on the full text
# --------------------------------------------------------------------------- #


def _s(units) -> str:
    return "".join(units)


def full_classes(data: Data, tag: str) -> pd.DataFrame:
    """Every form of B against full A (and A against full B): class, counterpart, operation, stem."""
    rows = []
    ma = np.isin(data.page, data.pool("A"))
    mb = np.isin(data.page, data.pool("B"))
    for dname, mx, my, lang in (("B_only", mb, ma, "B"), ("A_only", ma, mb, "A")):
        cl = F.classify_forms(data.counts(mx), data.counts(my), data.stem_of, data.scheme, MIN_COUNT, MAX_DIST)
        cl = cl[cl["cls"] != "shared"].copy()
        pages = pd.DataFrame({"f": data.code[mx], "p": data.page[mx], "s": data.section[mx]})
        np_ = pages.groupby("f")["p"].nunique()
        secs = pages.groupby("f")["s"].agg(lambda v: ";".join(sorted(set(v))))
        labels = data.label_counts.get(lang, {}) if data.label_counts else {}
        code_of = {f: i for i, f in enumerate(data.uniq)}
        for r in cl.itertuples(index=False):
            c = code_of[r.form]
            rows.append({"set": tag, "direction": dname, "form": _s(r.form), "cls": r.cls, "tokens": int(r.n_x),
                         "counterpart": _s(r.counterpart) if r.counterpart else "", "distance": r.distance,
                         "op": r.op, "stem": r.stem, "stem_new": r.stem_new, "n_pages": int(np_.get(c, 0)),
                         "sections": secs.get(c, ""), "labels": int(labels.get(r.form, (0, 0))[0])})
    return pd.DataFrame(rows).sort_values(["direction", "tokens", "form"], ascending=[True, False, True])


def top_new_stems(data: Data, tag: str, k: int = 30) -> pd.DataFrame:
    """The most frequent stems of B with no token in A (and of A with no token in B)."""
    rows = []
    lang = data.df["language"].to_numpy(dtype=object)
    for lg, other in (("B", "A"), ("A", "B")):
        mx, my = lang == lg, lang == other
        sy = set(data.stem[my].tolist())
        d = pd.DataFrame({"stem": data.stem[mx], "form": data.tokens[mx], "page": data.page[mx],
                          "sec": data.section[mx], "first": data.is_first[mx], "last": data.is_last[mx]})
        d = d[(d["stem"] != "") & ~d["stem"].isin(sy)]
        labels = data.label_counts.get(lg, {}) if data.label_counts else {}
        lab_by_stem: dict = {}
        for f, (nl, _) in labels.items():
            st = data.stem_of.get(f)
            if st is not None and nl:
                lab_by_stem[st] = lab_by_stem.get(st, 0) + nl
        g = d.groupby("stem").agg(tokens=("form", "size"), n_forms=("form", "nunique"), n_pages=("page", "nunique"),
                                  forms=("form", lambda v: ";".join(f"{a}:{b}" for a, b in Counter(v).most_common(6))),
                                  sections=("sec", lambda v: ";".join(f"{a}:{b}" for a, b in Counter(v).most_common())),
                                  pages=("page", lambda v: ";".join(f"{a}:{b}" for a, b in Counter(v).most_common(5))),
                                  first_share=("first", "mean"), last_share=("last", "mean"))
        g = g.sort_values(["tokens", "n_forms"], ascending=False).head(k).reset_index()
        g["labels"] = [lab_by_stem.get(s, 0) for s in g["stem"]]
        g.insert(0, "language", lg)
        g.insert(0, "set", tag)
        rows.append(g)
    return pd.concat(rows, ignore_index=True)


def cmd_lists(args) -> None:
    tag = args.set
    out = OUT / "extras" / f"real__{tag}"
    out.mkdir(parents=True, exist_ok=True)
    data = control_data(tag) if tag in CTL_SETS else real_data(tag)[0]
    full_classes(data, tag).round(4).to_csv(out / "classes.csv", index=False)
    top_new_stems(data, tag).round(4).to_csv(out / "top_new_stems.csv", index=False)
    print(f"lists {tag}: done", flush=True)


# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #


def _ci(s: pd.Series) -> dict:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if not len(s):
        return {"mean": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
    return {"mean": float(s.mean()), "lo": float(s.quantile(0.025)), "hi": float(s.quantile(0.975)), "n": int(len(s))}


def composition_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (st, src), g in runs.groupby(["set", "source"]):
        for cmp, dname in CMPS:
            for f in COMP_F:
                col = f"{cmp}__{dname}__comp__{f}"
                if col in g.columns:
                    rows.append({"set": st, "source": src, "cmp": cmp, "direction": dname, "metric": f, **_ci(g[col])})
    return pd.DataFrame(rows)


def profile_table(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (st, src), g in runs.groupby(["set", "source"]):
        for cmp, dname in CMPS:
            for cls in PROFILE_CLASSES:
                for f in PROFILE_F:
                    a, b = f"{cmp}__{dname}__{cls}__{f}", f"{cmp}__{dname}__{cls}_m__{f}"
                    if a not in g.columns:
                        continue
                    x, y = pd.to_numeric(g[a], errors="coerce"), pd.to_numeric(g[b], errors="coerce")
                    rows.append({"set": st, "source": src, "cmp": cmp, "direction": dname, "cls": cls, "metric": f,
                                 **{f"cls_{k}": v for k, v in _ci(x).items()},
                                 **{f"matched_{k}": v for k, v in _ci(y).items()},
                                 **{f"diff_{k}": v for k, v in _ci(x - y).items()}})
    return pd.DataFrame(rows)


ED_MARKER = {"ZL": r"ed", "GC": r"c[87]"}


def new_stem_marker(classes: pd.DataFrame) -> pd.DataFrame:
    """Share of the new-stem tokens (full text, per class) whose stem / form contains the ``ed`` glyph pair (v101 ``c8`` / ``c7``)."""
    rows = []
    for (st, dname), g in classes.groupby(["set", "direction"]):
        pat = ED_MARKER.get(st.split("_")[0])
        if pat is None:
            continue
        for grp, sel in (("only", g["cls"].isin(F.ONLY_CLASSES)), ("hapax", g["cls"] == "hapax")):
            x = g[sel]
            new = x[x["stem_new"] == True]  # noqa: E712
            tot, tn = float(x["tokens"].sum()), float(new["tokens"].sum())
            has_stem = new["stem"].astype(str).str.contains(pat, regex=True)
            has_form = new["form"].astype(str).str.contains(pat, regex=True)
            rows.append({"set": st, "direction": dname, "group": grp, "tokens": tot, "new_stem_tokens": tn,
                         "new_stem_share": tn / tot if tot else np.nan,
                         "stem_has_ed": float(new.loc[has_stem, "tokens"].sum() / tn) if tn else np.nan,
                         "form_has_ed": float(new.loc[has_form, "tokens"].sum() / tn) if tn else np.nan,
                         "n_new_stems": int(new["stem"].nunique())})
    return pd.DataFrame(rows)


def cmd_summarize(args) -> None:
    runs = []
    for d in sorted((OUT / "runs").glob("*__*")):
        f = d / "runs.csv"
        if f.exists():
            runs.append(pd.read_csv(f))
    runs = pd.concat(runs, ignore_index=True)
    comp = composition_table(runs)
    comp.round(5).to_csv(OUT / "composition.csv", index=False)
    prof = profile_table(runs)
    prof.round(5).to_csv(OUT / "class_profiles.csv", index=False)
    prof[prof["source"] == "control"].round(5).to_csv(OUT / "latin_reference.csv", index=False)
    for name in ("edy_rule", "edy_contexts", "classes", "top_new_stems"):
        parts = [pd.read_csv(f) for f in sorted((OUT / "extras").glob(f"*/{name}.csv"))]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(OUT / f"{name}.csv", index=False)
    f = OUT / "classes.csv"
    if f.exists():
        new_stem_marker(pd.read_csv(f)).round(5).to_csv(OUT / "new_stems_ed.csv", index=False)
    sizes = runs.groupby(["set", "source"]).agg(n_runs=("seed", "size"), N=("N", "first"),
                                                 mean_elapsed_s=("elapsed_s", "mean"))
    sizes.to_csv(OUT / "summary.csv")
    meta = {"sources": sorted(runs["source"].unique().tolist()), "sets": sorted(runs["set"].unique().tolist()),
            "n_runs": int(len(runs)), "min_count": MIN_COUNT, "max_dist": MAX_DIST, "n_shuf": N_SHUF, "n_perm": N_PERM,
            "min_ctx": MIN_CTX, "classes": list(F.CLASSES), "profile_metrics": list(PROFILE_F), "pairs": PAIRS,
            "contexts": list(CONTEXTS), "pliny_A": PLINY_A}
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    try:
        make_plots(comp, prof)
    except Exception as e:  # matplotlib optional
        print("plots skipped:", e)
    print("summarize: done", flush=True)


def make_plots(comp: pd.DataFrame, prof: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1. composition: stacked token shares of the classes, AB y_only vs BB per set / source
    cells = []
    for st in list(SETS) + list(CTL_SETS):
        for src in ("real", "control", "ts2020", "tmpl"):
            for cmp, dname in (("AB", "y_only"), ("BB", "x_only"), ("AB", "x_only")):
                g = comp[(comp["set"] == st) & (comp["source"] == src) & (comp["cmp"] == cmp) & (comp["direction"] == dname)]
                if g.empty:
                    continue
                cells.append((f"{st}\n{src}\n{cmp}:{dname}", g.set_index("metric")["mean"]))
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(cells)), 4.5))
    bottom = np.zeros(len(cells))
    for cls, col in zip(("share_d_variant", "share_other_le1", "share_no_counterpart", "share_hapax"), ("C3", "C1", "C0", "0.6")):
        vals = np.asarray([c[1].get(cls, 0) for c in cells])
        ax.bar(range(len(cells)), vals, bottom=bottom, color=col, label=cls.replace("share_", ""))
        bottom += vals
    ax.set_xticks(range(len(cells)))
    ax.set_xticklabels([c[0] for c in cells], rotation=90, fontsize=6)
    ax.set_ylabel("share of X tokens")
    ax.set_title("Forms only in X (B-only: AB y_only; halves: BB; A-only: AB x_only)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT / "composition.png", dpi=130)
    plt.close(fig)

    # 2. class profiles: class - matched shared, real B-only vs halves of B vs Latin
    metrics = ("paradigm_forms", "H_suffix_given_stem", "pages_per_token", "gries_dp", "section_entropy",
               "label_share", "first_share", "last_share", "par_first_share", "cooc_ratio", "bigram_ll", "seam_jsd")
    fig, axes = plt.subplots(3, 4, figsize=(14, 9))
    for ax, m in zip(axes.ravel(), metrics):
        k, labels = 0, []
        for st, src in (("ZL_split", "real"), ("GC_split", "real"), ("pliny", "control"), ("caesar_celsus", "control"),
                        ("ZL_split", "ts2020"), ("ZL_split", "tmpl")):
            for cmp, dname, col in (("AB", "y_only", "C3"), ("BB", "x_only", "C0")):
                for cls, mk in (("d_variant", "s"), ("other_le1", "^"), ("no_counterpart", "o"), ("hapax", "x")):
                    g = prof[(prof["set"] == st) & (prof["source"] == src) & (prof["cmp"] == cmp) & (prof["direction"] == dname)
                             & (prof["cls"] == cls) & (prof["metric"] == m)]
                    if g.empty or not np.isfinite(g["diff_mean"].iloc[0]):
                        continue
                    r = g.iloc[0]
                    ax.errorbar(k, r["diff_mean"], yerr=[[max(r["diff_mean"] - r["diff_lo"], 0)], [max(r["diff_hi"] - r["diff_mean"], 0)]],
                                fmt=mk, color=col, ms=4, capsize=2)
                    labels.append(f"{st[:2]}{src[0]}:{cmp}:{cls[:4]}")
                    k += 1
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=4)
        ax.set_title(m, fontsize=9)
        ax.grid(alpha=0.3)
    fig.suptitle("Class minus frequency-matched shared forms (red: A×B, blue: halves of B)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "class_profiles.png", dpi=130)
    plt.close(fig)

    # 3. -edy rule: P(second form | context value) per language, from edy_contexts.csv
    f = OUT / "edy_contexts.csv"
    if f.exists():
        ctxs = pd.read_csv(f)
        pairs = [p for p in ("ey_edy", "c9_c89", "y_dy", "ol_or") if p in set(ctxs["pair"])]
        fig, axes = plt.subplots(1, len(pairs), figsize=(4 * len(pairs), 3.6), squeeze=False)
        for ax, pair in zip(axes[0], pairs):
            for lg, col in (("A", "C0"), ("B", "C3")):
                g = ctxs[(ctxs["pair"] == pair) & (ctxs["language"] == lg) & (ctxs["source"] == "real")
                         & (ctxs["context"] == "prev_unit") & (ctxs["n"] >= MIN_CTX)]
                if g.empty:
                    continue
                ax.hist(g["p_s2"], bins=np.linspace(0, 1, 11), weights=g["n"], alpha=0.5, color=col, label=f"{lg} (n={int(g['n'].sum())})")
            ax.set_title(f"{pair}: P(2nd form | unit before), token-weighted", fontsize=8)
            ax.set_xlabel("P(second form)")
            ax.legend(fontsize=7)
            ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "edy_rule.png", dpi=130)
        plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("measure", "edy"):
        p = sub.add_parser(name)
        p.add_argument("--set", required=True, choices=SETS + CTL_SETS)
        p.add_argument("--source", default="real", choices=("real", "ts2020", "tmpl"))
        if name == "measure":
            p.add_argument("--seeds", default=None)
    p = sub.add_parser("lists")
    p.add_argument("--set", required=True, choices=SETS + CTL_SETS)
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    {"measure": cmd_measure, "edy": cmd_edy, "lists": cmd_lists, "summarize": cmd_summarize}[args.cmd](args)


if __name__ == "__main__":
    main()
