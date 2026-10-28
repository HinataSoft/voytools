"""Experiment 016: are the stems found only in Currier B (only in A) glyph variants of the other language's stems?

For two token sets X, Y at equal size (whole pages drawn to N tokens, N = half the smaller language of
the scope, as in 012) and four representations of a token - stem of decomposition (a), stem of (b),
whole token, whole token in BPE units (64 merges) - the script measures (``voytools.registers``):

1. counterparts: X-only stems (>= 2 tokens in X, absent from Y) and the share of their tokens whose
   nearest Y stem (Levenshtein on glyph units, exact up to 2 edits) is within 1 / 2 edits;
2. edit inventory: the operations of the counterpart pairs (``ins:ed@end`` ...), ranked, and how many
   of them explain 50 / 80 / 95 % of the tokens;
3. substitution: X-only stems whose script uses only the top-k operations (k = 1, 3, 5, 10, all) are
   rewritten to their counterpart; Jaccard / coverage of the rewritten X against Y;
4. the same on whole tokens and BPE units (no parser);
5. section check (full text): JSD of the section profiles of counterpart pairs vs random re-pairings
   and vs the same stem in A and B.

Comparisons in every draw: A vs B (``AB``), two halves of A (``AA``), two halves of B (``BB``), and
B-only stems against the shuffled glyphs of the A stems (``SHUF``, random-form null).  Both directions:
``x_only`` = X-only stems vs Y (A-only vs B for ``AB``), ``y_only`` = Y-only vs X (B-only vs A).  Edit
operations are always written in the X -> Y (A -> B) direction.

Sources: real (ZL / GC, split / join; scopes ``all`` and ``herbal``; 200 draws), ts2020 config 21 and the
005 template (A/B = layout labels; 30 seeds, one draw each), Latin controls: Pliny books 16-21 vs 23-30
("same language, other topic") and Caesar vs Celsus ("same language, other register"), 200 draws.

Usage (from the repository root)::

    python experiments/016_ab_stem_variants/run.py measure --set ZL_split                  # 200 draws
    python experiments/016_ab_stem_variants/run.py measure --set ZL_split --source ts2020  # 30 seeds
    python experiments/016_ab_stem_variants/run.py measure --set pliny                     # control
    python experiments/016_ab_stem_variants/run.py extras --set ZL_split                   # counterparts, inventory, sections
    python experiments/016_ab_stem_variants/run.py summarize
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

from voytools import registers as G  # noqa: E402
from voytools import run_seeds  # noqa: E402
from voytools._common import line_ids  # noqa: E402
from voytools.glyphs import collapse_eva  # noqa: E402
from voytools.morphology import apply_bpe, decompose_table, learn_bpe  # noqa: E402
from voytools.nulls import as_table  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M12 = _load("exp012", "012_ab_register")
M11 = M12.M11
M9 = M12.M9

OUT = ROOT / "results" / "016_ab_stem_variants"
CACHE = ROOT / "scratch" / "016_cache"
SETS = ("ZL_split", "GC_split", "ZL_join", "GC_join")
CTL_SETS = ("pliny", "caesar_celsus")
NULL_SETS = {"ts2020": ("ZL_split", "ZL_join"), "tmpl": ("ZL_split", "GC_split")}
SCOPES = ("all", "herbal")
REPS = ("a", "b", "tok", "bpe")
N_DRAW = 200
N_NULL = 30
N_INV = 20            # draws pooled for the edit inventory
MIN_COUNT = 2         # stems with >= 2 tokens (inventories, X-only, candidates)
MAX_DIST = 2
KS = (0, 1, 3, 5, 10, "all")
BPE_MERGES = 64
N_PERM = 500          # section check permutations
PLINY_A = ("16", "18", "21")
RATE_F = ("share_only", "n_only", "n_cand", "rate_le1", "rate_le2", "type_le1", "type_le2")
# Eva glyphs of the Currier / R&T markers: ed, e, qo, dy, ch/sh
CURRIER_UNITS = {"e", "ee", "d", "q"}


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


def _tuples(strings, preset) -> np.ndarray:
    cache: dict = {}
    out = np.empty(len(strings), dtype=object)
    for i, s in enumerate(strings):
        t = cache.get(s)
        if t is None:
            t = cache[s] = tuple(collapse_eva(s, preset)) if s else ()
        out[i] = t
    return out


def _object_array(items) -> np.ndarray:
    out = np.empty(len(items), dtype=object)
    for i, t in enumerate(items):
        out[i] = tuple(t)
    return out


class Data:
    """One corpus: page metadata, unit tuples of every token for each representation, scope pools."""

    def __init__(self, df: pd.DataFrame, ctx, rank: dict, reps=REPS, merges=None):
        self.df = df.reset_index(drop=True)
        self.page = self.df["folio"].to_numpy(dtype=object)
        meta = self.df.groupby("folio", sort=False).agg(lang=("language", "first"), hand=("hand", "first"),
                                                        sec=("illustration", "first"), n=("token", "size"))
        meta["order"] = [rank.get(p, i) for i, p in enumerate(meta.index)]
        self.meta = meta.sort_values("order")
        self.preset = ctx.gpreset
        self.units: dict = {}
        self.strings: dict = {}
        for var, (pre, suf) in M11.variants(ctx).items():
            if var in reps:
                d = decompose_table(self.df, prefixes=pre, suffixes=suf)
                self.strings[var] = d["stem"].to_numpy(dtype=object)
                self.units[var] = _tuples(self.strings[var], self.preset)
        tokens = self.df["token"].tolist()
        if "tok" in reps:
            self.strings["tok"] = np.asarray(tokens, dtype=object)
            self.units["tok"] = _tuples(tokens, self.preset)
        if "bpe" in reps:
            self.merges = merges if merges is not None else learn_bpe(tokens, BPE_MERGES, preset=self.preset)
            self.strings["bpe"] = np.asarray(tokens, dtype=object)
            self.units["bpe"] = _object_array(apply_bpe(tokens, self.merges, preset=self.preset))
        self.reps = [r for r in reps if r in self.units]

    def pool(self, scope: str, lang: str) -> list:
        m = self.meta
        sel = m["lang"] == lang
        if scope == "herbal":
            sel &= m["sec"] == "H"
        return m.index[sel].tolist()

    def size(self, pages) -> int:
        return int(self.meta.loc[list(pages), "n"].sum())

    def counts(self, rep: str, mask: np.ndarray) -> dict:
        c = Counter(self.units[rep][mask].tolist())
        c.pop((), None)
        return dict(c)


def scope_n(data: Data, scope: str) -> int:
    return min(data.size(data.pool(scope, "A")), data.size(data.pool(scope, "B"))) // 2


# --------------------------------------------------------------------------- #
# measures
# --------------------------------------------------------------------------- #


def one_direction(cx: dict, cy: dict, x_is_query: bool, rng=None, shuffled_cands: bool = False) -> dict:
    """Counterpart rates (+ inventory and substitution) of the X-only stems vs Y (or Y-only vs X).

    Operations are written X -> Y in both cases.  With ``shuffled_cands`` the candidates are the shuffled
    glyphs of the other side (rates only).
    """
    cq, cc = (cx, cy) if x_is_query else (cy, cx)
    if shuffled_cands:
        types = [s for s, c in cc.items() if c >= MIN_COUNT]
        shuf = G.shuffle_units(types, rng)
        cands: dict = {}
        for s, t in zip(types, shuf):
            cands[t] = cands.get(t, 0) + cc[s]
        r = G.counterpart_rates(cq, cands, MIN_COUNT, MAX_DIST, count_y_only=cc)
        return {f: r[f] for f in RATE_F}
    r = G.counterpart_rates(cq, cc, MIN_COUNT, MAX_DIST)
    out = {f: r[f] for f in RATE_F}
    scripts = G.edit_scripts(r["pairs"], "xy" if x_is_query else "yx")
    inv = G.edit_inventory(scripts, cq)
    cov = G.coverage_levels(inv)
    out.update({"n_ops": len(inv), "ops_50": cov[0.5], "ops_80": cov[0.8], "ops_95": cov[0.95],
                "top1_share": float(inv["share"].iloc[0]) if len(inv) else float("nan")})
    sub = G.substitution_sharing(cq, cc, r["pairs"], scripts, inv["op"].tolist(), KS, MIN_COUNT)
    for k, v in sub.items():
        out[f"jac_{k}"] = v["jaccard"]
        out[f"cov_{k}"] = v["cover_xy"]
    out["nrw_all"] = sub["all"]["n_rewritten"]
    out["_pairs"], out["_scripts"], out["_weights"] = r["pairs"], scripts, cq
    return out


def compare(data: Data, mx: np.ndarray, my: np.ndarray, prefix: str, rng=None, shuf: bool = False,
            keep: dict | None = None) -> dict:
    out = {}
    for rep in data.reps:
        cx, cy = data.counts(rep, mx), data.counts(rep, my)
        for dname, xq in (("x_only", True), ("y_only", False)):
            r = one_direction(cx, cy, xq)
            if keep is not None:
                keep[(prefix, rep, dname)] = (r.pop("_pairs"), r.pop("_scripts"), r.pop("_weights"))
            else:
                for k in ("_pairs", "_scripts", "_weights"):
                    r.pop(k)
            out.update({f"{prefix}__{rep}__{dname}__{k}": v for k, v in r.items()})
            if shuf:
                s = one_direction(cx, cy, xq, rng, shuffled_cands=True)
                out.update({f"{prefix.replace('AB', 'SHUF')}__{rep}__{dname}__{k}": v for k, v in s.items()})
    return out


def draw_measures(data: Data, seed: int, scopes=SCOPES, keep: dict | None = None) -> dict:
    """One equal-size draw per scope: A vs B (+ shuffled-form null), two halves of A, two halves of B."""
    rng = np.random.default_rng(seed)
    out = {}
    for scope in scopes:
        n = scope_n(data, scope)
        out[f"{scope}__N"] = n
        pa, pb = data.pool(scope, "A"), data.pool(scope, "B")
        ma, mb = G.draw_tokens(data.page, pa, n, rng), G.draw_tokens(data.page, pb, n, rng)
        out.update(compare(data, ma, mb, f"{scope}__AB", rng, shuf=True, keep=keep))
        for lab, pool, first in (("AA", pa, ma), ("BB", pb, mb)):
            second = G.draw_tokens(data.page, pool, n, rng, exclude=np.unique(data.page[first]))
            out.update(compare(data, first, second, f"{scope}__{lab}", keep=keep))
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


def control_data(name: str) -> Data:
    """Latin controls as two 'languages': Pliny books (16-21 = A, 23-30 = B) or Caesar (A) vs Celsus (B)."""
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
    rank = {p: i for i, p in enumerate(pd.unique(df["folio"].to_numpy()))}
    return Data(df, ctx, rank, reps=("a", "tok", "bpe"))


def real_data(tag: str) -> tuple:
    ctx = M9.Ctx(tag)
    return Data(ctx.real, ctx, ctx.rank), ctx


def cmd_measure(args) -> None:
    tag, source = args.set, args.source
    if tag in CTL_SETS:
        data = control_data(tag)
        scopes = ("all",)

        def one(seed: int) -> dict:
            return {"set": tag, "source": "control", "kind": "draw", **draw_measures(data, seed, scopes)}
        default = list(range(N_DRAW))
        extra = {"pliny_A": PLINY_A}
    elif source == "real":
        data, ctx = real_data(tag)
        scopes = SCOPES

        def one(seed: int) -> dict:
            return {"set": tag, "source": source, "kind": "draw", **draw_measures(data, seed, scopes)}
        default = list(range(N_DRAW))
        extra = {"bpe_merges": [list(m) for m in data.merges]}
    else:
        ctx = M9.Ctx(f"real:{tag}")
        scopes = SCOPES

        def one(seed: int) -> dict:
            data = Data(cached_corpus(ctx, source, tag, seed), ctx, ctx.rank)
            return {"set": tag, "source": source, "kind": "null", **draw_measures(data, seed, scopes)}
        default = list(range(N_NULL))
        extra = {"ts2020_config": M9.TS2020_CONFIG, "tmpl_seed0": M9.TMPL_SEED0}
    seeds = M9.parse_seeds(args.seeds, default)
    params = {"set": tag, "source": source, "scopes": scopes, "n_draw": N_DRAW, "n_null": N_NULL,
              "min_count": MIN_COUNT, "max_dist": MAX_DIST, "ks": [str(k) for k in KS], "bpe_merges_n": BPE_MERGES,
              **extra}
    run_seeds("016_ab_stem_variants", seeds, one, OUT / "runs" / f"{source}__{tag}", params=params,
              progress=M11.progress)


# --------------------------------------------------------------------------- #
# extras: counterpart table, pooled edit inventory, section check
# --------------------------------------------------------------------------- #


def _join(units) -> str:
    return "".join(units)


def _units_str(units) -> str:
    return " ".join(units)


def counterpart_table(data: Data, tag: str) -> pd.DataFrame:
    """Full A vs full B (unequal sizes): every language-specific stem, its counterpart, distance, operations."""
    rows = []
    pa, pb = data.pool("all", "A"), data.pool("all", "B")
    ma, mb = np.isin(data.page, pa), np.isin(data.page, pb)
    for rep in data.reps:
        ca, cb = data.counts(rep, ma), data.counts(rep, mb)
        for dname, cq, cc, direction in (("A_only", ca, cb, "xy"), ("B_only", cb, ca, "yx")):
            r = G.counterpart_rates(cq, cc, MIN_COUNT, MAX_DIST)
            scripts = G.edit_scripts(r["pairs"], direction)
            only = {s: c for s, c in cq.items() if c >= MIN_COUNT and not cc.get(s)}
            for s, c in sorted(only.items(), key=lambda kv: (-kv[1], kv[0])):
                p = r["pairs"].get(s)
                rows.append({"set": tag, "rep": rep, "direction": dname, "stem": _join(s), "units": _units_str(s),
                             "tokens": c, "counterpart": _join(p[0]) if p else "", "distance": p[1] if p else np.nan,
                             "counterpart_tokens": cc[p[0]] if p else np.nan,
                             "ops_A_to_B": ";".join(scripts[s]) if p else ""})
    return pd.DataFrame(rows)


def pooled_inventory(data: Data, tag: str, scopes, seed: int = 101) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Edit inventory pooled over N_INV equal-size draws for AB, AA, BB; coverage summary per cell."""
    keep: dict = {}
    for s in range(N_INV):
        local: dict = {}
        draw_measures(data, seed + s, scopes, keep=local)
        for key, (pairs, scripts, weights) in local.items():
            acc = keep.setdefault(key, ({}, {}))
            for x, ops in scripts.items():
                acc[0][(s, x)] = ops
                acc[1][(s, x)] = weights.get(x, 0)
    inv_rows, cov_rows = [], []
    for (prefix, rep, dname), (scripts, weights) in keep.items():
        scope, cmp = prefix.split("__")
        inv = G.edit_inventory(scripts, weights)
        cov = G.coverage_levels(inv)
        inv.insert(0, "direction", dname)
        inv.insert(0, "rep", rep)
        inv.insert(0, "cmp", cmp)
        inv.insert(0, "scope", scope)
        inv.insert(0, "set", tag)
        inv["currier"] = [is_currier_op(op) for op in inv["op"]] if data.preset else np.nan
        inv_rows.append(inv.head(60))
        cur = float(inv.loc[inv["currier"] == True, "tokens"].sum() / inv["tokens"].sum()) \
            if data.preset and len(inv) and inv["tokens"].sum() else np.nan  # noqa: E712
        cov_rows.append({"set": tag, "scope": scope, "cmp": cmp, "rep": rep, "direction": dname,
                         "n_pairs": len(scripts), "tokens": float(sum(weights.values())), "n_ops": len(inv),
                         "ops_50": cov[0.5], "ops_80": cov[0.8], "ops_95": cov[0.95],
                         "top1_op": inv["op"].iloc[0] if len(inv) else "", "top1_share": inv["share"].iloc[0] if len(inv) else np.nan,
                         "top3_ops": ";".join(inv["op"].head(3)), "currier_token_share": cur})
    return pd.concat(inv_rows, ignore_index=True), pd.DataFrame(cov_rows)


def is_currier_op(op: str) -> bool:
    """Insertion / deletion made only of e, ee, d, q (ed, eed, dy-, qo- ...) or a ch <-> sh substitution."""
    kind, rest = op.split(":", 1)
    units, _ = rest.rsplit("@", 1)
    if kind == "sub":
        return set(units.split(">")) == {"ch", "sh"}
    return all(u in CURRIER_UNITS for u in collapse_eva(units, "extended"))


def section_check(data: Data, tag: str, reps=("a", "tok")) -> pd.DataFrame:
    """JSD of the section profiles of counterpart pairs vs permuted pairs vs shared stems (same stem in A and B).

    Two profiles: ``sections`` = distribution over all illustration codes (A and B have different
    section inventories, so this JSD is dominated by the language) and ``herbal_binary`` = share of
    the stem's tokens in the herbal section vs the rest (the only section with real mass in both).
    """
    secs = sorted(data.meta["sec"].astype(str).unique())
    sec_tok = data.df["illustration"].astype(str).to_numpy()
    pa, pb = data.pool("all", "A"), data.pool("all", "B")
    ma, mb = np.isin(data.page, pa), np.isin(data.page, pb)
    rows = []
    for rep in reps:
        units = data.units[rep]

        def profiles(mask, cats, key):
            prof: dict = {}
            idx = {s: i for i, s in enumerate(cats)}
            for u, sec in zip(units[mask].tolist(), sec_tok[mask]):
                if u == ():
                    continue
                v = prof.get(u)
                if v is None:
                    v = prof[u] = np.zeros(len(cats))
                v[idx[key(sec)]] += 1
            return prof
        ca, cb = data.counts(rep, ma), data.counts(rep, mb)
        for profile, cats, key in (("sections", secs, lambda s: s),
                                   ("herbal_binary", ["H", "other"], lambda s: "H" if s == "H" else "other")):
            pA, pB = profiles(ma, cats, key), profiles(mb, cats, key)
            for dname, cq, cc, pq, pc in (("A_only", ca, cb, pA, pB), ("B_only", cb, ca, pB, pA)):
                r = G.counterpart_rates(cq, cc, MIN_COUNT, MAX_DIST)
                for dmax in (1, 2):
                    pairs = {x: y for x, (y, d) in r["pairs"].items() if d <= dmax}
                    t = G.profile_test(pairs, pq, pc, N_PERM, seed=1)
                    rows.append({"set": tag, "rep": rep, "profile": profile, "direction": dname,
                                 "kind": f"counterparts_le{dmax}", **t})
            shared = {s: s for s in ca if ca[s] >= MIN_COUNT and cb.get(s, 0) >= MIN_COUNT}
            t = G.profile_test(shared, pA, pB, N_PERM, seed=1)
            rows.append({"set": tag, "rep": rep, "profile": profile, "direction": "shared",
                         "kind": "same_stem_A_vs_B", **t})
    return pd.DataFrame(rows)


def cmd_extras(args) -> None:
    tag = args.set
    if tag in CTL_SETS:
        data = control_data(tag)
        scopes = ("all",)
    else:
        data, _ = real_data(tag)
        scopes = SCOPES
    out = OUT / "extras" / tag
    out.mkdir(parents=True, exist_ok=True)
    counterpart_table(data, tag).to_csv(out / "counterparts.csv", index=False)
    print(f"extras {tag}: counterparts", flush=True)
    inv, cov = pooled_inventory(data, tag, scopes)
    inv.round(6).to_csv(out / "edit_inventory.csv", index=False)
    cov.round(6).to_csv(out / "edit_coverage.csv", index=False)
    print(f"extras {tag}: inventory", flush=True)
    if tag not in CTL_SETS:
        section_check(data, tag).round(6).to_csv(out / "section_check.csv", index=False)
    print(f"extras {tag}: done", flush=True)


# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #


def _ci(s: pd.Series) -> dict:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if not len(s):
        return {"mean": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
    return {"mean": float(s.mean()), "lo": float(s.quantile(0.025)), "hi": float(s.quantile(0.975)), "n": int(len(s))}


def _long(runs: pd.DataFrame, metrics) -> pd.DataFrame:
    rows = []
    for (st, src), g in runs.groupby(["set", "source"]):
        for col in g.columns:
            parts = col.split("__")
            if len(parts) != 5 or parts[4] not in metrics:
                continue
            scope, cmp, rep, dname, metric = parts
            rows.append({"set": st, "source": src, "scope": scope, "cmp": cmp, "rep": rep, "direction": dname,
                         "metric": metric, **_ci(g[col])})
    return pd.DataFrame(rows)


def _diff_rows(runs: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Per-draw AB minus mean of the halves (AA, BB pooled over both directions), summarised."""
    rows = []
    for (st, src), g in runs.groupby(["set", "source"]):
        scopes = sorted({c.split("__")[0] for c in g.columns if c.endswith(f"__{metric}")})
        reps = sorted({c.split("__")[2] for c in g.columns if c.endswith(f"__{metric}")})
        for scope in scopes:
            for rep in reps:
                for dname in ("x_only", "y_only"):
                    ab = g.get(f"{scope}__AB__{rep}__{dname}__{metric}")
                    if ab is None:
                        continue
                    halves = pd.concat([g[f"{scope}__{c}__{rep}__{d}__{metric}"] for c in ("AA", "BB")
                                        for d in ("x_only", "y_only")], axis=1).mean(axis=1)
                    d = pd.to_numeric(ab, errors="coerce") - halves
                    rows.append({"set": st, "source": src, "scope": scope, "rep": rep, "direction": dname,
                                 "metric": metric, **{f"AB_{k}": v for k, v in _ci(ab).items()},
                                 **{f"halves_{k}": v for k, v in _ci(halves).items()},
                                 **{f"diff_{k}": v for k, v in _ci(d).items()}})
    return pd.DataFrame(rows)


def cmd_summarize(args) -> None:
    runs = []
    for d in sorted((OUT / "runs").glob("*__*")):
        f = d / "runs.csv"
        if f.exists():
            runs.append(pd.read_csv(f))
    runs = pd.concat(runs, ignore_index=True)
    rates = _long(runs, ("share_only", "rate_le1", "rate_le2", "type_le1", "type_le2", "n_only", "n_cand",
                         "ops_50", "ops_80", "ops_95", "n_ops", "top1_share"))
    rates.round(5).to_csv(OUT / "counterpart_rates.csv", index=False)
    diffs = pd.concat([_diff_rows(runs, m) for m in ("rate_le1", "rate_le2", "share_only")], ignore_index=True)
    diffs.round(5).to_csv(OUT / "counterpart_excess.csv", index=False)
    sub = _long(runs, tuple(f"jac_{k}" for k in KS) + tuple(f"cov_{k}" for k in KS) + ("nrw_all",))
    sub.round(5).to_csv(OUT / "substitution_sharing.csv", index=False)
    whole = pd.concat([rates[rates["rep"].isin(["tok", "bpe"])], sub[sub["rep"].isin(["tok", "bpe"])]],
                      ignore_index=True)
    whole.round(5).to_csv(OUT / "whole_tokens.csv", index=False)
    for name in ("counterparts", "edit_inventory", "edit_coverage", "section_check"):
        parts = [pd.read_csv(f) for f in sorted((OUT / "extras").glob(f"*/{name}.csv"))]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(OUT / f"{name}.csv", index=False)
    sizes = runs.groupby(["set", "source"]).agg(n_runs=("seed", "size"), N_all=("all__N", "first"))
    if "herbal__N" in runs.columns:
        sizes["N_herbal"] = runs.groupby(["set", "source"])["herbal__N"].first()
    sizes.to_csv(OUT / "summary.csv")
    meta = {"sources": sorted(runs["source"].unique().tolist()), "sets": sorted(runs["set"].unique().tolist()),
            "n_runs": int(len(runs)), "min_count": MIN_COUNT, "max_dist": MAX_DIST, "ks": [str(k) for k in KS],
            "bpe_merges": BPE_MERGES, "n_inv_draws": N_INV, "n_perm_sections": N_PERM, "pliny_A": PLINY_A}
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    try:
        make_plots(rates, diffs, sub)
    except Exception as e:  # matplotlib optional
        print("plots skipped:", e)
    print("summarize: done", flush=True)


def make_plots(rates: pd.DataFrame, diffs: pd.DataFrame, sub: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1. counterpart rates <= 1 / <= 2: AB vs halves vs shuffled, per set and representation (scope all)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, metric in zip(axes, ("rate_le1", "rate_le2")):
        d = rates[(rates["metric"] == metric) & (rates["scope"] == "all") & (rates["direction"] == "y_only")]
        labels, k = [], 0
        for st in list(SETS) + list(CTL_SETS):
            for src in ("real", "control", "ts2020", "tmpl"):
                for rep in ("a", "tok", "bpe"):
                    g = d[(d["set"] == st) & (d["source"] == src) & (d["rep"] == rep)]
                    if g.empty:
                        continue
                    for cmp, col, dx in (("AB", "C3", -0.2), ("BB", "C0", 0.0), ("SHUF", "0.6", 0.2)):
                        r = g[g["cmp"] == cmp]
                        if r.empty:
                            continue
                        r = r.iloc[0]
                        ax.errorbar(k + dx, r["mean"], yerr=[[r["mean"] - r["lo"]], [r["hi"] - r["mean"]]], fmt="o",
                                    color=col, ms=4, capsize=2, label=cmp if k == 0 else None)
                    labels.append(f"{st}\n{src}:{rep}")
                    k += 1
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_title(f"B-only tokens with a counterpart ({metric})")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[0].set_ylabel("share of B-only tokens")
    fig.tight_layout()
    fig.savefig(OUT / "counterpart_rates.png", dpi=130)
    plt.close(fig)

    # 2. sharing after substitution vs k
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, rep in zip(axes, ("a", "tok")):
        for st, col in zip(("ZL_split", "GC_split", "pliny", "caesar_celsus"), ("C3", "C1", "C2", "C4")):
            for cmp, ls in (("AB", "-"), ("AA", "--"), ("BB", ":")):
                ys, los, his = [], [], []
                for k in KS:
                    g = sub[(sub["set"] == st) & (sub["source"].isin(["real", "control"])) & (sub["scope"] == "all")
                            & (sub["rep"] == rep) & (sub["direction"] == "y_only") & (sub["cmp"] == cmp)
                            & (sub["metric"] == f"jac_{k}")]
                    ys.append(g["mean"].iloc[0] if len(g) else np.nan)
                    los.append(g["lo"].iloc[0] if len(g) else np.nan)
                    his.append(g["hi"].iloc[0] if len(g) else np.nan)
                x = range(len(KS))
                ax.plot(x, ys, ls, color=col, marker="o", ms=3, label=f"{st} {cmp}" if cmp != "BB" else None)
                if cmp == "AB":
                    ax.fill_between(x, los, his, color=col, alpha=0.12)
        ax.set_xticks(range(len(KS)))
        ax.set_xticklabels([str(k) for k in KS])
        ax.set_xlabel("k top operations allowed (B-only stems rewritten to A)")
        ax.set_title(f"Jaccard after substitution, rep {rep}")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Jaccard (stems >= 2 tokens)")
    axes[0].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "substitution_sharing.png", dpi=130)
    plt.close(fig)

    # 3. excess of AB over halves
    fig, ax = plt.subplots(figsize=(9, 4))
    d = diffs[(diffs["metric"] == "rate_le1") & (diffs["scope"] == "all") & (diffs["direction"] == "y_only")]
    d = d.sort_values(["set", "source", "rep"]).reset_index(drop=True)
    ax.axhline(0, color="k", lw=0.8)
    ax.errorbar(range(len(d)), d["diff_mean"], yerr=[d["diff_mean"] - d["diff_lo"], d["diff_hi"] - d["diff_mean"]],
                fmt="o", ms=4, capsize=2, color="C3")
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels([f"{r.set}\n{r.source}:{r.rep}" for r in d.itertuples()], rotation=90, fontsize=6)
    ax.set_ylabel("rate_le1: A×B − halves")
    ax.set_title("Counterparts within 1 edit: excess of A×B over halves of one language")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "counterpart_excess.png", dpi=130)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("measure")
    p.add_argument("--set", required=True, choices=SETS + CTL_SETS)
    p.add_argument("--source", default="real", choices=("real", "ts2020", "tmpl"))
    p.add_argument("--seeds", default=None)
    p = sub.add_parser("extras")
    p.add_argument("--set", required=True, choices=SETS + CTL_SETS)
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    {"measure": cmd_measure, "extras": cmd_extras, "summarize": cmd_summarize}[args.cmd](args)


if __name__ == "__main__":
    main()
