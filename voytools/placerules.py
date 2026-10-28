"""Place rules: how many concrete items carry the positional and word-seam information (experiment 014).

Three groups of tools:

1. **Item decomposition of mutual information.**  The plug-in MI between two labels
   I(X; Y) = sum_{x,y} p(x,y) log2 p(x,y) / (p(x) p(y)) is split into pointwise terms per pair
   (:func:`pointwise_contributions`), summed per value of X or Y (:func:`marginal_contributions`;
   these sums equal p(x) KL(p(y|x) || p(y)) and are never negative), and the number of items needed to
   cover 50 / 80 / 95 % of the total is read off (:func:`coverage`).  :func:`item_decomposition` does
   the three steps at once.  :func:`position_chain` splits I(Y; position class) along the chain
   place in line -> paragraph-first line -> relative position inside the line, so that "edge of the
   line" and "drift inside the line" (Feaster) get separate numbers; :func:`seam_pairs` gives the
   (label of token t, label of token t+1) pairs inside lines and :func:`seam_glyph_vs_affix` tests
   whether anything remains of the suffix -> next prefix dependence once the edge glyphs are known.

2. **Stability of a rule table across groups.**  :func:`conditional_jsd` compares two conditional
   tables p_A(y | x) and p_B(y | x) (weighted Jensen-Shannon divergence); :func:`stability` compares
   the divergence between two groups (A/B, scribes, sections) with the divergence between random
   halves of the same group (lines as units) and counts the "rules" (pointwise terms covering 80 % of
   the pooled information) whose direction agrees in both groups.

3. **Post-processor for a generated corpus.**  :func:`learn_place_rules` reads the place tables of
   a real token table and :func:`apply_place_rules` imposes them on another corpus (for example a
   ts2020 sample of :func:`voytools.selfcitation.selfcitation_like`): R1 the prefix distribution of
   the first token of a line and of the tokens of a paragraph-first line, R2 the suffix distribution
   of the last token of a line, R3 the prefix of every inner token conditional on the last glyph of
   the token before.  Prefixes and suffixes are rewritten with the *minimal* coupling that turns the
   generator's distribution in the context into the target distribution: a token keeps its affix
   with probability min(1, q(a) / p(a)) and otherwise draws a replacement from the excess
   (q - p)_+ ; stems (and the other affix) are kept.  :mod:`voytools.nulls` is not touched.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from ._common import glyph_lists, line_ids, resolve_preset
from .glyphs import collapse_eva
from .morphology import decompose_table
from .stats import jsd
from .variation import line_place, position_class

__all__ = [
    "pointwise_contributions",
    "marginal_contributions",
    "coverage",
    "item_decomposition",
    "mutual_information",
    "position_chain",
    "token_labels",
    "seam_pairs",
    "seam_glyph_vs_affix",
    "conditional_jsd",
    "stability",
    "PlaceRules",
    "learn_place_rules",
    "apply_place_rules",
]

LEVELS = (0.5, 0.8, 0.95)


# --------------------------------------------------------------------------- #
# codes and plug-in information
# --------------------------------------------------------------------------- #


def _codes(values) -> tuple[np.ndarray, list]:
    c, u = pd.factorize(pd.Series(list(values), dtype=object))
    return c.astype(np.int64), list(u)


def _joint(*arrays: np.ndarray) -> np.ndarray:
    key = arrays[0].astype(np.int64)
    for a in arrays[1:]:
        key = key * (int(a.max()) + 1 if len(a) else 1) + a
    return pd.factorize(key)[0].astype(np.int64)


def _H(c: np.ndarray) -> float:
    if len(c) == 0:
        return float("nan")
    p = np.bincount(c) / len(c)
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def mutual_information(x, y, z=None) -> float:
    """Plug-in I(X; Y | Z) in bits between label sequences (``z=None``: unconditional)."""
    xc, _ = _codes(x)
    yc, _ = _codes(y)
    if z is None:
        return _H(xc) + _H(yc) - _H(_joint(xc, yc))
    zc, _ = _codes(z)
    return _H(_joint(xc, zc)) + _H(_joint(yc, zc)) - _H(_joint(xc, yc, zc)) - _H(zc)


# --------------------------------------------------------------------------- #
# item decomposition
# --------------------------------------------------------------------------- #


def pointwise_contributions(x, y, x_name: str = "x", y_name: str = "y") -> pd.DataFrame:
    """Pointwise terms p(x,y) log2 p(x,y)/(p(x)p(y)) of the plug-in MI, one row per observed pair.

    Columns: ``x_name``, ``y_name``, ``n``, ``p_xy``, ``p_x``, ``p_y``, ``log_ratio``
    (log2 p(x,y)/(p(x)p(y)); > 0 = the pair is enriched, < 0 = depleted) and ``contribution``
    (bits; the column sums to I(X; Y)).  Sorted by ``contribution`` descending.
    """
    xs = pd.Series(list(x), dtype=object)
    ys = pd.Series(list(y), dtype=object)
    n = len(xs)
    if n == 0:
        return pd.DataFrame(columns=[x_name, y_name, "n", "p_xy", "p_x", "p_y", "log_ratio", "contribution"])
    pairs = pd.DataFrame({x_name: xs, y_name: ys}).groupby([x_name, y_name], sort=False).size().rename("n").reset_index()
    nx = xs.value_counts()
    ny = ys.value_counts()
    pairs["p_xy"] = pairs["n"] / n
    pairs["p_x"] = pairs[x_name].map(nx).to_numpy() / n
    pairs["p_y"] = pairs[y_name].map(ny).to_numpy() / n
    pairs["log_ratio"] = np.log2(pairs["p_xy"] / (pairs["p_x"] * pairs["p_y"]))
    pairs["contribution"] = pairs["p_xy"] * pairs["log_ratio"]
    return pairs.sort_values("contribution", ascending=False, kind="stable").reset_index(drop=True)


def marginal_contributions(pairs: pd.DataFrame, by: str) -> pd.DataFrame:
    """Sum of the pointwise terms per value of ``by`` (= p(x) KL(p(y|x) || p(y)), never negative).

    Columns: ``by``, ``n``, ``contribution``, ``share`` (of the total), ``rank``, ``cum_share``.
    """
    g = pairs.groupby(by, sort=False).agg(n=("n", "sum"), contribution=("contribution", "sum")).reset_index()
    g["contribution"] = g["contribution"].clip(lower=0.0)  # floating noise only; the sum is a KL divergence
    g = g.sort_values("contribution", ascending=False, kind="stable").reset_index(drop=True)
    total = float(g["contribution"].sum())
    g["share"] = g["contribution"] / total if total > 0 else np.nan
    g["rank"] = np.arange(1, len(g) + 1)
    g["cum_share"] = g["share"].cumsum()
    return g


def coverage(contributions, levels: Sequence[float] = LEVELS, total: float | None = None) -> dict:
    """Number of items (sorted by contribution, largest first) needed to reach each level of ``total``.

    ``total`` defaults to the sum of the contributions (= the MI when every pointwise term is given,
    = the KL sum for marginal contributions).  Negative terms (depleted pairs) are counted in
    ``n_negative``; because they lower the total, the cumulative sum of the positive terms reaches
    every level before any negative item is needed.  Returns ``{'n_items', 'n_positive',
    'n_negative', 'total', 'n_50', 'n_80', 'n_95', ...}`` (``n_<level*100>``); a level that cannot
    be reached (total <= 0) gives ``nan``.
    """
    c = np.sort(np.asarray(list(contributions), dtype=float))[::-1]
    tot = float(c.sum()) if total is None else float(total)
    out = {"n_items": int(len(c)), "n_positive": int((c > 0).sum()), "n_negative": int((c < 0).sum()),
           "total": tot}
    cum = np.cumsum(c)
    for lv in levels:
        key = f"n_{int(round(lv * 100))}"
        if tot <= 0 or len(c) == 0:
            out[key] = float("nan")
        else:
            hit = np.flatnonzero(cum >= lv * tot - 1e-12)
            out[key] = int(hit[0] + 1) if len(hit) else float("nan")
    return out


def item_decomposition(x, y, x_name: str = "x", y_name: str = "y", levels: Sequence[float] = LEVELS) -> dict:
    """Pointwise table, marginal tables and coverage of I(X; Y) in one call.

    Returns ``{'mi', 'pairs', 'by_x', 'by_y', 'coverage'}`` where ``coverage`` is a DataFrame with one
    row per ``items`` in (``'pair'``, ``x_name``, ``y_name``) and the columns of :func:`coverage`.
    """
    pairs = pointwise_contributions(x, y, x_name, y_name)
    by_x = marginal_contributions(pairs, x_name)
    by_y = marginal_contributions(pairs, y_name)
    mi = float(pairs["contribution"].sum())
    rows = [{"items": "pair", **coverage(pairs["contribution"], levels)},
            {"items": x_name, **coverage(by_x["contribution"], levels)},
            {"items": y_name, **coverage(by_y["contribution"], levels)}]
    return {"mi": mi, "pairs": pairs, "by_x": by_x, "by_y": by_y, "coverage": pd.DataFrame(rows)}


# --------------------------------------------------------------------------- #
# position chain and token labels
# --------------------------------------------------------------------------- #


def _lip(df: pd.DataFrame) -> np.ndarray:
    lip = pd.to_numeric(df["line_in_par"], errors="coerce").to_numpy(dtype=float)
    return np.where(np.isnan(lip), 2, np.minimum(lip, 2)).astype(np.int64)


def position_chain(df: pd.DataFrame, y, shuffles: int = 0, seed: int = 0) -> dict:
    """Chain decomposition of I(Y; position class) of experiment 007 into three terms (bits).

    position class = place in line (first / mid / last / single) x ``rel_pos`` in 5 bins x
    ``line_in_par`` (0 / 1 / 2+), so that exactly

    ``I_pos = I_place + I_lip_given_place + I_rel_given_place_lip``

    with ``I_place`` = I(Y; place) (the edges of the line), ``I_lip_given_place`` = I(Y; line in
    paragraph | place) (paragraph-first line) and ``I_rel_given_place_lip`` (drift inside the line).
    With ``shuffles > 0`` every term is corrected by its mean under ``shuffles`` global permutations of
    ``y`` (the plug-in bias; the corrected terms are returned, the raw ones with suffix ``_raw``).
    """
    yc, _ = _codes(y)
    place, _ = _codes(line_place(df))
    lip = _lip(df)
    rel = np.clip((df["rel_pos"].to_numpy(dtype=float) * 5).astype(np.int64), 0, 4)
    pos, _ = _codes(position_class(df))
    pl_lip = _joint(place, lip)

    def terms(yy: np.ndarray) -> dict:
        return {
            "I_pos": _H(yy) + _H(pos) - _H(_joint(yy, pos)),
            "I_place": _H(yy) + _H(place) - _H(_joint(yy, place)),
            "I_lip_given_place": _H(_joint(yy, place)) + _H(pl_lip) - _H(_joint(yy, pl_lip)) - _H(place),
            "I_rel_given_place_lip": _H(_joint(yy, pl_lip)) + _H(_joint(pl_lip, rel)) - _H(_joint(yy, pl_lip, rel))
            - _H(pl_lip),
        }

    raw = terms(yc)
    out = {f"{k}_raw": v for k, v in raw.items()}
    if shuffles > 0:
        rng = np.random.default_rng(seed)
        acc = {k: 0.0 for k in raw}
        for _ in range(shuffles):
            for k, v in terms(yc[rng.permutation(len(yc))]).items():
                acc[k] += v / shuffles
        out.update({k: raw[k] - acc[k] for k in raw})
        out.update({f"{k}_bias": acc[k] for k in raw})
    else:
        out.update(raw)
    return out


def token_labels(df: pd.DataFrame, prefixes: Sequence[str] | None = None, suffixes: Sequence[str] | None = None,
                 preset: Any = "auto") -> dict[str, np.ndarray]:
    """Label arrays of every token: ``prefix``, ``suffix``, ``stem``, ``first_glyph``, ``last_glyph``.

    Affixes come from :func:`voytools.morphology.decompose_table` with the given lists (``None`` =
    learned on ``df``); empty affix = ``''``.  ``attrs`` of the result dictionary (key ``'lists'``) hold
    the lists used.
    """
    preset = resolve_preset(df, preset)
    dfd = df if ("prefix" in df.columns and prefixes is None and suffixes is None) else \
        decompose_table(df, prefixes=prefixes, suffixes=suffixes, preset=preset)
    gl = glyph_lists(df, preset)
    out = {
        "prefix": dfd["prefix"].to_numpy(dtype=object),
        "suffix": dfd["suffix"].to_numpy(dtype=object),
        "stem": dfd["stem"].to_numpy(dtype=object),
        "first_glyph": np.asarray([g[0] if g else "" for g in gl], dtype=object),
        "last_glyph": np.asarray([g[-1] if g else "" for g in gl], dtype=object),
    }
    out["lists"] = {"prefixes": list(dfd.attrs.get("prefixes", prefixes or [])),
                    "suffixes": list(dfd.attrs.get("suffixes", suffixes or []))}
    return out


def seam_pairs(df: pd.DataFrame, left, right) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(left[t], right[t+1], t)`` for every pair of consecutive tokens inside one line."""
    line = line_ids(df)
    same = np.flatnonzero(line[:-1] == line[1:])
    left = np.asarray(list(left), dtype=object)
    right = np.asarray(list(right), dtype=object)
    return left[same], right[same + 1], same


def seam_glyph_vs_affix(df: pd.DataFrame, labels: Mapping[str, np.ndarray], shuffles: int = 50,
                        seed: int = 0) -> dict:
    """Is the word seam only a glyph pair?  Shuffle-corrected seam information (bits) inside lines:

    ``I_glyph`` = I(last glyph_t; first glyph_{t+1}), ``I_affix`` = I(suffix_t; prefix_{t+1}),
    ``I_affix_given_glyph`` = I(suffix_t; prefix_{t+1} | last glyph_t, first glyph_{t+1}) and
    ``I_glyph_given_affix`` = I(last glyph_t; first glyph_{t+1} | suffix_t, prefix_{t+1}).  Each is
    the plug-in value minus its mean under ``shuffles`` permutations of the tokens within lines (all
    labels travel with the token), the correction of :func:`voytools.stats.edge_glyph_mi`.  Raw values
    are returned with suffix ``_raw``.  ``share_affix_residual`` = ``I_affix_given_glyph / I_affix``.
    """
    line = line_ids(df)
    n = len(df)
    same = line[:-1] == line[1:]
    lg, _ = _codes(labels["last_glyph"])
    fg, _ = _codes(labels["first_glyph"])
    su, _ = _codes(labels["suffix"])
    pr, _ = _codes(labels["prefix"])

    def terms(order: np.ndarray | None) -> dict:
        a, b, c, d = (lg, fg, su, pr) if order is None else (lg[order], fg[order], su[order], pr[order])
        L_g, R_g = a[:-1][same], b[1:][same]
        L_a, R_a = c[:-1][same], d[1:][same]
        gpair = _joint(L_g, R_g)
        apair = _joint(L_a, R_a)
        return {
            "I_glyph": _H(L_g) + _H(R_g) - _H(gpair),
            "I_affix": _H(L_a) + _H(R_a) - _H(apair),
            "I_affix_given_glyph": _H(_joint(L_a, gpair)) + _H(_joint(R_a, gpair)) - _H(_joint(apair, gpair)) - _H(gpair),
            "I_glyph_given_affix": _H(_joint(L_g, apair)) + _H(_joint(R_g, apair)) - _H(_joint(gpair, apair)) - _H(apair),
        }

    raw = terms(None)
    rng = np.random.default_rng(seed)
    acc = {k: 0.0 for k in raw}
    for _ in range(shuffles):
        order = np.lexsort((rng.random(n), line))
        for k, v in terms(order).items():
            acc[k] += v / shuffles
    out = {f"{k}_raw": v for k, v in raw.items()}
    out.update({k: raw[k] - acc[k] for k in raw} if shuffles else raw)
    out["n_pairs"] = int(same.sum())
    out["share_affix_residual"] = out["I_affix_given_glyph"] / out["I_affix"] if out["I_affix"] else float("nan")
    return out


# --------------------------------------------------------------------------- #
# stability across groups
# --------------------------------------------------------------------------- #


def _cond_table(x: np.ndarray, y: np.ndarray, kx: int, ky: int) -> np.ndarray:
    t = np.zeros((kx, ky), dtype=float)
    np.add.at(t, (x, y), 1.0)
    return t


def conditional_jsd(xa, ya, xb, yb, min_count: int = 20) -> dict:
    """Weighted Jensen-Shannon divergence between two conditional tables p_A(y | x) and p_B(y | x).

    ``sum_x w(x) JSD(p_A(y|x), p_B(y|x))`` over the values ``x`` with at least ``min_count`` pairs in
    *both* groups, ``w`` = pooled share of ``x`` among those values.  For an unconditional table pass a
    constant ``x``.  Returns ``{'jsd', 'n_contexts', 'n_a', 'n_b'}``.
    """
    xs, _ = _codes(list(xa) + list(xb))
    ys, _ = _codes(list(ya) + list(yb))
    na = len(xa)
    kx, ky = int(xs.max()) + 1 if len(xs) else 0, int(ys.max()) + 1 if len(ys) else 0
    ta = _cond_table(xs[:na], ys[:na], kx, ky)
    tb = _cond_table(xs[na:], ys[na:], kx, ky)
    ra, rb = ta.sum(1), tb.sum(1)
    ok = (ra >= min_count) & (rb >= min_count)
    if not ok.any():
        return {"jsd": float("nan"), "n_contexts": 0, "n_a": int(na), "n_b": int(len(xb))}
    w = (ra + rb)[ok]
    w = w / w.sum()
    vals = np.asarray([jsd(ta[i], tb[i]) for i in np.flatnonzero(ok)])
    return {"jsd": float((w * vals).sum()), "n_contexts": int(ok.sum()), "n_a": int(ra[ok].sum()),
            "n_b": int(rb[ok].sum())}


def _half_split(units: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Boolean mask: rows whose unit (line) fell into the first random half of the units."""
    u = np.unique(units)
    pick = rng.permutation(len(u))[: len(u) // 2]
    return np.isin(units, u[pick])


def stability(x, y, groups, units, pairs: Sequence[tuple[str, str]] | None = None, n_rep: int = 200,
              min_count: int = 20, rule_level: float = 0.8, min_group: int = 1, seed: int = 0,
              contexts: Sequence[Any] | None = None) -> pd.DataFrame:
    """Divergence of a rule table p(y | x) between groups versus between random halves of a group.

    ``x``, ``y`` are the labels of the table's rows (one row per token or per seam pair), ``groups``
    their group label (Currier language, scribe, section ...), ``units`` the unit of the random split
    (line id).  With ``contexts`` the divergences are computed on the rows whose ``x`` is in
    ``contexts`` only (for example the line-first tokens), while the rules (below) are still read off
    the whole table, so that enrichment is measured against the other contexts.  For every pair of
    groups (default: all pairs of groups with at least ``min_group`` rows, in the order of first
    appearance) the function reports

    * ``jsd_between``: :func:`conditional_jsd` of the two whole groups;
    * ``jsd_between_half``: mean over ``n_rep`` draws of the divergence between a random half of A and
      a random half of B (the same sample sizes as the within-group halves);
    * ``jsd_within_a``, ``jsd_within_b``: mean divergence between two random halves of the same group;
      ``jsd_within_hi`` = 97.5 % quantile of the pooled within-group draws;
    * ``excess`` = ``jsd_between_half - mean(jsd_within_a, jsd_within_b)`` and ``ratio`` =
      ``jsd_between_half / mean(within)``: 1 = the groups differ no more than halves of one group;
    * ``n_rules``, ``n_rules_shared``, ``share_shared``: the pointwise terms (pairs (x, y)) covering
      ``rule_level`` of the pooled sum of absolute contributions of the two groups, and how many of
      them have the same sign of ``log_ratio`` (enriched / depleted) in both groups and in the pool.
    """
    x = np.asarray(list(x), dtype=object)
    y = np.asarray(list(y), dtype=object)
    g = np.asarray(list(groups), dtype=object)
    units = np.asarray(list(units))
    if pairs is None:
        labels = [lab for lab in pd.unique(g) if (g == lab).sum() >= min_group]
        pairs = [(a, b) for i, a in enumerate(labels) for b in labels[i + 1:]]
    rng = np.random.default_rng(seed)
    rows = []
    for a, b in pairs:
        ma, mb = g == a, g == b
        if ma.sum() < 2 * min_count or mb.sum() < 2 * min_count:
            continue
        xa, ya, ua = x[ma], y[ma], units[ma]
        xb, yb, ub = x[mb], y[mb], units[mb]
        ca = np.ones(len(xa), dtype=bool) if contexts is None else np.isin(xa, np.asarray(list(contexts), dtype=object))
        cb = np.ones(len(xb), dtype=bool) if contexts is None else np.isin(xb, np.asarray(list(contexts), dtype=object))
        between = conditional_jsd(xa[ca], ya[ca], xb[cb], yb[cb], min_count)
        bh, wa, wb = [], [], []
        for _ in range(n_rep):
            ha = _half_split(ua, rng)
            hb = _half_split(ub, rng)
            a1, a2, b1, b2 = ha & ca, ~ha & ca, hb & cb, ~hb & cb
            bh.append(conditional_jsd(xa[a1], ya[a1], xb[b1], yb[b1], min_count)["jsd"])
            wa.append(conditional_jsd(xa[a1], ya[a1], xa[a2], ya[a2], min_count)["jsd"])
            wb.append(conditional_jsd(xb[b1], yb[b1], xb[b2], yb[b2], min_count)["jsd"])
        bh, wa, wb = (np.asarray(v, dtype=float) for v in (bh, wa, wb))
        within = np.r_[wa, wb]
        w_mean = float(np.nanmean(within)) if np.isfinite(within).any() else float("nan")
        b_mean = float(np.nanmean(bh)) if np.isfinite(bh).any() else float("nan")
        # shared rules: pointwise terms of the pooled table that cover rule_level of sum |contribution|
        pooled = pointwise_contributions(np.r_[xa, xb], np.r_[ya, yb])
        pooled = pooled.assign(abs_c=pooled["contribution"].abs()).sort_values("abs_c", ascending=False)
        cum = pooled["abs_c"].cumsum() / pooled["abs_c"].sum()
        rules = pooled[cum.shift(fill_value=0.0) < rule_level]
        if contexts is not None:
            rules = rules[rules["x"].isin(list(contexts))]
        pa = pointwise_contributions(xa, ya).set_index(["x", "y"])["log_ratio"]
        pb = pointwise_contributions(xb, yb).set_index(["x", "y"])["log_ratio"]
        shared = 0
        for _, r in rules.iterrows():
            key = (r["x"], r["y"])
            la, lb = pa.get(key, np.nan), pb.get(key, np.nan)
            s = np.sign(r["log_ratio"])
            if np.isfinite(la) and np.isfinite(lb) and np.sign(la) == s and np.sign(lb) == s:
                shared += 1
        rows.append({"group_a": a, "group_b": b, "n_a": int(ma.sum()), "n_b": int(mb.sum()),
                     "n_contexts": between["n_contexts"], "jsd_between": between["jsd"],
                     "jsd_between_half": b_mean, "jsd_between_half_lo": float(np.nanquantile(bh, 0.025)),
                     "jsd_between_half_hi": float(np.nanquantile(bh, 0.975)),
                     "jsd_within_a": float(np.nanmean(wa)), "jsd_within_b": float(np.nanmean(wb)),
                     "jsd_within_hi": float(np.nanquantile(within, 0.975)),
                     "excess": b_mean - w_mean, "ratio": b_mean / w_mean if w_mean else float("nan"),
                     "p_within_ge_between": float((within >= b_mean).mean()),
                     "n_rules": int(len(rules)), "n_rules_shared": int(shared),
                     "share_shared": shared / len(rules) if len(rules) else float("nan")})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# post-processor: R1 line/paragraph starts, R2 line ends, R3 seam
# --------------------------------------------------------------------------- #


def _r1_context(df: pd.DataFrame) -> np.ndarray:
    """R1 contexts: ``'first|pf'`` (first token of a paragraph-first line), ``'first|other'`` (first token
    of another line), ``'inner|pf'`` (other tokens of a paragraph-first line), ``''`` (no rule)."""
    first = df["is_first"].to_numpy(dtype=bool)
    pf = _lip(df) == 0
    return np.where(first & pf, "first|pf", np.where(first, "first|other", np.where(pf, "inner|pf", ""))).astype(object)


def _r2_mask(df: pd.DataFrame) -> np.ndarray:
    return df["is_last"].to_numpy(dtype=bool) & (df["n_tokens"].to_numpy(dtype=np.int64) >= 2)


def _r3_context(df: pd.DataFrame, last_glyph: np.ndarray) -> np.ndarray:
    """R3 context of every token: last glyph of the previous token of the same line (``''`` for the first)."""
    line = line_ids(df)
    ctx = np.full(len(df), "", dtype=object)
    inner = np.flatnonzero(np.r_[False, line[1:] == line[:-1]])
    ctx[inner] = last_glyph[inner - 1]
    return ctx


def _dist(values: np.ndarray) -> pd.Series:
    s = pd.Series(values, dtype=object).value_counts()
    return s / s.sum()


@dataclass
class PlaceRules:
    """Learned place tables (see :func:`learn_place_rules`).

    ``r1``: context (``'first|pf'``, ``'first|other'``, ``'inner|pf'``) -> prefix distribution;
    ``r2``: suffix distribution of the last token of lines with >= 2 tokens; ``r3``: last glyph of the
    previous token -> prefix distribution of inner tokens (contexts with >= ``min_context`` pairs).
    """

    prefixes: list = field(default_factory=list)
    suffixes: list = field(default_factory=list)
    preset: Any = None
    r1: dict = field(default_factory=dict)
    r2: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    r3: dict = field(default_factory=dict)
    min_context: int = 20
    n_tokens: int = 0

    def table(self, rule: str) -> pd.DataFrame:
        """Long table ``(context, item, p)`` of ``'R1'``, ``'R2'`` or ``'R3'``."""
        if rule == "R1":
            src = self.r1
        elif rule == "R2":
            src = {"last": self.r2}
        elif rule == "R3":
            src = self.r3
        else:
            raise ValueError("rule must be 'R1', 'R2' or 'R3'")
        rows = [{"rule": rule, "context": c, "item": k, "p": float(v)} for c, s in src.items() for k, v in s.items()]
        return pd.DataFrame(rows, columns=["rule", "context", "item", "p"])


def learn_place_rules(df: pd.DataFrame, prefixes: Sequence[str] | None = None, suffixes: Sequence[str] | None = None,
                      preset: Any = "auto", min_context: int = 20) -> PlaceRules:
    """Read the place tables R1, R2, R3 of a (real) token table in text order.

    Affixes come from :func:`voytools.morphology.decompose_table` with the given lists (``None`` =
    learned on ``df``).  R3 keeps only last-glyph contexts with at least ``min_context`` inner pairs.
    """
    preset = resolve_preset(df, preset)
    lab = token_labels(df, prefixes, suffixes, preset)
    ctx1 = _r1_context(df)
    r1 = {c: _dist(lab["prefix"][ctx1 == c]) for c in ("first|pf", "first|other", "inner|pf") if (ctx1 == c).any()}
    m2 = _r2_mask(df)
    r2 = _dist(lab["suffix"][m2]) if m2.any() else pd.Series(dtype=float)
    ctx3 = _r3_context(df, lab["last_glyph"])
    r3 = {}
    for g in pd.unique(ctx3[ctx3 != ""]):
        m = ctx3 == g
        if m.sum() >= min_context:
            r3[g] = _dist(lab["prefix"][m])
    return PlaceRules(prefixes=list(lab["lists"]["prefixes"]), suffixes=list(lab["lists"]["suffixes"]), preset=preset,
                      r1=r1, r2=r2, r3=r3, min_context=min_context, n_tokens=len(df))


def _minimal_rewrite(current: np.ndarray, target: pd.Series, rng: np.random.Generator) -> np.ndarray:
    """Change as few labels as possible so that the labels' distribution becomes ``target`` (in expectation).

    A label ``a`` is kept with probability min(1, q(a) / p(a)) (``p`` = empirical distribution of
    ``current``); the others draw from the excess distribution (q - p)_+ / sum.
    """
    n = len(current)
    if n == 0 or len(target) == 0:
        return current.copy()
    p = _dist(current)
    q = target / target.sum()
    items = sorted(set(p.index) | set(q.index))
    pv = np.asarray([p.get(k, 0.0) for k in items])
    qv = np.asarray([q.get(k, 0.0) for k in items])
    keep_p = np.where(pv > 0, np.minimum(1.0, qv / np.where(pv > 0, pv, 1.0)), 0.0)
    excess = np.maximum(qv - pv, 0.0)
    out = current.copy()
    code = {k: i for i, k in enumerate(items)}
    idx = np.asarray([code[a] for a in current])
    keep = rng.random(n) < keep_p[idx]
    repl = np.flatnonzero(~keep)
    if len(repl):
        if excess.sum() <= 0:
            return out
        draw = rng.choice(len(items), size=len(repl), p=excess / excess.sum())
        out[repl] = np.asarray(items, dtype=object)[draw]
    return out


def _replace_tokens(df: pd.DataFrame, tokens: Sequence[str], preset: Any, source: str) -> pd.DataFrame:
    out = df.copy()
    toks = list(tokens)
    out["token"] = toks
    out["token_raw"] = toks
    out["len_chars"] = [len(t) for t in toks]
    out["len_glyphs"] = [len(collapse_eva(t, preset)) for t in toks]
    out = out.drop(columns=[c for c in ("prefix", "stem", "suffix") if c in out.columns])
    out.attrs.update(df.attrs)
    out.attrs["source"] = source
    return out


def apply_place_rules(df: pd.DataFrame, rules: PlaceRules, which: Sequence[str] = ("R1", "R2", "R3"),
                      seed: int = 0) -> pd.DataFrame:
    """Copy of ``df`` (a generated or real corpus in text order) with the place rules imposed.

    Tokens are decomposed with the rules' affix lists; ``'R1'`` rewrites the prefix of the first token
    of every line and of every token of a paragraph-first line towards the learned context
    distribution, ``'R2'`` the suffix of the last token of every line (>= 2 tokens), ``'R3'`` the prefix
    of every inner token towards the distribution learned for the last glyph of the previous token
    (contexts without a table are left alone).  Each rewrite is the minimal coupling of
    :func:`_minimal_rewrite`, applied per context; stems are never changed, so the last glyph of a
    token that has a stem is unchanged and R3 contexts can be read off the input.  Order: R1, R2, R3
    (R3 overrides R1 for inner tokens of paragraph-first lines).  The rewritten token is
    ``prefix + stem + suffix``; decomposing it again with the same lists can give another split when
    the new affix and the stem form a longer affix (``'ka' + 'm'`` -> ``'k' + 'am'``).  The result has the layout and
    metadata of ``df``; ``attrs['source']`` = ``'<source>+<rules>'`` and ``attrs['place_rules']`` lists
    the rules applied.  Deterministic for a given ``seed``.
    """
    which = tuple(which)
    bad = [w for w in which if w not in ("R1", "R2", "R3")]
    if bad:
        raise ValueError(f"unknown rules {bad}; use 'R1', 'R2', 'R3'")
    preset = rules.preset
    lab = token_labels(df, rules.prefixes, rules.suffixes, preset)
    pre, stem, suf = lab["prefix"].copy(), lab["stem"], lab["suffix"].copy()
    rng = np.random.default_rng(seed)
    n_changed = {}
    if "R1" in which:
        ctx = _r1_context(df)
        before = pre.copy()
        for c, target in rules.r1.items():
            m = np.flatnonzero(ctx == c)
            if len(m):
                pre[m] = _minimal_rewrite(pre[m], target, rng)
        n_changed["R1"] = int((pre != before).sum())
    if "R2" in which and len(rules.r2):
        m = np.flatnonzero(_r2_mask(df))
        before = suf.copy()
        if len(m):
            suf[m] = _minimal_rewrite(suf[m], rules.r2, rng)
        n_changed["R2"] = int((suf != before).sum())
    if "R3" in which:
        ctx = _r3_context(df, lab["last_glyph"])
        before = pre.copy()
        for g, target in rules.r3.items():
            m = np.flatnonzero(ctx == g)
            if len(m):
                pre[m] = _minimal_rewrite(pre[m], target, rng)
        n_changed["R3"] = int((pre != before).sum())
    tokens = [p + s + x for p, s, x in zip(pre, stem, suf)]
    out = _replace_tokens(df, tokens, preset, f"{df.attrs.get('source', 'corpus')}+{'+'.join(which)}")
    out.attrs["place_rules"] = list(which)
    out.attrs["place_rules_changed"] = n_changed
    return out


# --------------------------------------------------------------------------- #
# conditional post-processor (experiment 017): R2', R3'a/b/c, R3'', type licence
# --------------------------------------------------------------------------- #

__all__ += ["ConditionalRules", "learn_conditional_rules", "apply_conditional_rules", "CONDITIONAL_RULES"]

CONDITIONAL_RULES = ("R1", "R2", "R2p", "R3", "R3a", "R3b", "R3c", "R3pp")
_SEP = "\t"
_PREFIX_RULES = ("R3", "R3a", "R3b", "R3c")


def _key(*parts: np.ndarray) -> np.ndarray:
    """One string key per row from several object arrays (``'\\t'``-joined)."""
    out = np.asarray(parts[0], dtype=object).astype(str).astype(object)
    for p in parts[1:]:
        out = out + _SEP + np.asarray(p, dtype=object).astype(str).astype(object)
    return out


def _glyphs_cached(strings, preset: Any) -> list[tuple[str, ...]]:
    cache: dict = {}
    out = []
    for s in strings:
        g = cache.get(s)
        if g is None:
            g = collapse_eva(s, preset) if s else ()
            cache[s] = g
        out.append(g)
    return out


def _first_glyphs(strings, preset: Any) -> np.ndarray:
    return np.asarray([g[0] if g else "" for g in _glyphs_cached(strings, preset)], dtype=object)


def _last_glyphs(strings, preset: Any) -> np.ndarray:
    return np.asarray([g[-1] if g else "" for g in _glyphs_cached(strings, preset)], dtype=object)


def _coarse_position(df: pd.DataFrame) -> np.ndarray:
    """Coarse position class of R3'c: ``first`` / ``inner`` / ``last`` token x ``pf`` / ``np``
    (paragraph-first line or not)."""
    first = df["is_first"].to_numpy(dtype=bool)
    last = df["is_last"].to_numpy(dtype=bool)
    place = np.where(first, "first", np.where(last, "last", "inner")).astype(object)
    pf = np.where(_lip(df) == 0, "pf", "np").astype(object)
    return place + "." + pf


def _count_table(keys: np.ndarray, values: np.ndarray, mask: np.ndarray | None = None) -> dict:
    """``key -> Series of counts of value`` (long table grouped by key)."""
    if mask is not None:
        keys, values = keys[mask], values[mask]
    if len(keys) == 0:
        return {}
    t = pd.DataFrame({"k": list(keys), "v": list(values)}).groupby(["k", "v"], sort=False).size()
    out: dict = {}
    for (k, v), n in t.items():
        out.setdefault(k, {})[v] = int(n)
    return {k: pd.Series(d, dtype=float) for k, d in out.items()}


def _counts(values: np.ndarray) -> pd.Series:
    s = pd.Series(list(values), dtype=object).value_counts()
    return s.astype(float)


@dataclass
class ConditionalRules:
    """Learned place tables of experiment 017 (see :func:`learn_conditional_rules`).

    ``base`` holds R1-R3 of 014 (:class:`PlaceRules`, used unchanged for the sets ``+R1+R2`` and
    ``+R1+R2+R3``).  ``tables`` are count tables ``name -> {context key -> Series of counts}`` of the
    conditional prefix rules of inner tokens (tokens with a predecessor in the line): ``'R3'`` prefix |
    last glyph of the previous token, ``'R3a'`` prefix | (last glyph before, suffix), ``'R3b'`` prefix |
    (last glyph before, first glyph of the stem), ``'R3c'`` prefix | (last glyph before, suffix, coarse
    position), and of ``'R3pp'`` first glyph of the stem | last glyph before, learned on tokens without
    prefix.  ``uncond`` holds the unconditional fallbacks (``'R3'``: prefixes of inner tokens, ``'R3pp'``:
    first stem glyphs of inner tokens without prefix).  A context cell is used only with at least
    ``min_cell`` observations; otherwise the rewrite backs off along ``R3c -> R3a -> R3 -> uncond``
    (``R3b -> R3 -> uncond``, ``R3pp -> uncond``).  ``r2_inner`` is the last-glyph distribution of the tokens
    that are not last in their line and ``end_glyphs`` the glyphs enriched at the line end more than
    ``end_ratio`` times (R2').  ``vocab`` (token -> count of the real text) serves the type licence.
    """

    prefixes: list = field(default_factory=list)
    suffixes: list = field(default_factory=list)
    preset: Any = None
    base: PlaceRules = field(default_factory=PlaceRules)
    min_cell: int = 5
    tables: dict = field(default_factory=dict)
    uncond: dict = field(default_factory=dict)
    r2_inner: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    end_glyphs: list = field(default_factory=list)
    end_ratio: float = 5.0
    vocab: dict = field(default_factory=dict)
    licence_min_count: int = 2
    n_tokens: int = 0

    def supported(self, name: str, key: str) -> bool:
        """Does the cell ``key`` of table ``name`` have at least ``min_cell`` observations?"""
        s = self.tables.get(name, {}).get(key)
        return s is not None and float(s.sum()) >= self.min_cell

    def support(self) -> pd.DataFrame:
        """One row per table: number of contexts, contexts with >= ``min_cell`` observations and the share
        of the (real) tokens that sit in supported cells."""
        rows = []
        for name, tab in self.tables.items():
            ns = np.asarray([float(s.sum()) for s in tab.values()])
            ok = ns >= self.min_cell
            rows.append({"table": name, "n_contexts": int(len(ns)), "n_contexts_supported": int(ok.sum()),
                         "n_tokens": int(ns.sum()), "share_tokens_supported": float(ns[ok].sum() / ns.sum()) if ns.sum() else np.nan,
                         "min_cell": self.min_cell})
        return pd.DataFrame(rows)

    def table(self, name: str) -> pd.DataFrame:
        """Long table ``(rule, context, item, n, p)`` of one conditional table."""
        rows = []
        for c, s in self.tables.get(name, {}).items():
            tot = float(s.sum())
            for k, v in s.items():
                rows.append({"rule": name, "context": c.replace(_SEP, "|"), "item": k, "n": int(v), "p": float(v) / tot})
        return pd.DataFrame(rows, columns=["rule", "context", "item", "n", "p"])


def learn_conditional_rules(df: pd.DataFrame, prefixes: Sequence[str] | None = None,
                            suffixes: Sequence[str] | None = None, preset: Any = "auto", min_cell: int = 5,
                            min_context: int = 20, end_ratio: float = 5.0, end_min_count: int = 20,
                            licence_min_count: int = 2) -> ConditionalRules:
    """Read the conditional place tables of experiment 017 from a token table in text order.

    ``min_context`` is the R3 threshold of the 014 tables (``base``); ``min_cell`` the backoff threshold of
    the conditional tables.  ``end_ratio`` / ``end_min_count``: a glyph belongs to the line-end set when it
    occurs at least ``end_min_count`` times as the last glyph of a line-last token and its share there is
    more than ``end_ratio`` times its share among the other tokens.
    """
    preset = resolve_preset(df, preset)
    base = learn_place_rules(df, prefixes, suffixes, preset, min_context)
    lab = token_labels(df, base.prefixes, base.suffixes, preset)
    pre, stem, suf, lg = lab["prefix"], lab["stem"], lab["suffix"], lab["last_glyph"]
    sfg = _first_glyphs(stem, preset)
    ctx3 = _r3_context(df, lg)
    inner = ctx3 != ""
    pos = _coarse_position(df)
    tables = {
        "R3": _count_table(_key(ctx3), pre, inner),
        "R3a": _count_table(_key(ctx3, suf), pre, inner),
        "R3b": _count_table(_key(ctx3, sfg), pre, inner),
        "R3c": _count_table(_key(ctx3, suf, pos), pre, inner),
    }
    nopre = inner & (pre == "") & (stem != "")
    tables["R3pp"] = _count_table(_key(ctx3), sfg, nopre)
    uncond = {"R3": _counts(pre[inner]), "R3pp": _counts(sfg[nopre])}
    last = _r2_mask(df)
    not_last = ~df["is_last"].to_numpy(dtype=bool)
    r2_inner = _dist(lg[not_last]) if not_last.any() else pd.Series(dtype=float)
    pl_n = pd.Series(list(lg[last]), dtype=object).value_counts()
    pl = pl_n / pl_n.sum() if pl_n.sum() else pl_n
    end_glyphs = []
    for g, n in pl_n.items():
        if n < end_min_count:
            continue
        po = float(r2_inner.get(g, 0.0))
        if po == 0.0 or pl[g] / po > end_ratio:
            end_glyphs.append(g)
    vocab = {k: int(v) for k, v in df["token"].value_counts().items()}
    return ConditionalRules(prefixes=list(base.prefixes), suffixes=list(base.suffixes), preset=preset, base=base,
                            min_cell=min_cell, tables=tables, uncond=uncond, r2_inner=r2_inner,
                            end_glyphs=sorted(end_glyphs), end_ratio=end_ratio, vocab=vocab,
                            licence_min_count=licence_min_count, n_tokens=len(df))


def _minimal_rewrite_ex(current: np.ndarray, target: pd.Series, rng: np.random.Generator,
                        left: Sequence[str] | None = None, right: Sequence[str] | None = None,
                        accept=None, max_tries: int = 5) -> tuple[np.ndarray, int]:
    """:func:`_minimal_rewrite` with an optional type licence.

    When ``accept`` is given, a replacement label ``a`` of token ``i`` is accepted only if
    ``accept(left[i] + a + right[i])`` is true; otherwise it is drawn again (at most ``max_tries`` draws),
    and the original label is kept when every draw fails.  Returns ``(labels, n_rejected)``.
    """
    n = len(current)
    if n == 0 or len(target) == 0 or float(target.sum()) <= 0:
        return current.copy(), 0
    p = _dist(current)
    q = target / target.sum()
    items = sorted(set(p.index) | set(q.index))
    pv = np.asarray([p.get(k, 0.0) for k in items])
    qv = np.asarray([q.get(k, 0.0) for k in items])
    keep_p = np.where(pv > 0, np.minimum(1.0, qv / np.where(pv > 0, pv, 1.0)), 0.0)
    excess = np.maximum(qv - pv, 0.0)
    out = current.copy()
    if excess.sum() <= 0:
        return out, 0
    code = {k: i for i, k in enumerate(items)}
    idx = np.asarray([code[a] for a in current])
    keep = rng.random(n) < keep_p[idx]
    pending = np.flatnonzero(~keep)
    arr = np.asarray(items, dtype=object)
    pe = excess / excess.sum()
    tries = 1 if accept is None else max(1, int(max_tries))
    for _ in range(tries):
        if not len(pending):
            break
        draw = arr[rng.choice(len(items), size=len(pending), p=pe)]
        if accept is None:
            out[pending] = draw
            pending = pending[:0]
            break
        ok = np.asarray([accept(left[i] + a + right[i]) for i, a in zip(pending, draw)], dtype=bool)
        out[pending[ok]] = draw[ok]
        pending = pending[~ok]
    return out, int(len(pending))


def _backoff_levels(rules: ConditionalRules, rule: str, keys: Mapping[str, np.ndarray]) -> list[tuple[str, np.ndarray | None]]:
    """Backoff chain of a conditional rule: ``[(table name, key array), ..., ('uncond:<name>', None)]``."""
    chains = {"R3": ["R3"], "R3a": ["R3a", "R3"], "R3b": ["R3b", "R3"], "R3c": ["R3c", "R3a", "R3"], "R3pp": ["R3pp"]}
    out: list = [(name, keys[name]) for name in chains[rule]]
    out.append(("uncond:" + ("R3pp" if rule == "R3pp" else "R3"), None))
    return out


def _rewrite_with_backoff(rules: ConditionalRules, rule: str, idx: np.ndarray, current: np.ndarray,
                          keys: Mapping[str, np.ndarray], rng: np.random.Generator, left, right, accept,
                          max_tries: int) -> tuple[np.ndarray, dict]:
    """Rewrite ``current[idx]`` towards the most specific supported table of ``rule``; returns the new labels
    of the whole array and a report ``{'n_candidates', 'n_changed', 'n_rejected', 'support': {level: n}}``."""
    out = current.copy()
    level_of = np.full(len(idx), -1, dtype=np.int64)
    level_key = np.empty(len(idx), dtype=object)
    levels = _backoff_levels(rules, rule, keys)
    pending = np.arange(len(idx))
    for li, (name, karr) in enumerate(levels):
        if not len(pending):
            break
        if karr is None:
            level_of[pending] = li
            level_key[pending] = ""
            pending = pending[:0]
            break
        kk = karr[idx[pending]]
        ok = np.asarray([rules.supported(name, k) for k in kk], dtype=bool)
        level_of[pending[ok]] = li
        level_key[pending[ok]] = kk[ok]
        pending = pending[~ok]
    support = {}
    n_rej = 0
    for li, (name, karr) in enumerate(levels):
        sel = np.flatnonzero(level_of == li)
        support[name] = int(len(sel))
        if not len(sel):
            continue
        if karr is None:
            groups = {"": sel}
        else:
            groups = {}
            for j in sel:
                groups.setdefault(level_key[j], []).append(j)
            groups = {k: np.asarray(v, dtype=np.int64) for k, v in groups.items()}
        for k, g in groups.items():
            rows = idx[g]
            target = rules.uncond[name.split(":")[1]] if karr is None else rules.tables[name][k]
            new, rej = _minimal_rewrite_ex(out[rows], target, rng, [left[r] for r in rows] if accept else None,
                                           [right[r] for r in rows] if accept else None, accept, max_tries)
            out[rows] = new
            n_rej += rej
    return out, {"n_candidates": int(len(idx)), "n_changed": int((out[idx] != current[idx]).sum()),
                 "n_rejected": n_rej, "support": support}


def apply_conditional_rules(df: pd.DataFrame, rules: ConditionalRules, which: Sequence[str] = ("R1", "R2p", "R3a"),
                            seed: int = 0, licence: bool = False, max_tries: int = 5) -> pd.DataFrame:
    """Copy of ``df`` (a generated or real corpus in text order) with the conditional place rules imposed.

    ``which`` is a subset of :data:`CONDITIONAL_RULES`: ``'R1'`` and ``'R2'`` as in :func:`apply_place_rules`
    (tables of ``rules.base``); ``'R2p'`` = R2 plus, for the tokens that are not last in their line, the last
    glyph is rewritten when it belongs to ``rules.end_glyphs`` (kept with probability
    min(1, q(g) / p(g)), ``q`` = real share among non-last tokens, otherwise replaced by a draw from the real
    non-last distribution over the other glyphs; the glyph is changed in the suffix when the token has one,
    else in the stem); ``'R3'``, ``'R3a'``, ``'R3b'``, ``'R3c'`` (at most one) rewrite the prefix of every inner
    token towards the most specific supported table with backoff (see :class:`ConditionalRules`), per
    context the minimal coupling of :func:`_minimal_rewrite`; ``'R3pp'`` rewrites the first glyph of the
    stem of inner tokens without prefix (after the prefix rule).  Contexts (last glyph of the token before,
    own suffix, first glyph of the stem, coarse position) are read after R1 / R2 / R2' from the rewritten
    parts, so the prefix rules see the suffixes they condition on; the previous token's last glyph is read
    once before the prefix rules (a prefix rewrite changes it only for tokens with empty stem and suffix).

    With ``licence=True`` every rewrite is accepted only when the resulting token occurs at least
    ``rules.licence_min_count`` times in ``rules.vocab``; a rejected draw is repeated at most ``max_tries``
    times, then the original part is kept.  The result has the layout of ``df``; ``attrs['place_rules']``
    lists the rules, ``attrs['place_rules_changed']`` the number of changed tokens per rule,
    ``attrs['rewrite_report']`` per rule ``n_candidates`` / ``n_changed`` / ``n_rejected`` and the backoff
    ``support`` (tokens per table level), ``attrs['licence']`` the switch.  Deterministic for a given ``seed``.
    """
    which = tuple(which)
    bad = [w for w in which if w not in CONDITIONAL_RULES]
    if bad:
        raise ValueError(f"unknown rules {bad}; use {CONDITIONAL_RULES}")
    pr = [w for w in which if w in _PREFIX_RULES]
    if len(pr) > 1:
        raise ValueError(f"at most one prefix rule of {_PREFIX_RULES}, got {pr}")
    preset = rules.preset
    lab = token_labels(df, rules.prefixes, rules.suffixes, preset)
    pre, stem, suf = lab["prefix"].copy(), lab["stem"].copy(), lab["suffix"].copy()
    n = len(df)
    rng = np.random.default_rng(seed)
    accept = (lambda tok: rules.vocab.get(tok, 0) >= rules.licence_min_count) if licence else None
    report: dict = {}
    changed: dict = {}

    def tokens_now() -> list:
        return [p + s + x for p, s, x in zip(pre, stem, suf)]

    if "R1" in which:
        ctx = _r1_context(df)
        before = pre.copy()
        rej = 0
        for c, target in rules.base.r1.items():
            m = np.flatnonzero(ctx == c)
            if len(m):
                new, r = _minimal_rewrite_ex(pre[m], target, rng, [""] * len(m) if accept else None,
                                             [stem[i] + suf[i] for i in m] if accept else None, accept, max_tries)
                pre[m] = new
                rej += r
        changed["R1"] = int((pre != before).sum())
        report["R1"] = {"n_candidates": int((ctx != "").sum()), "n_changed": changed["R1"], "n_rejected": rej}
    if ("R2" in which or "R2p" in which) and len(rules.base.r2):
        m = np.flatnonzero(_r2_mask(df))
        before = suf.copy()
        rej = 0
        if len(m):
            new, rej = _minimal_rewrite_ex(suf[m], rules.base.r2, rng, [pre[i] + stem[i] for i in m] if accept else None,
                                           [""] * len(m) if accept else None, accept, max_tries)
            suf[m] = new
        changed["R2"] = int((suf != before).sum())
        report["R2"] = {"n_candidates": int(len(m)), "n_changed": changed["R2"], "n_rejected": rej}
    if "R2p" in which and len(rules.end_glyphs) and len(rules.r2_inner):
        toks = tokens_now()
        gl = _glyphs_cached(toks, preset)
        lg = np.asarray([g[-1] if g else "" for g in gl], dtype=object)
        not_last = ~df["is_last"].to_numpy(dtype=bool)
        E = list(rules.end_glyphs)
        cand = np.flatnonzero(not_last & np.isin(lg, E))
        p_all = _dist(lg[not_last]) if not_last.any() else pd.Series(dtype=float)
        keep_p = {g: min(1.0, float(rules.r2_inner.get(g, 0.0)) / float(p_all.get(g, 1.0))) if p_all.get(g, 0.0) > 0 else 0.0
                  for g in E}
        q_rest = rules.r2_inner.drop([g for g in E if g in rules.r2_inner.index])
        n_ch = rej = 0
        if len(cand) and len(q_rest) and float(q_rest.sum()) > 0:
            q_rest = q_rest / q_rest.sum()
            items = np.asarray(list(q_rest.index), dtype=object)
            keep = rng.random(len(cand)) < np.asarray([keep_p[g] for g in lg[cand]])
            for i in cand[~keep]:
                g_tok = gl[i]
                sg = collapse_eva(suf[i], preset) if suf[i] else ()
                in_suffix = bool(sg) and sg[-1] == g_tok[-1]
                base_part = suf[i] if in_suffix else stem[i]
                if not base_part:
                    continue
                parts = collapse_eva(base_part, preset)
                head = "".join(parts[:-1])
                left = (pre[i] + stem[i] + head) if in_suffix else (pre[i] + head)
                right = "" if in_suffix else suf[i]
                done = False
                for _ in range(max(1, max_tries) if accept else 1):
                    g_new = items[rng.choice(len(items), p=q_rest.to_numpy())]
                    if accept is None or accept(left + g_new + right):
                        done = True
                        break
                if not done:
                    rej += 1
                    continue
                if in_suffix:
                    suf[i] = head + g_new
                else:
                    stem[i] = head + g_new
                n_ch += 1
        changed["R2p"] = n_ch
        report["R2p"] = {"n_candidates": int(len(cand)), "n_changed": n_ch, "n_rejected": rej,
                         "end_glyphs": list(E)}
    if pr or "R3pp" in which:
        toks = tokens_now()
        lg = _last_glyphs(toks, preset)
        ctx3 = _r3_context(df, lg)
        inner = np.flatnonzero(ctx3 != "")
        sfg = _first_glyphs(stem, preset)
        keys = {"R3": _key(ctx3), "R3a": _key(ctx3, suf), "R3b": _key(ctx3, sfg), "R3c": _key(ctx3, suf, _coarse_position(df)),
                "R3pp": _key(ctx3)}
        if pr:
            rule = pr[0]
            left = [""] * n
            right = [stem[i] + suf[i] for i in range(n)] if accept else None
            new, rep = _rewrite_with_backoff(rules, rule, inner, pre, keys, rng, left, right, accept, max_tries)
            pre = new
            changed[rule] = rep["n_changed"]
            report[rule] = rep
        if "R3pp" in which:
            nopre = inner[(pre[inner] == "") & (stem[inner] != "")]
            sg = _glyphs_cached(stem, preset)
            tails = ["".join(g[1:]) for g in sg]
            left = [pre[i] for i in range(n)]
            right = [tails[i] + suf[i] for i in range(n)] if accept else None
            new, rep = _rewrite_with_backoff(rules, "R3pp", nopre, sfg, keys, rng, left, right, accept, max_tries)
            ch = np.flatnonzero(new != sfg)
            for i in ch:
                stem[i] = new[i] + tails[i]
            changed["R3pp"] = rep["n_changed"]
            report["R3pp"] = rep
    tokens = tokens_now()
    out = _replace_tokens(df, tokens, preset, f"{df.attrs.get('source', 'corpus')}+{'+'.join(which)}" + ("+L" if licence else ""))
    out.attrs["place_rules"] = list(which)
    out.attrs["place_rules_changed"] = changed
    out.attrs["rewrite_report"] = report
    out.attrs["licence"] = bool(licence)
    out.attrs["n_tokens_changed"] = int(sum(t != o for t, o in zip(tokens, df["token"].tolist())))
    return out
