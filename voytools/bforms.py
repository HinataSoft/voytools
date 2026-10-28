"""Forms found only in one Currier language (experiment 019): classes, profiles, the ``-ey`` / ``-edy`` choice.

Experiment 016 showed that about one fifth of the B tokens are whole forms that A does not have and
that are not glyph variants of A forms; the only systematic operation between the two vocabularies
is the insertion of ``d`` (``-ey`` -> ``-edy``, v101 ``8``).  This module splits that fifth into
classes and measures how each class behaves against frequency-matched shared forms:

* :func:`neighbours_within` - every candidate within a Levenshtein distance (same symmetric-delete
  index as :func:`voytools.registers.nearest_counterparts`, exact up to ``max_dist``);
* :func:`is_d_op` - is an edit operation one of the Currier "writing" operations (insertion /
  deletion of ``e``, ``ee``, ``d``, ``q`` or ``ch`` <-> ``sh``; v101 equivalents);
* :func:`classify_forms` - one row per form of X: ``shared`` (>= 1 token in Y), ``d_variant``
  (X-only, a Y counterpart reached by one d-operation), ``other_le1`` (X-only, a Y counterpart within
  one edit that is not a d-operation), ``no_counterpart`` (X-only, nothing within one edit), ``hapax``
  (one token in X, none in Y); with the stem of each form and whether the stem is new (no token in Y);
* :func:`composition` - token shares of the classes and of new stems among the X-only tokens;
* :func:`match_by_count` - shared forms drawn with the token-count distribution of a class
  (frequency-matched baseline);
* :func:`form_profile` - paradigm width, dispersion over pages (pages per token, Gries' DP, top-page
  share), section entropy, label share, line-edge shares, co-occurrence of two class forms in one
  line, glyph-bigram log-likelihood under a model of the shared forms, and the seam (first unit of
  the next token) of a set of forms;
* :func:`bigram_model`, :func:`bigram_loglik`, :func:`gries_dp`, :func:`entropy_bits` - helpers;
* :func:`suffix_pairs`, :func:`pair_choice` - tokens ending in one of two alternative suffixes
  (``-ey`` / ``-edy``), the probability of the second form given a context, the entropy of the choice
  with and without the context and its permutation null; :func:`both_forms` - stems with both forms.

All functions take plain arrays / mappings so that the same code runs on the real text, on halves of
one language, on generated nulls and on Latin controls.  Forms are tuples of units (glyphs of
``collapse_eva``, v101 characters, letters) as in :mod:`voytools.registers`.
"""

from __future__ import annotations

from collections import Counter
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .registers import _deletes, edit_script, levenshtein

Units = tuple

#: Currier "writing" operations per unit scheme: insertions / deletions made only of these units, or the
#: listed substitutions (both directions).  ``extended`` = Eva glyphs of the ``extended`` preset, ``v101``
#: = characters of the v101 alphabet (``c`` = e, ``C`` = ee, ``8`` / ``7`` = d, ``4`` = q, ``1`` / ``2`` = ch / sh).
D_OPS = {
    "extended": {"units": {"e", "ee", "d", "q"}, "subs": {frozenset({"ch", "sh"})}},
    "v101": {"units": {"c", "C", "8", "7", "4"}, "subs": {frozenset({"1", "2"})}},
}

CLASSES = ("shared", "d_variant", "other_le1", "no_counterpart", "hapax")
ONLY_CLASSES = ("d_variant", "other_le1", "no_counterpart")


# --------------------------------------------------------------------------- #
# counterparts and classes
# --------------------------------------------------------------------------- #


def neighbours_within(queries: Iterable[Units], candidates: Iterable[Units], max_dist: int = 2) -> dict:
    """Every candidate within ``max_dist`` edits of every query: ``{query: [(candidate, distance), ...]}``.

    Exact for distances <= ``max_dist`` (symmetric-delete index); the list is sorted by distance,
    then by the candidate tuple.

    >>> neighbours_within([("o", "k", "y")], [("o", "k", "e", "y"), ("d", "y"), ("o", "k", "y", "d")], 1)
    {('o', 'k', 'y'): [(('o', 'k', 'e', 'y'), 1), (('o', 'k', 'y', 'd'), 1)]}
    """
    index: dict = {}
    for c in dict.fromkeys(tuple(c) for c in candidates):
        for d in _deletes(c, max_dist):
            index.setdefault(d, []).append(c)
    out = {}
    for q in queries:
        q = tuple(q)
        seen: set = set()
        hits = []
        for d in _deletes(q, max_dist):
            for c in index.get(d, ()):
                if c in seen:
                    continue
                seen.add(c)
                dist = 0 if c == q else levenshtein(q, c, max_dist)
                if dist <= max_dist:
                    hits.append((c, dist))
        hits.sort(key=lambda t: (t[1], t[0]))
        out[q] = hits
    return out


def is_d_op(op: str, scheme: str | None = "extended") -> bool:
    """Is the edit operation (``kind:units@position`` of :func:`voytools.registers.edit_script`) a Currier operation?

    Insertion or deletion made only of the scheme's d-units (``e``, ``ee``, ``d``, ``q``; v101 ``c``,
    ``C``, ``8``, ``7``, ``4``) or one of its substitutions (``ch`` <-> ``sh``; v101 ``1`` <-> ``2``).
    ``scheme=None`` (letters): never.

    >>> is_d_op("ins:d@mid"), is_d_op("sub:ch>sh@begin"), is_d_op("ins:l@begin"), is_d_op("ins:8@mid", "v101")
    (True, True, False, True)
    """
    if scheme is None or scheme not in D_OPS:
        return False
    spec = D_OPS[scheme]
    kind, rest = op.split(":", 1)
    units, _ = rest.rsplit("@", 1)
    if kind == "sub":
        return frozenset(units.split(">")) in spec["subs"]
    if scheme == "extended":
        from .glyphs import collapse_eva
        parts = collapse_eva(units, "extended")
    else:
        parts = tuple(units)
    return all(u in spec["units"] for u in parts)


def classify_forms(
    count_x: Mapping[Units, float], count_y: Mapping[Units, float], stem_of: Mapping[Units, str],
    scheme: str | None = "extended", min_count: int = 2, max_dist: int = 2,
) -> pd.DataFrame:
    """Classes of the forms of X against the vocabulary of Y.

    ``shared``: >= 1 token in Y; ``hapax``: 1 token in X, none in Y; the X-only forms (>= ``min_count``
    in X, none in Y) are ``d_variant`` if a Y form with >= ``min_count`` tokens is reached by a single
    d-operation (:func:`is_d_op`; merged insertions such as ``ins:ee`` count as one operation, hence
    ``max_dist`` 2), else ``other_le1`` if some Y form is within one edit, else ``no_counterpart``.
    Among several candidates the one with more Y tokens wins.  The edit script is written Y -> X
    (A -> B for the B-only forms).

    ``stem_of`` maps every X form to its stem (decomposition); ``stem_new`` is True when no Y token
    has that stem (summed over the Y forms present in ``stem_of``; Y forms missing from ``stem_of``
    are ignored).  Returns one row per X form: ``form`` (tuple), ``n_x``, ``n_y``, ``cls``,
    ``counterpart`` (tuple or None), ``distance``, ``op``, ``stem``, ``stem_new``.
    """
    cand = {s: c for s, c in count_y.items() if c >= min_count}
    only = [s for s, c in count_x.items() if c >= min_count and not count_y.get(s)]
    near = neighbours_within(only, cand, max_dist) if only and cand else {s: [] for s in only}
    stem_y: dict = {}
    for s, c in count_y.items():
        st = stem_of.get(s)
        if st is not None:
            stem_y[st] = stem_y.get(st, 0) + c
    rows = []
    for s, c in count_x.items():
        ny = count_y.get(s, 0)
        cp, dist, op, cls = None, np.nan, "", "shared"
        if not ny:
            if c < min_count:
                cls = "hapax"
            else:
                best_d, best_1 = None, None
                for cnd, d in near.get(s, []):
                    script = edit_script(cnd, s)
                    if len(script) == 1 and is_d_op(script[0], scheme):
                        if best_d is None or cand[cnd] > cand[best_d[0]]:
                            best_d = (cnd, d, script[0])
                    elif d <= 1 and (best_1 is None or cand[cnd] > cand[best_1[0]]):
                        best_1 = (cnd, d, ";".join(script))
                if best_d is not None:
                    cp, dist, op = best_d
                    cls = "d_variant"
                elif best_1 is not None:
                    cp, dist, op = best_1
                    cls = "other_le1"
                else:
                    cls = "no_counterpart"
        st = stem_of.get(s, "")
        rows.append({"form": s, "n_x": c, "n_y": ny, "cls": cls, "counterpart": cp, "distance": dist, "op": op,
                     "stem": st, "stem_new": bool(st != "" and not stem_y.get(st))})
    return pd.DataFrame(rows, columns=["form", "n_x", "n_y", "cls", "counterpart", "distance", "op", "stem", "stem_new"])


def composition(classes: pd.DataFrame) -> dict:
    """Token shares of the classes of :func:`classify_forms` and of new stems among the X-only tokens.

    ``share_<cls>`` = tokens of the class / all X tokens; ``share_only`` = the three X-only classes
    together; ``only_<cls>`` = tokens of the class / X-only tokens; ``only_new_stem`` = X-only tokens
    (three classes) whose stem is new / X-only tokens, ``only_new_form`` its complement;
    ``hapax_new_stem`` the same among hapaxes; ``n_forms_<cls>`` type counts.
    """
    tot = float(classes["n_x"].sum())
    out: dict = {"n_tokens": tot, "n_forms": int(len(classes))}
    by = classes.groupby("cls")["n_x"].sum()
    only = float(sum(by.get(c, 0) for c in ONLY_CLASSES))
    out["share_only"] = only / tot if tot else np.nan
    for c in CLASSES:
        t = float(by.get(c, 0))
        out[f"share_{c}"] = t / tot if tot else np.nan
        out[f"n_forms_{c}"] = int((classes["cls"] == c).sum())
    for c in ONLY_CLASSES:
        out[f"only_{c}"] = float(by.get(c, 0)) / only if only else np.nan
    sel = classes["cls"].isin(ONLY_CLASSES)
    new = float(classes.loc[sel & classes["stem_new"], "n_x"].sum())
    out["only_new_stem"] = new / only if only else np.nan
    out["only_new_form"] = 1 - out["only_new_stem"] if only else np.nan
    for c in ONLY_CLASSES:
        s = classes["cls"] == c
        t = float(classes.loc[s, "n_x"].sum())
        out[f"{c}_new_stem"] = float(classes.loc[s & classes["stem_new"], "n_x"].sum()) / t if t else np.nan
    h = classes["cls"] == "hapax"
    out["hapax_new_stem"] = float(classes.loc[h & classes["stem_new"], "n_x"].sum()) / float(classes.loc[h, "n_x"].sum()) \
        if h.any() else np.nan
    return out


def match_by_count(classes: pd.DataFrame, cls: str, rng: np.random.Generator, pool_cls: str = "shared",
                   tol: float = 0.25) -> tuple[list, list]:
    """Forms of ``cls`` and shared forms matched to them by token count (without replacement).

    For every form of the class (random order) a shared form with the same ``n_x`` is drawn; when none
    is left, the shared form with the nearest count, accepted only if it differs by at most
    ``tol * n_x`` (exact for counts up to 3 with the default); forms without a match are dropped, so
    the two lists are shorter than the class when the pool runs out (hapaxes: shared forms with one
    token).  Returns ``(class_forms, matched_forms)`` as lists of tuples of equal length.
    """
    target = classes.loc[classes["cls"] == cls, ["form", "n_x"]]
    pool = classes.loc[classes["cls"] == pool_cls, ["form", "n_x"]]
    if target.empty or pool.empty:
        return [], []
    by_count: dict = {}
    for f, c in zip(pool["form"], pool["n_x"]):
        by_count.setdefault(int(c), []).append(f)
    for lst in by_count.values():
        rng.shuffle(lst)
    counts = sorted(by_count)
    chosen, forms = [], []
    order = rng.permutation(len(target))
    for i in order:
        f, c = target["form"].iloc[i], int(target["n_x"].iloc[i])
        if not counts:
            break
        if c in by_count and by_count[c]:
            k = c
        else:
            k = min(counts, key=lambda v: (abs(v - c), rng.random()))
            if abs(k - c) > tol * c:
                continue
        chosen.append(by_count[k].pop())
        forms.append(f)
        if not by_count[k]:
            del by_count[k]
            counts = sorted(by_count)
    return forms, chosen


# --------------------------------------------------------------------------- #
# profiles
# --------------------------------------------------------------------------- #


def entropy_bits(counts) -> float:
    """Shannon entropy in bits of a count vector (0 for an empty or single-valued vector)."""
    c = np.asarray(list(counts), dtype=float)
    c = c[c > 0]
    if c.size == 0:
        return 0.0
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def gries_dp(form_counts: Mapping, sizes: Mapping) -> float:
    """Gries' deviation of proportions of one form over parts (pages): 0.5 * sum |v_i - s_i|.

    ``form_counts`` = tokens of the form per part, ``sizes`` = tokens of every part; 0 = perfectly
    even, 1 = concentrated in a vanishing part.

    >>> round(gries_dp({"a": 5, "b": 5}, {"a": 50, "b": 50}), 3)
    0.0
    """
    tot_s = float(sum(sizes.values()))
    tot_f = float(sum(form_counts.values()))
    if not tot_s or not tot_f:
        return float("nan")
    dp = 0.0
    for part, n in sizes.items():
        dp += abs(form_counts.get(part, 0) / tot_f - n / tot_s)
    return 0.5 * dp


def bigram_model(forms: Mapping[Units, float], vocab: Sequence, alpha: float = 0.5) -> dict:
    """Unit-bigram model (with start ``<`` and end ``>`` symbols) from weighted forms, add-``alpha`` smoothing.

    Returns ``{"counts": {prev: Counter}, "vocab": [...], "alpha": alpha}`` for :func:`bigram_loglik`.
    """
    counts: dict = {}
    for f, w in forms.items():
        seq = ("<",) + tuple(f) + (">",)
        for a, b in zip(seq[:-1], seq[1:]):
            counts.setdefault(a, Counter())[b] += w
    return {"counts": counts, "vocab": list(vocab) + [">"], "alpha": float(alpha)}


def bigram_loglik(form: Units, model: dict) -> float:
    """Mean log2 probability per transition of a form under :func:`bigram_model`."""
    seq = ("<",) + tuple(form) + (">",)
    V = len(model["vocab"])
    a = model["alpha"]
    ll = 0.0
    for p, n in zip(seq[:-1], seq[1:]):
        c = model["counts"].get(p, {})
        tot = float(sum(c.values()))
        ll += np.log2((c.get(n, 0) + a) / (tot + a * V))
    return float(ll / (len(seq) - 1))


def form_profile(
    forms: Sequence[Units], tok: Mapping[str, np.ndarray], label_counts: Mapping[Units, tuple] | None = None,
    model: dict | None = None, rng: np.random.Generator | None = None, n_shuf: int = 20,
    seam_ref: np.ndarray | None = None,
) -> dict:
    """Profile of a set of forms inside one token sample.

    ``tok`` holds, for every token of the sample, arrays ``form`` (tuple), ``stem``, ``suffix``,
    ``page``, ``section``, ``line`` (line id), ``is_first``, ``is_last``, ``is_par_first`` (bool) and
    ``next_first`` (first unit of the next token in the line, ``""`` at the line end).
    ``label_counts[form] = (n_in_labels, n_in_paragraphs)`` from the full text (label share);
    ``model`` = :func:`bigram_model` of the shared forms (log-likelihood); ``seam_ref`` = boolean mask
    of the reference tokens (the shared forms) for the seam comparison.  Optional ``form_code``
    (integer code per token) and ``form_uniq`` (code -> form) speed up repeated calls.

    Metrics (type means over forms unless noted): ``paradigm_forms`` = distinct sample forms sharing
    the form's stem; ``H_suffix_given_stem`` = conditional entropy over the sample tokens whose stem
    is one of the forms' stems (token-weighted, bits); ``pages_per_token``, ``gries_dp``,
    ``top_page_share``; ``section_entropy`` (bits); ``label_share`` (tokens in labels / all, pooled),
    ``label_any`` (share of forms with a label occurrence); ``first_share``, ``last_share``,
    ``par_first_share`` (pooled over tokens); ``cooc_ratio`` = within-line pairs of forms' tokens /
    mean of ``n_shuf`` shuffles of the tokens' membership; ``bigram_ll`` (mean per transition);
    ``seam_jsd`` = JSD (bits) between the next-first-unit distribution after the forms' tokens and after
    the reference tokens; ``n_forms``, ``n_tokens``.
    """
    fset = set(forms)
    out = {"n_forms": len(fset), "n_tokens": 0}
    if not fset:
        return out
    codes, uniq = tok.get("form_code"), tok.get("form_uniq")
    if codes is None:
        codes, uniq = pd.factorize(pd.Series(tok["form"], dtype=object))
        uniq = list(uniq)
    code_of = {f: i for i, f in enumerate(uniq)}
    sel_codes = np.asarray([code_of[f] for f in fset if f in code_of], dtype=np.int64)
    mask = np.isin(codes, sel_codes)
    n = int(mask.sum())
    out["n_tokens"] = n
    if n == 0:
        return out
    # paradigm: distinct sample forms per stem, H(suffix | stem) over the tokens of the forms' stems
    stem_code = pd.Series(tok["stem"], dtype=object)
    stem_of_code = stem_code.groupby(codes).first()
    forms_per_stem = stem_of_code.groupby(stem_of_code.values).size()
    out["paradigm_forms"] = float(forms_per_stem.reindex(stem_of_code.loc[sel_codes].to_numpy()).mean())
    stems = set(stem_of_code.loc[sel_codes].tolist())
    smask = stem_code.isin(stems).to_numpy()
    if smask.any():
        d = pd.DataFrame({"s": tok["stem"][smask], "a": tok["suffix"][smask]})
        tot = len(d)
        h = 0.0
        for _, g in d.groupby("s"):
            h += len(g) / tot * entropy_bits(g["a"].value_counts().to_numpy())
        out["H_suffix_given_stem"] = h
    else:
        out["H_suffix_given_stem"] = np.nan
    # dispersion over pages and sections (type means)
    sizes = pd.Series(tok["page"], dtype=object).value_counts()
    d = pd.DataFrame({"f": codes[mask], "p": tok["page"][mask], "s": tok["section"][mask]})
    nf = d.groupby("f").size()
    fp = d.groupby(["f", "p"]).size()
    npages = fp.groupby(level=0).size()
    top = fp.groupby(level=0).max()
    out["pages_per_token"] = float((npages / nf).mean())
    out["top_page_share"] = float((top / nf).mean())
    share_f = fp / nf.reindex(fp.index.get_level_values(0)).to_numpy()
    share_p = (sizes / sizes.sum()).reindex(fp.index.get_level_values(1)).to_numpy()
    # DP = 0.5 * (sum over pages with the form |v - s| + sum over the other pages s)
    dev = pd.Series(np.abs(share_f.to_numpy() - share_p) - share_p, index=fp.index).groupby(level=0).sum()
    out["gries_dp"] = float((0.5 * (dev + 1.0)).mean())
    fs = d.groupby(["f", "s"]).size()
    out["section_entropy"] = float(fs.groupby(level=0).apply(lambda v: entropy_bits(v.to_numpy())).mean())
    if label_counts is not None:
        lab = [label_counts.get(f, (0, 0)) for f in fset]
        tl, tp = float(sum(a for a, _ in lab)), float(sum(b for _, b in lab))
        out["label_share"] = tl / (tl + tp) if tl + tp else np.nan
        out["label_any"] = float(np.mean([a > 0 for a, _ in lab]))
    for name in ("is_first", "is_last", "is_par_first"):
        out[name.replace("is_", "") + "_share"] = float(tok[name][mask].mean())
    out["cooc_ratio"] = _cooc_ratio(mask, tok["line"], rng, n_shuf)
    if model is not None:
        out["bigram_ll"] = float(np.mean([bigram_loglik(f, model) for f in fset]))
    if seam_ref is not None:
        nxt = tok["next_first"]
        inner = nxt != ""
        a = Counter(nxt[mask & inner].tolist())
        b = Counter(nxt[seam_ref & inner].tolist())
        keys = sorted(set(a) | set(b))
        from .registers import jsd
        out["seam_jsd"] = jsd([a.get(k, 0) for k in keys], [b.get(k, 0) for k in keys]) if keys else np.nan
    return out


def _cooc_ratio(mask: np.ndarray, lines: np.ndarray, rng, n_shuf: int) -> float:
    def pairs(m):
        c = Counter(lines[m].tolist())
        return sum(k * (k - 1) / 2 for k in c.values())
    obs = pairs(mask)
    if obs == 0 or rng is None:
        return np.nan
    exp = []
    for _ in range(n_shuf):
        exp.append(pairs(rng.permutation(mask)))
    e = float(np.mean(exp))
    return obs / e if e else np.nan


# --------------------------------------------------------------------------- #
# the -ey / -edy choice
# --------------------------------------------------------------------------- #


def suffix_pairs(tokens: Sequence[str], s1: str, s2: str, exclude_base_end: str | None = None) -> pd.DataFrame:
    """Tokens ending in ``s1`` or ``s2`` (``s2`` checked first, so ``-edy`` beats ``-ey``), with base and choice.

    Returns a DataFrame indexed like ``tokens`` (positions of the matching tokens): ``base`` = token
    without the suffix, ``choice`` = 0 (``s1``) / 1 (``s2``).  Tokens whose base ends in
    ``exclude_base_end`` are dropped (``-y`` / ``-dy`` without the ``-ey`` / ``-edy`` pairs); a base
    must be non-empty.

    >>> suffix_pairs(["chey", "chedy", "dy", "ol"], "ey", "edy")["choice"].tolist()
    [0, 1]
    """
    idx, base, choice = [], [], []
    for i, t in enumerate(tokens):
        t = str(t)
        if t.endswith(s2):
            if len(t) <= len(s2):
                continue
            b, c = t[: -len(s2)], 1
        elif t.endswith(s1) and len(t) > len(s1):
            b, c = t[: -len(s1)], 0
        else:
            continue
        if exclude_base_end and b.endswith(exclude_base_end):
            continue
        idx.append(i)
        base.append(b)
        choice.append(c)
    return pd.DataFrame({"base": base, "choice": choice}, index=idx)


def pair_choice(choice: np.ndarray, context: np.ndarray, min_n: int = 5, n_perm: int = 200, seed: int = 0) -> dict:
    """How far does a context determine the choice between two forms?

    ``choice`` = 0 / 1 per token, ``context`` = the context value per token.  Returns ``p_global``
    (share of choice 1), ``H`` (entropy of the choice, bits), ``H_given`` (conditional entropy given
    the context), ``reduction`` = H - H_given, ``perm_mean`` / ``perm_sd`` / ``z`` / ``p`` of the
    reduction under ``n_perm`` permutations of the choices among the tokens, ``n_contexts``
    (with >= ``min_n`` tokens), ``share_determined`` = share of tokens in contexts with P in
    [0, 0.1] or [0.9, 1], ``share_middle`` = share in contexts with P in [0.3, 0.7] (contexts with
    >= ``min_n`` tokens), ``n``.
    """
    ch = np.asarray(choice, dtype=int)
    ctx = np.asarray(context, dtype=object)
    n = len(ch)
    out: dict = {"n": int(n)}
    if n == 0:
        return {**out, "p_global": np.nan, "H": np.nan, "H_given": np.nan, "reduction": np.nan, "perm_mean": np.nan,
                "perm_sd": np.nan, "z": np.nan, "p": np.nan, "n_contexts": 0, "share_determined": np.nan,
                "share_middle": np.nan}
    codes, inv = np.unique(ctx, return_inverse=True)

    def h_given(c):
        tot = np.bincount(inv, minlength=len(codes)).astype(float)
        ones = np.bincount(inv, weights=c, minlength=len(codes))
        p = np.divide(ones, tot, out=np.zeros_like(tot), where=tot > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            h = -(np.where(p > 0, p * np.log2(p), 0) + np.where(p < 1, (1 - p) * np.log2(1 - p), 0))
        return float((tot / n * h).sum()), tot, p
    hg, tot, p = h_given(ch)
    pg = float(ch.mean())
    H = entropy_bits([pg * n, (1 - pg) * n])
    rng = np.random.default_rng(seed)
    sims = np.asarray([H - h_given(rng.permutation(ch))[0] for _ in range(n_perm)]) if n_perm else np.asarray([])
    red = H - hg
    sd = float(sims.std()) if sims.size else np.nan
    big = tot >= min_n
    w = tot[big] / tot[big].sum() if big.any() else np.asarray([])
    det = float((w * ((p[big] <= 0.1) | (p[big] >= 0.9))).sum()) if big.any() else np.nan
    mid = float((w * ((p[big] >= 0.3) & (p[big] <= 0.7))).sum()) if big.any() else np.nan
    out.update({"p_global": pg, "H": H, "H_given": hg, "reduction": red,
                "perm_mean": float(sims.mean()) if sims.size else np.nan, "perm_sd": sd,
                "z": (red - sims.mean()) / sd if sims.size and sd > 0 else np.nan,
                "p": float((1 + (sims >= red - 1e-12).sum()) / (1 + len(sims))) if sims.size else np.nan,
                "n_contexts": int(big.sum()), "share_determined": det, "share_middle": mid})
    return out


def both_forms(base: Sequence[str], choice: Sequence[int], page: Sequence, min_n: int = 2) -> dict:
    """Bases (stems) that occur with both alternative forms, and with both on one page.

    Over bases with >= ``min_n`` tokens: ``n_bases``, ``share_both`` (bases with both choices),
    ``share_both_same_page`` (of those, bases with both choices on at least one page),
    ``token_share_both`` (tokens of bases with both forms / tokens of all bases).

    >>> both_forms(["a", "a", "b", "b"], [0, 1, 0, 0], ["p1", "p1", "p1", "p2"])["share_both"]
    0.5
    """
    d = pd.DataFrame({"b": list(base), "c": list(choice), "p": list(page)})
    g = d.groupby("b")
    n = g.size()
    keep = n.index[n >= min_n]
    if not len(keep):
        return {"n_bases": 0, "share_both": np.nan, "share_both_same_page": np.nan, "token_share_both": np.nan}
    d = d[d["b"].isin(keep)]
    both = d.groupby("b")["c"].nunique() == 2
    bb = both.index[both]
    same = 0
    for b in bb:
        x = d[d["b"] == b].groupby("p")["c"].nunique()
        if (x == 2).any():
            same += 1
    return {"n_bases": int(len(keep)), "share_both": float(both.mean()),
            "share_both_same_page": same / len(bb) if len(bb) else np.nan,
            "token_share_both": float(d["b"].isin(bb).mean())}
