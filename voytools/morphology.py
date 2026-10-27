"""Word-structure analysis: affixes, prefix-stem-suffix decomposition, paradigms and BPE.

Voynichese words show a strong slot structure (Stolfi's "crust-mantle-core" model;
Zattera 2022 gives a glyph-level slot grammar).  This module offers simple, transparent
tools in that spirit; none of them claims to identify true morphemes.

* :func:`affix_inventory` lists initial/final glyph n-grams with their productivity
  (how many different remainders they combine with).
* :func:`decompose` / :func:`decompose_table` split a word greedily into
  ``prefix + stem + suffix`` (longest match first, prefix before suffix, stem never empty)
  using an affix list that is either given or learned from the data
  (:func:`learn_affixes`: the most productive initial/final n-grams).
* :func:`paradigm_table` and :func:`paradigm_metrics` summarise how freely stems combine with
  suffixes and prefixes, including a mutual-information measure with a permutation baseline.
* :func:`learn_bpe`, :func:`apply_bpe` and :func:`bpe_dependence_curve` implement
  within-token byte-pair encoding and the "dependence gap" ``H(U) - H(U_{t+1}|U_t)`` used by
  Rozanova & Temerev (2026): merging frequent glyph pairs into units makes the unit stream
  more or less predictable, and the point where the gap is smallest is read as the "natural
  unit scale".  Merges never cross token boundaries.

Default affix lists
-------------------
``decompose_table(df)`` with no lists learns them from ``df`` (the 12 most productive
prefixes and 14 most productive suffixes of up to 4 glyphs, occurring at least 50 times).
For ZL3b paragraph text (basic preset) this yields the lists recorded in
:data:`REFERENCE_PREFIXES` and :data:`REFERENCE_SUFFIXES`; pass ``prefixes='reference'``
to use those fixed lists, which is required for fair comparisons between corpora (learn the
lists on the real text and pass them explicitly when decomposing a null corpus).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from ._common import glyph_lists, iter_groups, resolve_preset
from .glyphs import collapse_eva
from .stats import conditional_entropy, entropy, jsd

__all__ = [
    "REFERENCE_PREFIXES",
    "REFERENCE_SUFFIXES",
    "NULL_LABEL",
    "affix_inventory",
    "learn_affixes",
    "decompose",
    "decompose_table",
    "paradigm_table",
    "ParadigmMetrics",
    "paradigm_metrics",
    "inventory_jaccard",
    "productivity_slope",
    "stem_dispersion",
    "slot_independent_sample",
    "learn_bpe",
    "apply_bpe",
    "bpe_dependence_curve",
]

NULL_LABEL = "-"
"""Label used for an empty prefix/suffix in printed tables."""

# The lists learned by decompose_table() from the ZL3b paragraph text (generic type P, clean tokens,
# 'basic' glyph preset); used for prefixes='reference' / suffixes='reference'.
REFERENCE_PREFIXES: tuple[str, ...] = ("che", "cho", "qo", "ch", "ol", "sh", "d", "l", "o", "s", "t", "y")
REFERENCE_SUFFIXES: tuple[str, ...] = (
    "aiin", "ain", "chy", "dy", "al", "am", "ar", "ol", "or", "d", "l", "r", "s", "y",
)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _units(df: pd.DataFrame, preset: Any, unit: str) -> list[tuple[str, ...]]:
    if unit == "glyph":
        return glyph_lists(df, preset)
    if unit == "char":
        return [tuple(t) for t in df["token"].tolist()]
    raise ValueError("unit must be 'glyph' or 'char'")


# --------------------------------------------------------------------------- #
# affix inventory
# --------------------------------------------------------------------------- #


def affix_inventory(
    df: pd.DataFrame,
    side: str = "prefix",
    max_len: int = 3,
    min_count: int = 20,
    preset: Any = "auto",
    unit: str = "glyph",
    by: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Inventory of initial (``side='prefix'``) or final (``'suffix'``) glyph n-grams.

    Columns
    -------
    ``affix``, ``n_glyphs``
        the n-gram and its length in glyphs;
    ``count``, ``share``
        tokens that begin/end with it, and their share of all tokens;
    ``n_types``
        distinct token types with the affix;
    ``n_remainders``
        distinct non-empty remainders after removing it (the productivity measure asked
        for; note that it is almost the same as ``n_types`` because different types leave
        different remainders);
    ``n_free``
        remainders that also occur as complete tokens ("attested stems") - a stricter,
        more informative productivity measure, used by :func:`learn_affixes`;
    ``exact``
        tokens that consist of the affix alone;
    ``ext_share``
        share of the affix's tokens that continue with the single most frequent
        one-glyph-longer affix (close to 1 means the affix is not an independent unit,
        e.g. ``q`` is nearly always ``qo``).

    Rows with ``count < min_count`` are dropped; sorted by count.  With ``by`` extra
    leading columns hold the group.
    """
    if side not in ("prefix", "suffix"):
        raise ValueError("side must be 'prefix' or 'suffix'")
    preset = resolve_preset(df, preset)
    frames = []
    cols = [] if by is None else ([by] if isinstance(by, str) else list(by))
    for key, sub in iter_groups(df, by):
        units = _units(sub, preset, unit)
        freq: Counter = Counter(units)
        total = sum(freq.values())
        free = {"".join(g) for g in freq}
        count: Counter = Counter()
        types: Counter = Counter()
        exact: Counter = Counter()
        rem: dict[tuple, set] = {}
        free_rem: dict[tuple, set] = {}
        ext: dict[tuple, Counter] = {}
        for g, c in freq.items():
            L = len(g)
            for n in range(1, min(max_len + 1, L) + 1):
                if side == "prefix":
                    aff, r = g[:n], g[n:]
                else:
                    aff, r = g[-n:], g[:-n]
                if n <= max_len:
                    count[aff] += c
                    types[aff] += 1
                    if r:
                        rem.setdefault(aff, set()).add(r)
                        if "".join(r) in free:
                            free_rem.setdefault(aff, set()).add(r)
                    else:
                        exact[aff] += c
                if n >= 2:  # aff is a one-glyph extension of the (n-1)-glyph affix
                    shorter = aff[: n - 1] if side == "prefix" else aff[1:]
                    ext.setdefault(shorter, Counter())[aff] += c
        rows = []
        for a, c in count.items():
            if c < min_count:
                continue
            top = ext.get(a)
            rows.append(
                {
                    "affix": "".join(a),
                    "n_glyphs": len(a),
                    "count": c,
                    "share": c / total if total else float("nan"),
                    "n_types": types[a],
                    "n_remainders": len(rem.get(a, ())),
                    "n_free": len(free_rem.get(a, ())),
                    "exact": exact.get(a, 0),
                    "ext_share": max(top.values()) / c if top else 0.0,
                }
            )
        if not rows:
            continue
        f = pd.DataFrame(rows).sort_values(["count", "affix"], ascending=[False, True], kind="stable")
        for c_name, k in zip(reversed(cols), reversed(key)):
            f.insert(0, c_name, k)
        frames.append(f)
    if not frames:
        return pd.DataFrame(
            columns=[*cols, "affix", "n_glyphs", "count", "share", "n_types", "n_remainders", "n_free", "exact", "ext_share"]
        )
    return pd.concat(frames, ignore_index=True)


def learn_affixes(
    df: pd.DataFrame,
    side: str = "prefix",
    top_n: int = 12,
    max_len: int = 4,
    min_count: int = 50,
    preset: Any = "auto",
    max_dominance: float = 0.9,
) -> list[str]:
    """Learn a list of productive affixes from the data.

    Candidates are the initial (final) glyph n-grams of up to ``max_len`` glyphs that occur at
    least ``min_count`` times.  They are ranked by ``n_free`` (number of attested stems, see
    :func:`affix_inventory`), then by count.  A candidate is skipped if one longer
    candidate (still within ``max_len``) accounts for at least ``max_dominance`` of its
    tokens (so ``q`` gives way to ``qo``, and ``iin`` to ``aiin``).  The ``top_n`` survivors are
    returned sorted from longest to shortest, ready for greedy longest-match decomposition.
    """
    inv = affix_inventory(df, side, max_len, min_count, preset)
    if not len(inv):
        return []
    inv = inv.sort_values(["n_free", "count", "affix"], ascending=[False, False, True], kind="stable")
    dominated = (inv["ext_share"] >= max_dominance) & (inv["n_glyphs"] < max_len)
    chosen = inv.loc[~dominated.to_numpy(dtype=bool), "affix"].head(top_n).tolist()
    p = resolve_preset(df, preset)
    return sorted(chosen, key=lambda a: (-len(collapse_eva(a, p)), a))


# --------------------------------------------------------------------------- #
# decomposition
# --------------------------------------------------------------------------- #


def _prep_affixes(affixes: Iterable[str | Sequence[str]], preset: Any) -> list[tuple[str, ...]]:
    out = []
    for a in affixes:
        g = collapse_eva(a, preset) if isinstance(a, str) else tuple(a)
        if g:
            out.append(g)
    return sorted(set(out), key=lambda g: (-len(g), g))


def _decompose(g: tuple[str, ...], pre: list[tuple[str, ...]], suf: list[tuple[str, ...]], min_stem: int):
    p: tuple[str, ...] = ()
    for a in pre:
        if len(g) - len(a) >= min_stem and g[: len(a)] == a:
            p = a
            break
    rest = g[len(p):]
    s: tuple[str, ...] = ()
    for a in suf:
        if len(rest) - len(a) >= min_stem and rest[-len(a):] == a:
            s = a
            break
    stem = rest[: len(rest) - len(s)] if s else rest
    return p, stem, s


def decompose(
    token_glyphs: Sequence[str] | str,
    prefixes: Iterable[str | Sequence[str]],
    suffixes: Iterable[str | Sequence[str]],
    preset: Any = "basic",
    min_stem: int = 1,
) -> tuple[str, str, str]:
    """Greedy longest-match split of one word into ``(prefix, stem, suffix)``.

    The prefix is chosen first (longest matching prefix that leaves at least ``min_stem``
    glyphs), then the suffix from what remains.  Affixes given as strings are converted to
    glyphs with ``preset`` so that ``'ch'`` is one glyph; ``token_glyphs`` may be a glyph
    sequence or a plain string (collapsed with ``preset``).  Empty prefix/suffix are ``''``.

    >>> decompose("qokeedy", ["qo", "o"], ["edy", "dy", "y"])
    ('qo', 'ke', 'edy')
    >>> decompose("daiin", ["qo"], ["aiin", "iin"])
    ('', 'd', 'aiin')
    >>> decompose("or", ["o"], ["or", "r"])
    ('o', 'r', '')
    """
    g = collapse_eva(token_glyphs, preset) if isinstance(token_glyphs, str) else tuple(token_glyphs)
    p, stem, s = _decompose(g, _prep_affixes(prefixes, preset), _prep_affixes(suffixes, preset), min_stem)
    return "".join(p), "".join(stem), "".join(s)


def decompose_table(
    df: pd.DataFrame,
    prefixes: Iterable[str] | str | None = None,
    suffixes: Iterable[str] | str | None = None,
    preset: Any = "auto",
    *,
    n_prefixes: int = 12,
    n_suffixes: int = 14,
    max_len: int = 4,
    min_count: int = 50,
    min_stem: int = 1,
    learn_on: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return a copy of the token table with the columns ``prefix``, ``stem``, ``suffix``.

    Parameters
    ----------
    prefixes, suffixes:
        Explicit affix lists (Eva strings), ``'reference'`` for the built-in lists
        (:data:`REFERENCE_PREFIXES`, :data:`REFERENCE_SUFFIXES`), or ``None`` (default) to
        learn them with :func:`learn_affixes` from ``learn_on`` (default: ``df`` itself).
    n_prefixes, n_suffixes, max_len, min_count:
        Parameters of the learning step.
    min_stem:
        Minimum number of glyphs left for the stem.

    The lists actually used are stored in ``result.attrs['prefixes']`` / ``['suffixes']``
    so that a null corpus can be decomposed with exactly the same lists.
    """
    preset = resolve_preset(df, preset)
    src = df if learn_on is None else learn_on

    def resolve(spec, side, n):
        if spec is None:
            return learn_affixes(src, side, n, max_len, min_count, preset)
        if isinstance(spec, str):
            if spec == "reference":
                return list(REFERENCE_PREFIXES if side == "prefix" else REFERENCE_SUFFIXES)
            raise ValueError("affix spec must be None, 'reference' or a list of strings")
        return list(spec)

    pre_list = resolve(prefixes, "prefix", n_prefixes)
    suf_list = resolve(suffixes, "suffix", n_suffixes)
    pre = _prep_affixes(pre_list, preset)
    suf = _prep_affixes(suf_list, preset)
    glyphs = glyph_lists(df, preset)
    cache: dict[tuple[str, ...], tuple[str, str, str]] = {}
    P, S, X = [], [], []
    for g in glyphs:
        r = cache.get(g)
        if r is None:
            p, stem, s = _decompose(g, pre, suf, min_stem)
            r = cache[g] = ("".join(p), "".join(stem), "".join(s))
        P.append(r[0])
        S.append(r[1])
        X.append(r[2])
    out = df.copy()
    out["prefix"] = P
    out["stem"] = S
    out["suffix"] = X
    out.attrs.update(df.attrs)
    out.attrs["prefixes"] = pre_list
    out.attrs["suffixes"] = suf_list
    return out


# --------------------------------------------------------------------------- #
# paradigms
# --------------------------------------------------------------------------- #


def paradigm_table(
    dfd: pd.DataFrame,
    counts: str = "tokens",
    min_stem_count: int = 1,
    top_stems: int | None = None,
    top_suffixes: int | None = None,
    axis: str = "suffix",
) -> pd.DataFrame:
    """Contingency table of stems (rows) by suffixes (columns) from a decomposed table.

    Parameters
    ----------
    counts:
        ``'tokens'`` counts occurrences, ``'types'`` counts distinct word forms.
    min_stem_count, top_stems, top_suffixes:
        Keep only stems with at least this many tokens, the ``top_stems`` most frequent
        stems and the ``top_suffixes`` most frequent suffixes.
    axis:
        ``'suffix'`` (default) or ``'prefix'`` for the column variable.

    The empty affix is labelled :data:`NULL_LABEL` (``'-'``).
    """
    if axis not in ("suffix", "prefix"):
        raise ValueError("axis must be 'suffix' or 'prefix'")
    if counts not in ("tokens", "types"):
        raise ValueError("counts must be 'tokens' or 'types'")
    stems = dfd["stem"].tolist()
    cols = [a or NULL_LABEL for a in dfd[axis].tolist()]
    if counts == "types":
        seen = set()
        pairs = []
        for s, c, t in zip(stems, cols, dfd["token"].tolist()):
            if (s, c, t) not in seen:
                seen.add((s, c, t))
                pairs.append((s, c))
    else:
        pairs = list(zip(stems, cols))
    if not pairs:
        return pd.DataFrame()
    stem_names, stem_codes = np.unique(np.asarray([p[0] for p in pairs], dtype=object), return_inverse=True)
    col_names, col_codes = np.unique(np.asarray([p[1] for p in pairs], dtype=object), return_inverse=True)
    mat = np.zeros((len(stem_names), len(col_names)), dtype=np.int64)
    np.add.at(mat, (stem_codes, col_codes), 1)
    table = pd.DataFrame(mat, index=stem_names.tolist(), columns=col_names.tolist())
    table.index.name = "stem"
    table = table[table.sum(axis=1) >= min_stem_count]
    order = table.sum(axis=1).sort_values(ascending=False, kind="stable").index
    table = table.loc[order]
    if top_stems is not None:
        table = table.head(top_stems)
    colorder = table.sum(axis=0).sort_values(ascending=False, kind="stable").index
    table = table[list(colorder)]
    if top_suffixes is not None:
        table = table.iloc[:, :top_suffixes]
    return table


@dataclass
class ParadigmMetrics:
    """Result of :func:`paradigm_metrics`.

    ``per_stem`` has one row per stem with at least ``min_stem_count`` tokens
    (``count``, ``n_types``, ``n_suffixes``, ``n_prefixes``, ``H_suffix``, ``H_prefix``);
    ``summary`` holds the global numbers (see :func:`paradigm_metrics`).
    """

    per_stem: pd.DataFrame
    summary: dict[str, float]


def _codes(values: Sequence[str]) -> tuple[np.ndarray, int]:
    uniq, inv = np.unique(np.asarray(list(values), dtype=object), return_inverse=True)
    return inv.astype(np.int64), len(uniq)


def _h_keys(*arrays: np.ndarray, sizes: Sequence[int]) -> float:
    key = np.zeros(len(arrays[0]), dtype=np.int64)
    for a, k in zip(arrays, sizes):
        key = key * k + a
    _, c = np.unique(key, return_counts=True)
    return entropy(c)


def _cmi(p: np.ndarray, s: np.ndarray, z: np.ndarray, kp: int, ks: int, kz: int) -> float:
    """Plug-in conditional mutual information I(P;S|Z) in bits."""
    return (
        _h_keys(p, z, sizes=(kp, kz))
        + _h_keys(s, z, sizes=(ks, kz))
        - _h_keys(p, s, z, sizes=(kp, ks, kz))
        - _h_keys(z, sizes=(kz,))
    )


def paradigm_metrics(
    dfd: pd.DataFrame,
    min_stem_count: int = 5,
    shuffles: int = 20,
    seed: int = 0,
) -> ParadigmMetrics:
    """Stem-level productivity and prefix-suffix dependence of a decomposed table.

    ``summary`` keys
    ----------------
    ``n_tokens``, ``n_stems`` (all stems), ``n_stems_k`` (stems with >= ``min_stem_count``
    tokens), ``token_coverage_k``, ``mean_suffixes_per_stem`` and ``mean_prefixes_per_stem``
    (over stems >= k, unweighted), ``frac_stems_multi_suffix`` (share of those stems seen with
    more than one suffix), ``H_suffix_given_stem`` (token-weighted mean over stems >= k),
    ``mi_prefix_suffix`` (unconditional I(P;S) over all tokens) and
    ``mi_prefix_suffix_shuffled`` (mean over ``shuffles`` permutations of the suffix labels),
    ``cmi_prefix_suffix_given_stem`` (token-weighted I(P;S|stem) over stems >= k),
    ``cmi_shuffled`` (prefix labels permuted within stems) and ``cmi_excess``.

    The shuffled values estimate the upward bias of the plug-in estimators; the excess is
    the part attributable to genuine prefix-suffix dependence given the stem.
    Note that these are properties of the *vocabulary with frequencies*; they do not change when
    tokens are merely reordered (e.g. a within-line shuffle).
    """
    stem = dfd["stem"].tolist()
    pre = dfd["prefix"].tolist()
    suf = dfd["suffix"].tolist()
    tok = dfd["token"].tolist()
    n = len(stem)
    zc, kz = _codes(stem)
    pc, kp = _codes(pre)
    sc, ks = _codes(suf)
    counts = np.bincount(zc, minlength=kz)
    keep_stem = counts >= min_stem_count
    keep = keep_stem[zc]

    # per stem
    rows = []
    stems_u = np.unique(np.asarray(stem, dtype=object))
    by_stem: dict[int, list[int]] = {}
    for i in np.flatnonzero(keep):
        by_stem.setdefault(int(zc[i]), []).append(int(i))
    for z, idxs in by_stem.items():
        ii = np.asarray(idxs)
        sc_ = np.bincount(sc[ii])
        pc_ = np.bincount(pc[ii])
        rows.append(
            {
                "stem": stems_u[z],
                "count": len(ii),
                "n_types": len({tok[i] for i in idxs}),
                "n_suffixes": int((sc_ > 0).sum()),
                "n_prefixes": int((pc_ > 0).sum()),
                "H_suffix": entropy(sc_),
                "H_prefix": entropy(pc_),
            }
        )
    per = pd.DataFrame(rows, columns=["stem", "count", "n_types", "n_suffixes", "n_prefixes", "H_suffix", "H_prefix"])
    if len(per):
        per = per.sort_values(["count", "stem"], ascending=[False, True], kind="stable").reset_index(drop=True)

    rng = np.random.default_rng(seed)
    mi_u = _h_keys(pc, sizes=(kp,)) + _h_keys(sc, sizes=(ks,)) - _h_keys(pc, sc, sizes=(kp, ks)) if n else float("nan")
    mi_u_sh = []
    for _ in range(shuffles):
        mi_u_sh.append(_h_keys(pc, sizes=(kp,)) + _h_keys(sc, sizes=(ks,)) - _h_keys(pc, rng.permutation(sc), sizes=(kp, ks)))
    zk, pk, sk = zc[keep], pc[keep], sc[keep]
    if len(zk):
        cmi = _cmi(pk, sk, zk, kp, ks, kz)
        base = np.argsort(zk, kind="stable")
        sh = []
        for _ in range(shuffles):
            order = np.lexsort((rng.random(len(zk)), zk))
            pnew = np.empty_like(pk)
            pnew[base] = pk[order]
            sh.append(_cmi(pnew, sk, zk, kp, ks, kz))
        cmi_sh = float(np.mean(sh)) if sh else float("nan")
    else:
        cmi, cmi_sh = float("nan"), float("nan")
    summary = {
        "n_tokens": n,
        "n_stems": int(len(stems_u)),
        "n_stems_k": int(len(per)),
        "token_coverage_k": float(keep.sum() / n) if n else float("nan"),
        "mean_suffixes_per_stem": float(per["n_suffixes"].mean()) if len(per) else float("nan"),
        "mean_prefixes_per_stem": float(per["n_prefixes"].mean()) if len(per) else float("nan"),
        "frac_stems_multi_suffix": float((per["n_suffixes"] > 1).mean()) if len(per) else float("nan"),
        "H_suffix_given_stem": float(np.average(per["H_suffix"], weights=per["count"])) if len(per) else float("nan"),
        "mi_prefix_suffix": float(mi_u),
        "mi_prefix_suffix_shuffled": float(np.mean(mi_u_sh)) if mi_u_sh else float("nan"),
        "cmi_prefix_suffix_given_stem": float(cmi),
        "cmi_shuffled": cmi_sh,
        "cmi_excess": float(cmi - cmi_sh) if shuffles else float("nan"),
    }
    return ParadigmMetrics(per, summary)


def inventory_jaccard(
    dfd: pd.DataFrame,
    by: str,
    column: str = "suffix",
    min_count: int = 3,
    groups: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Pairwise similarity of the inventories of ``column`` between groups of a decomposed table.

    For every group of ``by`` (e.g. ``'illustration'`` or ``'language'``; rows with a missing
    value are dropped) the inventory is the set of values of ``column`` (``'stem'``,
    ``'suffix'``, ``'prefix'`` ...; the empty affix counts as a value) that occur at least
    ``min_count`` times in that group.  Returns one row per unordered pair of groups with
    ``group_a``, ``group_b``, ``n_tokens_a/b``, ``n_inv_a/b`` (inventory sizes), ``n_shared``,
    ``jaccard`` (|A & B| / |A | B|, NaN if both are empty) and ``jsd`` (Jensen-Shannon
    divergence in bits of the full frequency distributions, no ``min_count`` threshold).
    ``groups`` restricts and orders the groups (default: all, sorted).
    """
    sub = dfd[dfd[by].notna()]
    values = sub[column].astype(str).to_numpy()
    keys = sub[by].astype(str).to_numpy()
    names = sorted(set(keys.tolist())) if groups is None else [g for g in groups if g in set(keys.tolist())]
    counts = {g: Counter(values[keys == g].tolist()) for g in names}
    support = sorted(set(values.tolist()))
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ca, cb = counts[a], counts[b]
            ia = {v for v, c in ca.items() if c >= min_count}
            ib = {v for v, c in cb.items() if c >= min_count}
            union = ia | ib
            rows.append(
                {
                    "group_a": a,
                    "group_b": b,
                    "n_tokens_a": int(sum(ca.values())),
                    "n_tokens_b": int(sum(cb.values())),
                    "n_inv_a": len(ia),
                    "n_inv_b": len(ib),
                    "n_shared": len(ia & ib),
                    "jaccard": len(ia & ib) / len(union) if union else float("nan"),
                    "jsd": jsd([ca.get(v, 0) for v in support], [cb.get(v, 0) for v in support]),
                }
            )
    cols = ["group_a", "group_b", "n_tokens_a", "n_tokens_b", "n_inv_a", "n_inv_b", "n_shared", "jaccard", "jsd"]
    return pd.DataFrame(rows, columns=cols)


def productivity_slope(dfd: pd.DataFrame, min_stem_count: int = 5, axis: str = "suffix") -> dict[str, float]:
    """OLS fit of ``log(number of distinct affixes of a stem)`` on ``log(stem frequency)``.

    Uses stems with at least ``min_stem_count`` tokens; ``axis`` is ``'suffix'`` or
    ``'prefix'``.  Returns ``slope``, ``intercept`` (natural logs), ``r2`` and ``n_stems``.
    A stem that takes affixes freely gains new ones as it gets more frequent (slope near the
    value of a random pairing); a closed paradigm saturates (smaller slope).  NaN when fewer
    than two distinct frequencies remain.
    """
    if axis not in ("suffix", "prefix"):
        raise ValueError("axis must be 'suffix' or 'prefix'")
    g = pd.DataFrame({"stem": dfd["stem"].to_numpy(), "a": dfd[axis].to_numpy()}).groupby("stem")["a"]
    freq = g.size()
    n_aff = g.nunique()
    keep = freq >= min_stem_count
    x = np.log(freq[keep].to_numpy(dtype=float))
    y = np.log(n_aff[keep].to_numpy(dtype=float))
    out = {"slope": float("nan"), "intercept": float("nan"), "r2": float("nan"), "n_stems": int(keep.sum())}
    if len(np.unique(x)) < 2:
        return out
    A = np.column_stack([x, np.ones_like(x)])
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - (slope * x + intercept)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    out.update(slope=float(slope), intercept=float(intercept),
               r2=1.0 - float((resid**2).sum()) / ss_tot if ss_tot > 0 else float("nan"))
    return out


def stem_dispersion(dfd: pd.DataFrame, min_stem_count: int = 5, per_stem: bool = False) -> Any:
    """How evenly the occurrences of a stem are spread through the text.

    For every stem with at least ``min_stem_count`` tokens, the median distance (in tokens,
    in the row order of ``dfd``, which must be text order) between consecutive occurrences
    is divided by ``N / k`` (``N`` tokens in the table, ``k`` occurrences of the stem), the
    gap expected if the occurrences were evenly spaced.  Returns the median of this ratio over
    stems (``per_stem=True``: a DataFrame with ``stem``, ``count``, ``median_gap``, ``ratio``).

    Evenly spaced occurrences give 1; occurrences at independent random positions give gaps
    that are roughly geometric, so the ratio is about ``ln 2 ~ 0.69``; occurrences clustered in
    bursts (local copying, topic runs) give much smaller values.
    """
    stems = dfd["stem"].to_numpy()
    n = len(stems)
    codes, k = _codes(stems) if n else (np.zeros(0, dtype=np.int64), 0)
    cnt = np.bincount(codes, minlength=k)
    names = np.unique(np.asarray(stems, dtype=object)) if n else np.zeros(0, dtype=object)
    order = np.argsort(codes, kind="stable")  # positions grouped by stem, ascending within stem
    starts = np.concatenate([[0], np.cumsum(cnt)])
    rows = []
    for z in np.flatnonzero(cnt >= max(min_stem_count, 2)):
        pos = order[starts[z] : starts[z + 1]]
        gap = float(np.median(np.diff(pos)))
        rows.append({"stem": names[z], "count": int(cnt[z]), "median_gap": gap, "ratio": gap / (n / cnt[z])})
    per = pd.DataFrame(rows, columns=["stem", "count", "median_gap", "ratio"])
    if per_stem:
        return per.sort_values(["count", "stem"], ascending=[False, True], kind="stable").reset_index(drop=True)
    return float(per["ratio"].median()) if len(per) else float("nan")


def slot_independent_sample(dfd: pd.DataFrame, seed: int = 0, preset: Any = "auto") -> pd.DataFrame:
    """Slot-independence null: prefix, stem and suffix of every token drawn independently.

    For each row of the decomposed table ``dfd`` a prefix, a stem and a suffix are drawn
    independently and with replacement from the real token-level marginals (i.e. each slot
    value with probability proportional to its frequency) and concatenated into ``token``.
    Token count, line layout and all positional metadata are kept.  The ``prefix`` /
    ``stem`` / ``suffix`` columns are removed so that the sample goes through the same
    :func:`decompose_table` call (with the real affix lists) as the real text; the greedy
    re-parse may split a concatenation differently from the draw.
    """
    preset = resolve_preset(dfd, preset)
    rng = np.random.default_rng(seed)
    n = len(dfd)
    parts = [dfd[c].astype(str).to_numpy()[rng.integers(0, n, size=n)] if n else np.zeros(0, dtype=object)
             for c in ("prefix", "stem", "suffix")]
    toks = [p + s + x for p, s, x in zip(*parts)]
    out = dfd.drop(columns=["prefix", "stem", "suffix"]).copy()
    out["token"] = toks
    if "token_raw" in out.columns:
        out["token_raw"] = toks
    if "len_chars" in out.columns:
        out["len_chars"] = [len(t) for t in toks]
    if "len_glyphs" in out.columns:
        out["len_glyphs"] = [len(collapse_eva(t, preset)) for t in toks]
    out.attrs.update(dfd.attrs)
    for key in ("prefixes", "suffixes"):
        out.attrs.pop(key, None)
    out.attrs["source"] = "slot_independent"
    return out



# --------------------------------------------------------------------------- #
# byte-pair encoding (within tokens)
# --------------------------------------------------------------------------- #


def _as_units(token: str | Sequence[str], preset: Any) -> tuple[str, ...]:
    if isinstance(token, str):
        return collapse_eva(token, preset)
    return tuple(token)


def _word_freq(tokens: Iterable[Any], preset: Any) -> Counter:
    items = tokens
    if is_nested_tokens(tokens):
        items = (t for line in tokens for t in line)
    return Counter(_as_units(t, preset) for t in items)


def is_nested_tokens(tokens: Any) -> bool:
    """True if ``tokens`` is a list of lines (lists) rather than a flat list of tokens."""
    if isinstance(tokens, (str, bytes)):
        return False
    for first in tokens:
        # a glyph tuple (tuple of str) is a token, a list/tuple of tokens/glyph-tuples is a line
        if isinstance(first, list):
            return True
        if isinstance(first, tuple):
            return len(first) > 0 and not isinstance(first[0], str)
        return False
    return False


def _bpe_train(
    word_freq: Counter,
    merges: int,
    min_pair_count: int,
    snapshots: Iterable[int] = (),
):
    segs: dict[tuple[str, ...], list[str]] = {w: list(w) for w in word_freq}
    pair_counts: Counter = Counter()
    where: dict[tuple[str, str], set] = {}
    for w, units in segs.items():
        f = word_freq[w]
        for pair in zip(units, units[1:]):
            pair_counts[pair] += f
            where.setdefault(pair, set()).add(w)
    done: list[tuple[str, str]] = []
    snaps: dict[int, dict[tuple[str, ...], tuple[str, ...]]] = {}
    want = set(snapshots)
    if 0 in want:
        snaps[0] = {w: tuple(u) for w, u in segs.items()}
    for step in range(1, merges + 1):
        best = None
        best_key = None
        for pair, c in pair_counts.items():
            if c >= min_pair_count and c > 0:
                k = (c, pair)
                if best_key is None or k > best_key:
                    best_key, best = k, pair
        if best is None:
            break
        done.append(best)
        left, right = best
        merged = left + right
        for w in list(where.get(best, ())):
            units = segs[w]
            f = word_freq[w]
            for pair in zip(units, units[1:]):
                pair_counts[pair] -= f
                s = where.get(pair)
                if s is not None:
                    s.discard(w)
            out: list[str] = []
            i = 0
            while i < len(units):
                if i + 1 < len(units) and units[i] == left and units[i + 1] == right:
                    out.append(merged)
                    i += 2
                else:
                    out.append(units[i])
                    i += 1
            segs[w] = out
            for pair in zip(out, out[1:]):
                pair_counts[pair] += f
                where.setdefault(pair, set()).add(w)
        for pair in [p for p, c in pair_counts.items() if c <= 0]:
            del pair_counts[pair]
        if step in want:
            snaps[step] = {w: tuple(u) for w, u in segs.items()}
    if len(done) not in snaps:
        snaps[len(done)] = {w: tuple(u) for w, u in segs.items()}
    return done, snaps, segs


def learn_bpe(
    token_list: Iterable[Any],
    merges: int = 64,
    preset: Any = None,
    min_pair_count: int = 2,
) -> list[tuple[str, str]]:
    """Learn within-token byte-pair-encoding merges.

    Parameters
    ----------
    token_list:
        Token strings (or glyph tuples), one entry per token *occurrence* (so that merge
        decisions are frequency-weighted), or a list of lines of such tokens.
    merges:
        Maximum number of merges.
    preset:
        How strings are split into initial units: ``None`` = single characters, ``'basic'`` /
        ``'extended'`` = Eva glyphs with composites collapsed.
    min_pair_count:
        Stop when the best pair occurs fewer times than this.

    Returns
    -------
    list of ``(left, right)`` in merge order.  The most frequent pair is merged first; ties
    go to the lexicographically larger pair (as in the Rozanova & Temerev scripts).  Pairs
    are only counted inside one token, so merges never cross word boundaries.
    """
    wf = _word_freq(list(token_list), preset)
    done, _, _ = _bpe_train(wf, merges, min_pair_count)
    return done


def _apply_merges(units: tuple[str, ...], merges: Sequence[tuple[str, str]]) -> tuple[str, ...]:
    cur = list(units)
    for left, right in merges:
        if len(cur) < 2:
            break
        out: list[str] = []
        i = 0
        changed = False
        while i < len(cur):
            if i + 1 < len(cur) and cur[i] == left and cur[i + 1] == right:
                out.append(left + right)
                i += 2
                changed = True
            else:
                out.append(cur[i])
                i += 1
        if changed:
            cur = out
    return tuple(cur)


def apply_bpe(
    tokens: Iterable[Any],
    merges: Sequence[tuple[str, str]],
    preset: Any = None,
    n_merges: int | None = None,
) -> list[tuple[str, ...]]:
    """Segment every token with the first ``n_merges`` (default: all) learned merges.

    Returns a list of unit tuples aligned with ``tokens`` (flat list of tokens; for a
    list of lines the result is flattened line by line).
    """
    use = list(merges if n_merges is None else merges[:n_merges])
    tokens = list(tokens)
    flat = [t for line in tokens for t in line] if is_nested_tokens(tokens) else tokens
    cache: dict[tuple[str, ...], tuple[str, ...]] = {}
    out = []
    for t in flat:
        u = _as_units(t, preset)
        r = cache.get(u)
        if r is None:
            r = cache[u] = _apply_merges(u, use)
        out.append(r)
    return out


def bpe_dependence_curve(
    tokens: Iterable[Any],
    merge_checkpoints: Sequence[int] = (0, 16, 32, 64, 128, 256),
    preset: Any = None,
    train: Iterable[Any] | None = None,
    min_pair_count: int = 2,
) -> pd.DataFrame:
    """Dependence gap ``H(U) - H(U_{t+1}|U_t)`` of the BPE-segmented text at several merge counts.

    The text is the concatenation of the units of the tokens of each line (spaces removed, as
    in Rozanova & Temerev); pass a list of lines (lists of tokens) to keep line boundaries,
    or a flat token list for one stream.  Merges are learned on ``train`` (default: the
    same tokens) and the first ``k`` of them are applied for every checkpoint ``k``.

    Columns: ``merges`` (actually applied; fewer than requested if the pair statistics ran
    out), ``n_units``, ``n_symbols``, ``mean_glyphs_per_unit``, ``H1`` = H(U),
    ``H2`` = H(U_{t+1}|U_t), ``gap`` = H1 - H2 and ``bits_per_glyph`` = H2 * n_units / n_glyphs.
    """
    cps = sorted(set(int(k) for k in merge_checkpoints))
    if not cps:
        raise ValueError("merge_checkpoints must not be empty")
    tokens = list(tokens)
    lines = [list(l) for l in tokens] if is_nested_tokens(tokens) else [tokens]
    line_units = [[_as_units(t, preset) for t in l] for l in lines]
    eval_wf = Counter(u for l in line_units for u in l)
    train_wf = eval_wf if train is None else _word_freq(list(train), preset)
    done, snaps, _ = _bpe_train(train_wf, max(cps), min_pair_count, snapshots=cps)
    n_glyphs = sum(len(u) for l in line_units for u in l)
    rows = []
    for k in cps:
        k_eff = min(k, len(done))
        if train is None:
            seg = snaps[k_eff]
            segmented = {u: seg[u] for u in eval_wf}
        else:
            merges_k = done[:k_eff]
            segmented = {u: _apply_merges(u, merges_k) for u in eval_wf}
        stream = [[unit for u in l for unit in segmented[u]] for l in line_units]
        flat = [x for s in stream for x in s]
        if flat:
            _, counts = np.unique(np.asarray(flat, dtype=object), return_counts=True)
        else:
            counts = np.zeros(0)
        h1 = entropy(counts)
        h2 = conditional_entropy(stream, 1) if flat else 0.0
        rows.append(
            {
                "merges": k_eff,
                "n_units": len(flat),
                "n_symbols": int(len(counts)),
                "mean_glyphs_per_unit": n_glyphs / len(flat) if flat else float("nan"),
                "H1": h1,
                "H2": h2,
                "gap": h1 - h2,
                "bits_per_glyph": h2 * len(flat) / n_glyphs if n_glyphs else float("nan"),
            }
        )
    return pd.DataFrame(rows)
