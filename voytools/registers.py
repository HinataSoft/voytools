"""Comparing two registers of the text (Currier A and B): shared vocabulary, affix mapping, gradient.

Experiment 012.  The functions take plain arrays (stems, affixes, page ids, labels) so that the same
code measures the real A/B split, halves of one language, label permutations and generated nulls:

* :func:`draw_tokens` - a random set of whole pages of a pool, cut to exactly ``n`` tokens
  (equal-size comparisons; vocabulary measures depend strongly on size);
* :func:`stem_sharing` - Jaccard similarity and token coverage of the stem inventories;
* :func:`transfer_matrix`, :func:`mapping_fit`, :func:`max_weight_matching` - the within-stem affix
  transfer matrix M[a, b] = sum over shared stems of p(a | stem, X) p(b | stem, Y), its sparseness,
  diagonal mass, the best one-to-one relabelling (Hungarian algorithm, no scipy) and the JSD of the
  affix distributions left after it;
* :func:`jsd_contributions` - per-value terms of the Jensen-Shannon divergence;
* :func:`page_scores`, :func:`spearman`, :func:`order_drift` - projection of every page on the A-B
  axis (leave-one-page-out log-ratio score) and drift along the page order.

Experiment 016 (are the stems found only in one language glyph variants of the other's stems?):

* :func:`levenshtein`, :func:`edit_script`, :func:`nearest_counterparts` - edit distance on unit
  tuples, one optimal edit script (``ins:ed@end``), nearest neighbour search within 2 edits;
* :func:`counterpart_rates` - X-only stems and the share of their tokens with a counterpart in Y;
* :func:`edit_scripts`, :func:`edit_inventory`, :func:`coverage_levels` - ranked inventory of the
  operations and how many of them explain 50 / 80 / 95 % of the tokens;
* :func:`substitution_sharing` - stem sharing after rewriting with the top-k operations;
* :func:`shuffle_units` - random-form null; :func:`profile_test` - do counterpart pairs share a
  section profile more than random pairs?
"""

from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# equal-size samples of pages
# --------------------------------------------------------------------------- #


def draw_tokens(
    page: np.ndarray, pool: Sequence, n: int, rng: np.random.Generator, exclude: Sequence = ()
) -> np.ndarray:
    """Boolean token mask: whole pages of ``pool`` (minus ``exclude``) in random order until ``n`` tokens.

    The last page drawn is cut to its first tokens so that exactly ``n`` tokens are selected (fewer if
    the pool is too small).  ``page`` holds the page id of every token (tokens of a page contiguous
    or not; the cut keeps the page's first tokens in table order).
    """
    page = np.asarray(page)
    ex = set(exclude)
    cand = [p for p in pool if p not in ex]
    order = [cand[i] for i in rng.permutation(len(cand))]
    mask = np.zeros(len(page), dtype=bool)
    got = 0
    for p in order:
        if got >= n:
            break
        idx = np.flatnonzero(page == p)
        take = idx[: n - got]
        mask[take] = True
        got += len(take)
    return mask


def split_pages(pages: Sequence, sizes: Sequence[int]) -> tuple[list, list]:
    """Contiguous halves of a page sequence by token mass (first pages until half the tokens, rest)."""
    sizes = np.asarray(sizes, dtype=float)
    cum = np.cumsum(sizes)
    k = int(np.searchsorted(cum, cum[-1] / 2)) + 1 if len(sizes) else 0
    pages = list(pages)
    return pages[:k], pages[k:]


# --------------------------------------------------------------------------- #
# vocabulary sharing
# --------------------------------------------------------------------------- #


def stem_sharing(stems_x: Sequence[str], stems_y: Sequence[str], min_count: int = 3) -> dict:
    """Overlap of two stem inventories (stems with >= ``min_count`` tokens on their side).

    ``jaccard`` = |X & Y| / |X | Y| of the inventories; ``cover_xy`` = share of X tokens whose stem is
    in the inventory of Y (and ``cover_yx`` vice versa); ``n_x``, ``n_y`` = inventory sizes.

    >>> stem_sharing(list("aaabbbccc"), list("aaabbbddd"))["jaccard"]
    0.5
    """
    cx = pd.Series(list(stems_x), dtype=object).value_counts()
    cy = pd.Series(list(stems_y), dtype=object).value_counts()
    ix = set(cx.index[cx >= min_count])
    iy = set(cy.index[cy >= min_count])
    union = ix | iy
    tx, ty = cx.sum(), cy.sum()
    return {"jaccard": len(ix & iy) / len(union) if union else float("nan"),
            "cover_xy": float(cx[cx.index.isin(list(iy))].sum() / tx) if tx else float("nan"),
            "cover_yx": float(cy[cy.index.isin(list(ix))].sum() / ty) if ty else float("nan"),
            "n_x": len(ix), "n_y": len(iy)}


# --------------------------------------------------------------------------- #
# affix mapping
# --------------------------------------------------------------------------- #


def transfer_matrix(
    stem_x: Sequence[str], aff_x: Sequence[str], stem_y: Sequence[str], aff_y: Sequence[str],
    min_count: int = 5, labels: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, int]:
    """Within-stem transfer M[a, b] = sum_stem p(a | stem, X) p(b | stem, Y) over stems with >= ``min_count`` tokens in both.

    Returns the (unnormalised, sum = number of stems) matrix with rows = affixes of X, columns =
    affixes of Y (both indexed by ``labels`` if given, else the union of affixes seen) and the number
    of stems used.
    """
    x = pd.DataFrame({"stem": list(stem_x), "aff": list(aff_x)})
    y = pd.DataFrame({"stem": list(stem_y), "aff": list(aff_y)})
    nx, ny = x["stem"].value_counts(), y["stem"].value_counts()
    shared = sorted(set(nx.index[nx >= min_count]) & set(ny.index[ny >= min_count]))
    if labels is None:
        labels = sorted(set(x["aff"]) | set(y["aff"]))
    labels = list(labels)
    if not shared:
        return pd.DataFrame(0.0, index=labels, columns=labels), 0
    px = pd.crosstab(x["stem"], x["aff"]).reindex(index=shared, columns=labels, fill_value=0).to_numpy(dtype=float, copy=True)
    py = pd.crosstab(y["stem"], y["aff"]).reindex(index=shared, columns=labels, fill_value=0).to_numpy(dtype=float, copy=True)
    px /= px.sum(axis=1, keepdims=True)
    py /= py.sum(axis=1, keepdims=True)
    return pd.DataFrame(px.T @ py, index=labels, columns=labels), len(shared)


def max_weight_matching(w: np.ndarray) -> np.ndarray:
    """One-to-one assignment of rows to columns maximising the summed weight (Hungarian algorithm).

    Square or rectangular ``w`` (rows <= columns after an internal transpose); returns ``col[i]`` for
    every row (``-1`` for rows left unmatched when there are more rows than columns).

    >>> max_weight_matching(np.array([[1, 5], [4, 1]])).tolist()
    [1, 0]
    """
    w = np.asarray(w, dtype=float)
    transpose = w.shape[0] > w.shape[1]
    a = -(w.T if transpose else w)
    n, m = a.shape
    INF = float("inf")
    u, v = np.zeros(n + 1), np.zeros(m + 1)
    p = np.zeros(m + 1, dtype=np.int64)   # p[j] = row (1-based) matched to column j
    way = np.zeros(m + 1, dtype=np.int64)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = np.full(m + 1, INF)
        used = np.zeros(m + 1, dtype=bool)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta, j1 = INF, -1
            for j in range(1, m + 1):
                if not used[j]:
                    cur = a[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j], way[j] = cur, j0
                    if minv[j] < delta:
                        delta, j1 = minv[j], j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    row_to_col = np.full(n, -1, dtype=np.int64)
    for j in range(1, m + 1):
        if p[j]:
            row_to_col[p[j] - 1] = j - 1
    if not transpose:
        return row_to_col
    out = np.full(w.shape[0], -1, dtype=np.int64)
    for col, row in enumerate(row_to_col):  # transposed: "rows" were the original columns
        if row >= 0:
            out[row] = col
    return out


def jsd(p: np.ndarray, q: np.ndarray) -> float:
    """Jensen-Shannon divergence (bits) of two count or probability vectors."""
    return float(jsd_contributions(p, q).sum())


def jsd_contributions(p, q) -> np.ndarray:
    """Per-value terms of JSD(p, q) in bits (they sum to the divergence); inputs are counts or probabilities."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum() if p.sum() else p
    q = q / q.sum() if q.sum() else q
    m = (p + q) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(p > 0, p * np.log2(p / m), 0.0)
        b = np.where(q > 0, q * np.log2(q / m), 0.0)
    return (a + b) / 2


def mapping_fit(M: pd.DataFrame, count_x: pd.Series, count_y: pd.Series) -> dict:
    """Sparseness and fit of a transfer matrix (:func:`transfer_matrix`).

    ``row_entropy``: mass-weighted mean entropy (bits) of the row-normalised matrix (0 = every affix
    of X goes to one affix of Y); ``diag_mass``: share of the mass on the diagonal ("unchanged");
    ``matched_mass``: share on the best one-to-one assignment; ``n_moved``: affixes the assignment
    maps elsewhere; ``jsd_identity`` / ``jsd_matched``: JSD of the affix distributions of X and Y
    (``count_x``, ``count_y``, indexed like M) before and after relabelling X by the assignment.
    The assignment itself is returned under ``matching`` (dict a -> b).
    """
    w = M.to_numpy(dtype=float, copy=True)
    tot = w.sum()
    labels = list(M.index)
    if tot <= 0:
        return {"row_entropy": float("nan"), "diag_mass": float("nan"), "matched_mass": float("nan"),
                "n_moved": 0, "jsd_identity": float("nan"), "jsd_matched": float("nan"), "matching": {}}
    rs = w.sum(axis=1)
    ent = 0.0
    for r, s in zip(w, rs):
        if s > 0:
            pr = r[r > 0] / s
            ent += s / tot * float(-(pr * np.log2(pr)).sum())
    col = max_weight_matching(w)
    matched = float(sum(w[i, j] for i, j in enumerate(col) if j >= 0) / tot)
    cx = count_x.reindex(labels, fill_value=0).to_numpy(dtype=float, copy=True)
    cy = count_y.reindex(labels, fill_value=0).to_numpy(dtype=float, copy=True)
    relab = np.zeros_like(cx)
    for i, j in enumerate(col):
        relab[j if j >= 0 else i] += cx[i]
    matching = {labels[i]: labels[j] for i, j in enumerate(col) if j >= 0}
    return {"row_entropy": ent, "diag_mass": float(np.trace(w) / tot), "matched_mass": matched,
            "n_moved": int(sum(1 for i, j in enumerate(col) if j >= 0 and j != i and rs[i] > 0)),
            "jsd_identity": jsd(cx, cy), "jsd_matched": jsd(relab, cy), "matching": matching}


# --------------------------------------------------------------------------- #
# gradient along the page order
# --------------------------------------------------------------------------- #


def page_scores(
    page: Sequence, value: Sequence[str], label: Sequence[str], a: str = "A", b: str = "B",
    alpha: float = 0.5, min_tokens: int = 1,
) -> pd.DataFrame:
    """Leave-one-page-out projection of every page on the A-B axis.

    S(page) = sum_v p_page(v) log2(p_B(v) / p_A(v)), with p_A, p_B estimated from the tokens labelled
    ``a`` / ``b`` on all *other* pages (add-``alpha`` smoothing over the values seen anywhere).
    Positive = B-like.  Returns one row per page (first-appearance order) with ``page``, ``label``
    (majority), ``n_tokens``, ``score``; pages with fewer than ``min_tokens`` tokens are dropped.
    """
    d = pd.DataFrame({"page": list(page), "v": list(value), "lab": list(label)})
    vals = sorted(d["v"].unique())
    tab = pd.crosstab(d["page"], d["v"]).reindex(columns=vals, fill_value=0)
    lab = d.groupby("page", sort=False)["lab"].agg(lambda s: s.value_counts().index[0])
    pages = list(pd.unique(d["page"]))
    tab = tab.reindex(pages)
    lab = lab.reindex(pages)
    ca = tab[lab == a].sum(axis=0).to_numpy(dtype=float, copy=True)
    cb = tab[lab == b].sum(axis=0).to_numpy(dtype=float, copy=True)
    rows = []
    for p in pages:
        c = tab.loc[p].to_numpy(dtype=float, copy=True)
        n = c.sum()
        if n < min_tokens:
            continue
        ra = ca - (c if lab[p] == a else 0)
        rb = cb - (c if lab[p] == b else 0)
        pa = (ra + alpha) / (ra.sum() + alpha * len(vals))
        pb = (rb + alpha) / (rb.sum() + alpha * len(vals))
        rows.append({"page": p, "label": lab[p], "n_tokens": int(n), "score": float((c / n * np.log2(pb / pa)).sum())})
    return pd.DataFrame(rows)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rank correlation (average ranks for ties).

    >>> spearman([1, 2, 3, 4], [10, 20, 30, 40])
    1.0
    """
    rx = pd.Series(list(x), dtype=float).rank().to_numpy()
    ry = pd.Series(list(y), dtype=float).rank().to_numpy()
    if len(rx) < 3 or rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def order_drift(score: Sequence[float], order: Sequence[float], n_perm: int = 200, seed: int = 0) -> dict:
    """Spearman correlation of page scores with page order and its permutation p-value (two-sided).

    The null permutes the scores among the pages (``n_perm`` times).
    """
    s = np.asarray(score, dtype=float)
    o = np.asarray(order, dtype=float)
    rho = spearman(s, o)
    rng = np.random.default_rng(seed)
    sims = np.asarray([spearman(s[rng.permutation(len(s))], o) for _ in range(n_perm)])
    p = float((1 + (np.abs(sims) >= abs(rho) - 1e-12).sum()) / (1 + n_perm)) if np.isfinite(rho) else float("nan")
    return {"rho": rho, "p": p, "n_pages": int(len(s)), "perm_sd": float(np.nanstd(sims))}


def boundary_jumps(score: Sequence[float], label: Sequence[str]) -> dict:
    """Mean |score difference| of consecutive pages across a label change vs within the same label.

    Returns ``jump_boundary``, ``jump_within``, their ratio and the number of boundaries.
    """
    s = np.asarray(score, dtype=float)
    lab = np.asarray(label, dtype=object)
    if len(s) < 2:
        return {"jump_boundary": float("nan"), "jump_within": float("nan"), "ratio": float("nan"), "n_boundaries": 0}
    d = np.abs(np.diff(s))
    change = lab[1:] != lab[:-1]
    jb = float(d[change].mean()) if change.any() else float("nan")
    jw = float(d[~change].mean()) if (~change).any() else float("nan")
    return {"jump_boundary": jb, "jump_within": jw, "ratio": jb / jw if jw else float("nan"),
            "n_boundaries": int(change.sum())}


# --------------------------------------------------------------------------- #
# experiment 016: counterparts of language-specific stems (edit distance on glyph units)
# --------------------------------------------------------------------------- #
#
# Stems are handled as tuples of units (glyphs of ``collapse_eva``, v101 characters, BPE units,
# Latin letters) so that one edit = one unit.  ``nearest_counterparts`` uses a symmetric-delete
# index (every candidate and query indexed by all strings obtained by deleting up to ``max_dist``
# units; two sequences within Levenshtein distance d always share such a reduced string), so the
# search is exact for distances <= ``max_dist`` and far cheaper than all pairs.

Units = tuple


def levenshtein(a: Sequence, b: Sequence, max_dist: int | None = None) -> int:
    """Levenshtein distance of two unit sequences (insertion, deletion, substitution cost 1).

    With ``max_dist`` the function returns ``max_dist + 1`` as soon as the distance is known to
    exceed it.

    >>> levenshtein("chol", "shol"), levenshtein(("q", "o", "k"), ("o", "k"))
    (1, 1)
    """
    a, b = tuple(a), tuple(b)
    if max_dist is not None and abs(len(a) - len(b)) > max_dist:
        return max_dist + 1
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        if max_dist is not None and min(cur) > max_dist:
            return max_dist + 1
        prev = cur
    return prev[-1]


def edit_script(a: Sequence, b: Sequence) -> list[str]:
    """Edit operations of one optimal alignment turning ``a`` into ``b``.

    Every operation is a string ``kind:units@position`` with ``kind`` = ``ins`` / ``del`` / ``sub``
    (``sub`` units written ``x>y``) and ``position`` relative to ``a``: ``begin`` (first unit or in
    front of it), ``end`` (last unit or after it), ``mid`` otherwise.  Adjacent insertions and
    adjacent deletions are merged into one operation (``ins:ed@end``), so a script can be shorter
    than the Levenshtein distance.  Ties in the alignment prefer substitution, then deletion.

    >>> edit_script(("o", "k", "y"), ("q", "o", "k", "e", "d", "y"))
    ['ins:q@begin', 'ins:ed@mid']
    """
    a, b = tuple(a), tuple(b)
    n, m = len(a), len(b)
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        D[i][0] = i
    for j in range(1, m + 1):
        D[0][j] = j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            D[i][j] = min(D[i - 1][j] + 1, D[i][j - 1] + 1, D[i - 1][j - 1] + (a[i - 1] != b[j - 1]))
    raw: list[tuple[str, int, str]] = []  # (kind, index in a, units)
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and D[i][j] == D[i - 1][j - 1] + (a[i - 1] != b[j - 1]):
            if a[i - 1] != b[j - 1]:
                raw.append(("sub", i - 1, f"{a[i - 1]}>{b[j - 1]}"))
            i, j = i - 1, j - 1
        elif i > 0 and D[i][j] == D[i - 1][j] + 1:
            raw.append(("del", i - 1, a[i - 1]))
            i -= 1
        else:
            raw.append(("ins", i, b[j - 1]))
            j -= 1
    raw.reverse()
    merged: list[list] = []
    for kind, idx, u in raw:
        if merged and kind in ("ins", "del") and merged[-1][0] == kind and (
                (kind == "ins" and merged[-1][1] == idx) or (kind == "del" and merged[-1][2] == idx - 1)):
            merged[-1][2] = idx
            merged[-1][3] += u
        else:
            merged.append([kind, idx, idx, u])
    out = []
    for kind, start, end, u in merged:
        if kind == "ins":
            pos = "begin" if start == 0 else "end" if start == n else "mid"
        else:
            pos = "begin" if start == 0 else "end" if end == n - 1 else "mid"
        out.append(f"{kind}:{u}@{pos}")
    return out


def _deletes(units: Units, max_dist: int) -> set:
    out = {units}
    frontier = {units}
    for _ in range(max_dist):
        nxt = set()
        for u in frontier:
            for i in range(len(u)):
                nxt.add(u[:i] + u[i + 1:])
        nxt -= out
        out |= nxt
        frontier = nxt
    return out


def nearest_counterparts(
    queries: Iterable[Units], candidates: Iterable[Units], max_dist: int = 2,
    weights: Mapping[Units, float] | None = None,
) -> dict:
    """Nearest candidate (Levenshtein on units, exact up to ``max_dist``) of every query.

    Returns ``{query: (candidate, distance)}`` for the queries that have a candidate within
    ``max_dist`` and ``{query: None}`` for the others.  Ties go to the candidate with the larger
    ``weights`` value (e.g. token count), then to the lexicographically larger tuple.

    >>> nearest_counterparts([("s", "h", "o", "l")], [("c", "h", "o", "l"), ("d", "y")])
    {('s', 'h', 'o', 'l'): (('c', 'h', 'o', 'l'), 1)}
    """
    weights = dict(weights or {})
    index: dict = {}
    for c in dict.fromkeys(tuple(c) for c in candidates):
        for d in _deletes(c, max_dist):
            index.setdefault(d, []).append(c)
    out = {}
    for q in queries:
        q = tuple(q)
        seen = set()
        best, bd = None, max_dist + 1
        for d in _deletes(q, max_dist):
            for c in index.get(d, ()):
                if c in seen:
                    continue
                seen.add(c)
                dist = 0 if c == q else levenshtein(q, c, max_dist)
                if dist < bd or (dist == bd and best is not None and (weights.get(c, 0), c) > (weights.get(best, 0), best)):
                    best, bd = c, dist
        out[q] = (best, bd) if best is not None else None
    return out


def shuffle_units(stems: Iterable[Units], rng: np.random.Generator) -> list:
    """Length-preserving pooled shuffle: the units of all stems permuted and re-cut to the original lengths.

    A "random form" null: unit frequencies and stem lengths are kept, the composition of the
    stems is destroyed.
    """
    stems = [tuple(s) for s in stems]
    pool = [u for s in stems for u in s]
    pool = [pool[i] for i in rng.permutation(len(pool))]
    out, k = [], 0
    for s in stems:
        out.append(tuple(pool[k: k + len(s)]))
        k += len(s)
    return out


def counterpart_rates(
    count_x: Mapping[Units, float], count_y: Mapping[Units, float], min_count: int = 2, max_dist: int = 2,
    count_y_only: Mapping[Units, float] | None = None,
) -> dict:
    """Stems specific to X and their nearest counterparts in Y.

    X-only stems = stems with >= ``min_count`` tokens in X and no token in Y (``count_y_only`` if given,
    e.g. the real inventory when ``count_y`` is a shuffled one); candidates = stems with
    >= ``min_count`` tokens in ``count_y``.  Returns ``share_only`` (X tokens in X-only stems / all X tokens),
    ``n_only``, ``n_cand``, ``rate_le1`` / ``rate_le2`` (share of the X-only *tokens* whose stem has a
    counterpart within 1 / 2 edits), ``type_le1`` / ``type_le2`` (the same over stems) and ``pairs``
    (``{x: (y, distance)}`` for the X-only stems with a counterpart).
    """
    ref = count_y if count_y_only is None else count_y_only
    only = {s: c for s, c in count_x.items() if c >= min_count and not ref.get(s)}
    cand = {s: c for s, c in count_y.items() if c >= min_count}
    tot_x = float(sum(count_x.values()))
    tot_only = float(sum(only.values()))
    near = nearest_counterparts(only, cand, max_dist, weights=cand) if only and cand else {s: None for s in only}
    pairs = {s: r for s, r in near.items() if r is not None}
    le1 = sum(only[s] for s, (_, d) in pairs.items() if d <= 1)
    le2 = sum(only[s] for s, (_, d) in pairs.items() if d <= 2)
    return {"share_only": tot_only / tot_x if tot_x else float("nan"), "n_only": len(only), "n_cand": len(cand),
            "rate_le1": le1 / tot_only if tot_only else float("nan"),
            "rate_le2": le2 / tot_only if tot_only else float("nan"),
            "type_le1": sum(1 for _, d in pairs.values() if d <= 1) / len(only) if only else float("nan"),
            "type_le2": sum(1 for _, d in pairs.values() if d <= 2) / len(only) if only else float("nan"),
            "pairs": pairs}


def edit_scripts(pairs: Mapping[Units, tuple], direction: str = "xy") -> dict:
    """Edit script (:func:`edit_script`) of every counterpart pair ``{x: (y, d)}``.

    ``direction='xy'``: operations turn x into y; ``'yx'``: y into x.  Returns ``{x: [ops]}``.
    """
    out = {}
    for x, (y, _) in pairs.items():
        out[x] = edit_script(x, y) if direction == "xy" else edit_script(y, x)
    return out


def edit_inventory(scripts: Mapping[Units, Sequence[str]], weights: Mapping[Units, float]) -> pd.DataFrame:
    """Inventory of edit operations over counterpart pairs, ranked by token weight.

    One row per operation: ``op``, ``kind``, ``units``, ``position``, ``tokens`` (weight of the pairs
    containing it), ``n_pairs``, ``share`` (tokens / total weight of the pairs), ``rank`` and
    ``explained`` = share of the pairs' weight whose script uses only the operations ranked so far
    (how many operations cover 50 / 80 / 95 % of the tokens: :func:`coverage_levels`).
    """
    tok: dict = {}
    npair: dict = {}
    total = 0.0
    for x, ops in scripts.items():
        w = float(weights.get(x, 0))
        total += w
        for op in set(ops):
            tok[op] = tok.get(op, 0.0) + w
            npair[op] = npair.get(op, 0) + 1
    rows = sorted(tok.items(), key=lambda kv: (-kv[1], kv[0]))
    allowed: set = set()
    out = []
    for r, (op, w) in enumerate(rows, 1):
        allowed.add(op)
        expl = sum(float(weights.get(x, 0)) for x, ops in scripts.items() if set(ops) <= allowed)
        kind, rest = op.split(":", 1)
        units, pos = rest.rsplit("@", 1)
        out.append({"op": op, "kind": kind, "units": units, "position": pos, "tokens": w, "n_pairs": npair[op],
                    "share": w / total if total else float("nan"), "rank": r,
                    "explained": expl / total if total else float("nan")})
    return pd.DataFrame(out, columns=["op", "kind", "units", "position", "tokens", "n_pairs", "share", "rank",
                                      "explained"])


def coverage_levels(inventory: pd.DataFrame, levels: Sequence[float] = (0.5, 0.8, 0.95)) -> dict:
    """Number of top-ranked operations needed to explain each share of the tokens (``nan`` if never)."""
    out = {}
    for lv in levels:
        hit = inventory.index[inventory["explained"] >= lv - 1e-12] if len(inventory) else []
        out[lv] = int(inventory.loc[hit[0], "rank"]) if len(hit) else float("nan")
    return out


def substitution_sharing(
    count_x: Mapping[Units, float], count_y: Mapping[Units, float], pairs: Mapping[Units, tuple],
    scripts: Mapping[Units, Sequence[str]], ops_ranked: Sequence[str], ks: Sequence, min_count: int = 2,
) -> dict:
    """Stem sharing after rewriting the X-only stems whose edit script uses only the top-``k`` operations.

    For every ``k`` in ``ks`` (an int or ``'all'``) the X-only stems of ``pairs`` whose script is a
    subset of ``ops_ranked[:k]`` are replaced by their Y counterpart (counts added); the rest of X is
    unchanged.  Returns ``{k: {'jaccard', 'cover_xy', 'cover_yx', 'n_rewritten'}}`` of the rewritten
    X against Y (:func:`stem_sharing` with ``min_count``); ``k = 0`` is the sharing without rewriting.
    """
    out = {}
    for k in ks:
        allowed = set(ops_ranked if k == "all" else ops_ranked[:int(k)])
        new: dict = dict(count_x)
        n_rw = 0
        for x, (y, _) in pairs.items():
            if set(scripts[x]) <= allowed and x in new:
                new[y] = new.get(y, 0) + new.pop(x)
                n_rw += 1
        r = _sharing_counts(new, count_y, min_count)
        r["n_rewritten"] = n_rw
        out[k] = r
    return out


def _sharing_counts(cx: Mapping, cy: Mapping, min_count: int) -> dict:
    ix = {s for s, c in cx.items() if c >= min_count}
    iy = {s for s, c in cy.items() if c >= min_count}
    union = ix | iy
    tx, ty = float(sum(cx.values())), float(sum(cy.values()))
    return {"jaccard": len(ix & iy) / len(union) if union else float("nan"),
            "cover_xy": sum(c for s, c in cx.items() if s in iy) / tx if tx else float("nan"),
            "cover_yx": sum(c for s, c in cy.items() if s in ix) / ty if ty else float("nan")}


def profile_test(
    pairs: Mapping, prof_x: Mapping[Units, np.ndarray], prof_y: Mapping[Units, np.ndarray],
    n_perm: int = 200, seed: int = 0,
) -> dict:
    """Do counterpart pairs share a distribution (e.g. over sections) more than random re-pairings?

    ``pairs`` = ``{x: y}`` or ``{x: (y, d)}``; ``prof_x[x]`` / ``prof_y[y]`` = count vectors over the same
    categories.  Returns the mean JSD (bits) of the pairs, the mean and sd of ``n_perm`` permutations of
    the y side among the pairs, the one-sided p-value (real <= permuted), z-score and ``n_pairs``.
    """
    xs = list(pairs)
    ys = [pairs[x][0] if isinstance(pairs[x], tuple) and len(pairs[x]) == 2 and isinstance(pairs[x][1], (int, np.integer))
          else pairs[x] for x in xs]
    if not xs:
        return {"jsd_mean": float("nan"), "perm_mean": float("nan"), "perm_sd": float("nan"), "p": float("nan"),
                "z": float("nan"), "n_pairs": 0}
    px = [np.asarray(prof_x[x], dtype=float) for x in xs]
    py = [np.asarray(prof_y[y], dtype=float) for y in ys]
    real = float(np.mean([jsd(a, b) for a, b in zip(px, py)]))
    rng = np.random.default_rng(seed)
    sims = []
    for _ in range(n_perm):
        perm = rng.permutation(len(xs))
        sims.append(float(np.mean([jsd(px[i], py[perm[i]]) for i in range(len(xs))])))
    sims = np.asarray(sims)
    sd = float(sims.std())
    return {"jsd_mean": real, "perm_mean": float(sims.mean()), "perm_sd": sd,
            "p": float((1 + (sims <= real + 1e-12).sum()) / (1 + n_perm)),
            "z": (real - sims.mean()) / sd if sd > 0 else float("nan"), "n_pairs": len(xs)}
