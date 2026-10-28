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
"""

from __future__ import annotations

from typing import Sequence

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
