"""Locality of the page vocabulary: smooth distance decay versus blocks (experiment 018).

The unit is a *page* (folio with paragraph text).  From the stem counts per page
(:func:`stem_page_counts`) three page similarity measures are built, each as a pair of matrices
``(N, W)`` whose pooled ratio ``sum(N) / sum(W)`` over a set of page pairs is the similarity of that
set (:func:`similarity_matrices`):

* ``'rare'`` - sharing of rare stems (2-5 occurrences, experiment 015): ``N`` = observed number of
  occurrence pairs between the two pages, ``W`` = its expectation when every occurrence sits on a
  random size-weighted page (``2 w_i w_j`` times the number of pairs); the pooled ratio is the
  observed/expected ratio of 015;
* ``'cos'`` - cosine of the ``log(1+f)`` stem vectors without the top function stems (006),
  ``W = 1``;
* ``'jac'`` - Jaccard of the sets of non-hapax stems (without the top function stems), ``W = 1``.

For ``'rare'`` the page weights are the size-weighted expectation, so page size is accounted for;
``'cos'`` and ``'jac'`` keep the page-size dependence of 006 (small pages are less similar to
everything), which page permutation nulls share.

A *smooth distance model* ``f(d)`` is the pooled similarity per page distance (1..``dmax``, then one
constant), made monotonically decreasing by the pool-adjacent-violators algorithm
(:func:`smooth_fit`); the residual is ``R = N - f(D) W`` (:func:`residuals`).  *Blocks* are the
optimal segmentation of the page sequence into ``K`` contiguous blocks maximising the residual sum
inside blocks (dynamic programming over one-sided deviance block scores, :func:`segment`,
:func:`block_scores`); the *gain* of a segmentation is the pooled
residual inside the blocks (:func:`block_gain`), cross-validated by splitting the stems into two
halves (:func:`split_stems`, :func:`crossval_gain`).  Change points are compared with section /
language / hand / quire boundaries (:func:`match_boundaries`, :func:`random_boundary_baseline`).
Physical order: distance-matched excess of same-bifolio / same-quire pairs
(:func:`matched_excess`), nearest-neighbour anomalies (:func:`neighbours`) and the quire x quire
residual matrix with a spectral seriation (:func:`quire_matrix`, :func:`spectral_order`,
:func:`spearman`).  No scipy.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

DMAX = 40
DISTANCE_BINS = (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 41)   # lower edges; the last bin is > DMAX
MEASURES = ("rare", "cos", "jac")
NEAR_EDGES = (2, 5, 16)


# --------------------------------------------------------------------------- #
# counts and similarity matrices
# --------------------------------------------------------------------------- #


def stem_page_counts(dfd: pd.DataFrame, folios: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """Dense ``pages x stems`` count matrix of the non-empty stems of ``dfd`` on the pages ``folios``.

    Rows follow ``folios`` (pages absent from ``dfd`` are zero rows); tokens on other pages are
    ignored.  Returns ``(counts, stems)``.
    """
    keep = dfd["stem"].astype(str) != ""
    fol = dfd.loc[keep, "folio"].astype(str).to_numpy()
    stem = dfd.loc[keep, "stem"].astype(str).to_numpy()
    page_idx = {p: i for i, p in enumerate(folios)}
    ok = np.array([f in page_idx for f in fol], dtype=bool)
    rows = np.array([page_idx[f] for f in fol[ok]], dtype=np.int64)
    codes, stems = pd.factorize(stem[ok], sort=True)
    counts = np.zeros((len(folios), len(stems)), dtype=np.int64)
    np.add.at(counts, (rows, codes), 1)
    return counts, np.asarray(stems, dtype=object)


def rare_pair_matrices(counts: np.ndarray, weight: np.ndarray, min_count: int = 2, max_count: int = 5,
                       columns: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Observed and expected numbers of occurrence pairs of rare stems between every two pages.

    ``counts`` is ``pages x stems``; stems with total count in ``[min_count, max_count]`` (restricted to
    ``columns`` if given) are used.  ``N[i, j] = sum_s c_is c_js`` (i != j); ``W[i, j] = 2 w_i w_j P``
    with ``P`` the number of occurrence pairs and ``w`` the page weights (size-weighted page
    permutation of 015: an unordered pair lands on pages ``{i, j}`` with probability ``2 w_i w_j``).
    Diagonals are zero.
    """
    tot = counts.sum(axis=0)
    sel = (tot >= min_count) & (tot <= max_count)
    if columns is not None:
        m = np.zeros(counts.shape[1], dtype=bool)
        m[np.asarray(columns, dtype=np.int64)] = True
        sel &= m
    c = counts[:, sel].astype(float)
    n_obs = c @ c.T
    k = tot[sel].astype(float)
    pairs = float((k * (k - 1) / 2).sum())
    w = np.asarray(weight, dtype=float)
    w_exp = 2.0 * np.outer(w, w) * pairs
    np.fill_diagonal(n_obs, 0.0)
    np.fill_diagonal(w_exp, 0.0)
    return n_obs, w_exp


def _drop_top(counts: np.ndarray, drop_top: int) -> np.ndarray:
    tot = counts.sum(axis=0)
    keep = np.ones(counts.shape[1], dtype=bool)
    if drop_top > 0 and counts.shape[1] > drop_top:
        order = np.argsort(-tot, kind="stable")
        keep[order[:drop_top]] = False
    return keep


def cosine_matrix(counts: np.ndarray, drop_top: int = 30, columns: np.ndarray | None = None) -> np.ndarray:
    """Cosine of the ``log(1+count)`` page vectors without the ``drop_top`` most frequent stems (006); diagonal 0."""
    keep = _drop_top(counts, drop_top)
    if columns is not None:
        m = np.zeros(counts.shape[1], dtype=bool)
        m[np.asarray(columns, dtype=np.int64)] = True
        keep &= m
    x = np.log1p(counts[:, keep].astype(float))
    norm = np.sqrt((x * x).sum(axis=1))
    norm[norm == 0] = 1.0
    s = (x @ x.T) / np.outer(norm, norm)
    np.fill_diagonal(s, 0.0)
    return s


def jaccard_matrix(counts: np.ndarray, drop_top: int = 30, min_count: int = 2,
                   columns: np.ndarray | None = None) -> np.ndarray:
    """Jaccard index of the sets of stems with corpus count >= ``min_count`` (without the top stems); diagonal 0."""
    keep = _drop_top(counts, drop_top) & (counts.sum(axis=0) >= min_count)
    if columns is not None:
        m = np.zeros(counts.shape[1], dtype=bool)
        m[np.asarray(columns, dtype=np.int64)] = True
        keep &= m
    b = (counts[:, keep] > 0).astype(float)
    inter = b @ b.T
    size = b.sum(axis=1)
    union = size[:, None] + size[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        j = np.where(union > 0, inter / union, 0.0)
    np.fill_diagonal(j, 0.0)
    return j


def similarity_matrices(measure: str, counts: np.ndarray, weight: np.ndarray, columns: np.ndarray | None = None,
                        drop_top: int = 30, min_count: int = 2, max_count: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """``(N, W)`` of one measure (:data:`MEASURES`); ``W`` is 1 off the diagonal for ``'cos'`` / ``'jac'``."""
    if measure == "rare":
        return rare_pair_matrices(counts, weight, min_count, max_count, columns)
    if measure == "cos":
        s = cosine_matrix(counts, drop_top, columns)
    elif measure == "jac":
        s = jaccard_matrix(counts, drop_top, min_count, columns)
    else:
        raise ValueError(f"unknown measure {measure!r}")
    w = np.ones_like(s)
    np.fill_diagonal(w, 0.0)
    return s, w


def split_stems(n_stems: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Random split of the stem columns ``0..n_stems-1`` into two halves (index arrays)."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_stems)
    h = n_stems // 2
    return np.sort(perm[:h]), np.sort(perm[h:])


# --------------------------------------------------------------------------- #
# distance model and residuals
# --------------------------------------------------------------------------- #


def distance_matrix(rank: Sequence[int]) -> np.ndarray:
    """``|rank_i - rank_j|`` for the pages."""
    r = np.asarray(rank, dtype=np.int64)
    return np.abs(r[:, None] - r[None, :])


def isotonic_decreasing(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted least-squares fit of a non-increasing sequence (pool-adjacent-violators)."""
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    blocks: list[list[float]] = []  # [value, weight, size]
    for yi, wi in zip(y, w):
        blocks.append([yi, wi, 1])
        while len(blocks) > 1 and blocks[-2][0] < blocks[-1][0]:
            v2, w2, n2 = blocks.pop()
            v1, w1, n1 = blocks.pop()
            tot = w1 + w2
            v = (v1 * w1 + v2 * w2) / tot if tot > 0 else (v1 + v2) / 2
            blocks.append([v, tot, n1 + n2])
    return np.concatenate([[v] * n for v, _, n in blocks])


def pooled_profile(n_mat: np.ndarray, w_mat: np.ndarray, dist: np.ndarray, mask: np.ndarray | None = None,
                   dmax: int = DMAX) -> pd.DataFrame:
    """Pooled similarity per distance ``1..dmax`` and ``dmax+1`` (all larger distances), over the pairs in ``mask``.

    Columns ``d``, ``num``, ``den``, ``ratio`` (``num / den``; NaN where ``den`` is 0).  Each unordered
    pair is counted once.
    """
    iu = np.triu_indices(n_mat.shape[0], 1)
    n, w, d = n_mat[iu], w_mat[iu], dist[iu]
    if mask is not None:
        m = mask[iu]
        n, w, d = n[m], w[m], d[m]
    dd = np.minimum(d, dmax + 1)
    num = np.bincount(dd, weights=n, minlength=dmax + 2)[1:]
    den = np.bincount(dd, weights=w, minlength=dmax + 2)[1:]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(den > 0, num / den, np.nan)
    return pd.DataFrame({"d": np.arange(1, dmax + 2), "num": num, "den": den, "ratio": ratio})


def smooth_fit(n_mat: np.ndarray, w_mat: np.ndarray, dist: np.ndarray, mask: np.ndarray | None = None,
               dmax: int = DMAX, bins: Sequence[int] = DISTANCE_BINS, isotonic: bool = True) -> np.ndarray:
    """Non-increasing distance model ``f`` indexed by distance ``0..max(dist)`` (``f[0]`` unused, = ``f[1]``).

    The pooled profile (:func:`pooled_profile`) is first pooled into the distance bins ``bins``
    (lower edges; the last bin holds everything beyond ``dmax``), then made non-increasing by the
    isotonic fit weighted by the pooled ``den``; ``f`` is piecewise constant on the bins.  Binning
    keeps the fit stable when the pairs per single distance are few (a half of the rare stems gives
    only a handful of pairs at ``d = 1``), which otherwise biases the isotonic fit upwards at short
    distances.  With ``isotonic=False`` the binned pooled ratios are used as they are (the residual
    then sums to zero inside every bin: the unbiased evaluation model of :func:`crossval_gain`).
    """
    prof = pooled_profile(n_mat, w_mat, dist, mask, dmax)
    edges = np.asarray(bins, dtype=np.int64)
    d = prof["d"].to_numpy()
    b = np.searchsorted(edges, d, side="right") - 1
    num = np.bincount(b, weights=prof["num"].to_numpy(), minlength=len(edges))
    den = np.bincount(b, weights=prof["den"].to_numpy(), minlength=len(edges))
    with np.errstate(invalid="ignore", divide="ignore"):
        y = np.where(den > 0, num / den, 0.0)
    fit = isotonic_decreasing(y, den) if isotonic else np.asarray(y, dtype=float)
    top = int(dist.max()) if dist.size else 0
    f = np.empty(max(top, dmax + 1) + 1)
    dd = np.arange(1, len(f))
    f[1:] = fit[np.searchsorted(edges, np.minimum(dd, dmax + 1), side="right") - 1]
    f[0] = f[1]
    return f


def residuals(n_mat: np.ndarray, w_mat: np.ndarray, dist: np.ndarray, mask: np.ndarray | None = None,
              dmax: int = DMAX, f: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(R, Wm, f)``: residual ``N - f(D) W`` and the weights, both zero outside ``mask`` and on the diagonal.

    ``f`` (:func:`smooth_fit`) is fitted on ``(n_mat, w_mat)`` unless given.
    """
    if f is None:
        f = smooth_fit(n_mat, w_mat, dist, mask, dmax)
    w = w_mat.copy()
    if mask is not None:
        w = np.where(mask, w, 0.0)
    np.fill_diagonal(w, 0.0)
    r = np.where(w > 0, n_mat - f[dist] * w_mat, 0.0)
    return r, w, f


# --------------------------------------------------------------------------- #
# segmentation
# --------------------------------------------------------------------------- #


def _block_sums(r: np.ndarray) -> np.ndarray:
    """``S[a, b]`` = sum of ``r[i, j]`` over ``a <= i < j < b`` (symmetric ``r`` with zero diagonal)."""
    n = r.shape[0]
    p = np.zeros((n + 1, n + 1))
    p[1:, 1:] = r.cumsum(axis=0).cumsum(axis=1)
    diag = np.diag(p)
    return (diag[None, :] - p - p.T + diag[:, None]) / 2.0


def block_scores(num: np.ndarray, den: np.ndarray, kind: str) -> np.ndarray:
    """``S[a, b]``: one-sided deviance evidence that the block ``a..b-1`` is more similar inside than the model.

    ``kind='poisson'`` (rare stem pairs): ``num`` = observed pair counts ``O``, ``den`` = expected
    ``f(d) W``; score ``O log(O / E) - (O - E)`` when ``O > E``, else 0 (profile log-likelihood ratio
    of a block-specific rate).  ``kind='gauss'`` (cosine / Jaccard): ``num`` = residuals, ``den`` =
    pair weights; score ``(sum R)^2 / (2 sigma^2 sum W)`` when ``sum R > 0``, else 0, with ``sigma^2``
    the pooled residual variance.  Blocks without pairs score 0.
    """
    t_num, t_den = _block_sums(num), _block_sums(den)
    ok = t_den > 0
    if kind == "poisson":
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(ok & (t_num > 0), t_num / np.where(ok, t_den, 1.0), 1.0)
            sc = np.where(ok & (t_num > t_den), t_num * np.log(ratio) - (t_num - t_den), 0.0)
        return sc
    if kind == "gauss":
        iu = np.triu_indices(num.shape[0], 1)
        w = den[iu]
        sigma2 = float((num[iu] ** 2 * (w > 0)).sum() / max(w.sum(), 1e-12))
        sigma2 = max(sigma2, 1e-12)
        with np.errstate(invalid="ignore", divide="ignore"):
            sc = np.where(ok & (t_num > 0), t_num ** 2 / (2.0 * sigma2 * np.where(ok, t_den, 1.0)), 0.0)
        return sc
    raise ValueError("kind must be 'poisson' or 'gauss'")


def segment(num: np.ndarray, den: np.ndarray, k_max: int, kind: str = "gauss") -> dict[int, list[int]]:
    """Optimal segmentations into ``K = 1..k_max`` contiguous blocks maximising the summed block evidence.

    Dynamic programming over :func:`block_scores` (``O(K n^2)``).  Returns ``{K: boundaries}`` with
    ``boundaries`` the start indices of blocks 2..K (sorted, between 1 and n-1); ``K = 1`` has no
    boundary.
    """
    n = num.shape[0]
    k_max = max(1, min(k_max, n))
    s = block_scores(num, den, kind)
    neg = -np.inf
    best = np.full((k_max + 1, n + 1), neg)
    arg = np.zeros((k_max + 1, n + 1), dtype=np.int64)
    best[1, 1:] = s[0, 1:]
    for k in range(2, k_max + 1):
        for b in range(k, n + 1):
            cand = best[k - 1, k - 1:b] + s[k - 1:b, b]
            i = int(np.argmax(cand))
            best[k, b], arg[k, b] = cand[i], k - 1 + i
    out: dict[int, list[int]] = {}
    for k in range(1, k_max + 1):
        bounds, b = [], n
        for kk in range(k, 1, -1):
            b = int(arg[kk, b])
            bounds.append(b)
        out[k] = sorted(bounds)
    return out


def segmentation_inputs(r: np.ndarray, w_model: np.ndarray, measure: str) -> tuple[np.ndarray, np.ndarray, str]:
    """``(num, den, kind)`` for :func:`segment`: observed / expected pair counts (Poisson) for ``'rare'``,
    residuals / weights (Gaussian) otherwise.  ``w_model`` is ``f(d) W`` for ``'rare'`` and ``W`` otherwise."""
    if measure == "rare":
        return r + w_model, w_model, "poisson"
    return r, w_model, "gauss"


def block_labels(bounds: Sequence[int], n: int) -> np.ndarray:
    """Block index of every position for block start positions ``bounds``."""
    lab = np.zeros(n, dtype=np.int64)
    for b in bounds:
        lab[b:] += 1
    return lab


def block_gain(r: np.ndarray, w: np.ndarray, bounds: Sequence[int]) -> float:
    """Pooled residual ``sum R / sum W`` over the page pairs inside the blocks (NaN if no weight)."""
    lab = block_labels(bounds, r.shape[0])
    inside = lab[:, None] == lab[None, :]
    den = float(w[inside].sum())
    return float(r[inside].sum() / den) if den > 0 else float("nan")


def gain_curve(segs: dict[int, list[int]], r_eval: np.ndarray, w_eval: np.ndarray, k_max: int) -> np.ndarray:
    """Gains for ``K = 1..k_max`` of the segmentations ``segs``: pooled residual measured on ``(r_eval, w_eval)``."""
    return np.array([block_gain(r_eval, w_eval, segs[k]) if k in segs else np.nan for k in range(1, k_max + 1)])


def crossval_gain(measure: str, counts: np.ndarray, weight: np.ndarray, dist: np.ndarray, mask: np.ndarray | None,
                  seed: int, k_max: int = 30, dmax: int = DMAX, f: np.ndarray | None = None, **kw) -> np.ndarray:
    """Cross-validated gain curve: stems split in two halves (``seed``), blocks on half 1, gain on half 2.

    Blocks are searched on the residuals of half 1 against the distance model ``f`` of all stems
    (fitted here if not given).  The gain is the pooled residual of half 2 inside the blocks against
    the *binned pooled profile of half 2 itself* (:func:`smooth_fit` with ``isotonic=False``), which
    sums to zero within every distance bin: any set of blocks that is independent of half 2 has
    expected gain 0, whatever the distance structure.  (A model shared between the halves, or an
    isotonic fit of the noisy half, both bias the gain downwards at short distances, where almost
    every pair lies inside a block.)  ``kw`` go to :func:`similarity_matrices`.
    """
    if f is None:
        n_all, w_all = similarity_matrices(measure, counts, weight, None, **kw)
        f = smooth_fit(n_all, w_all, dist, mask, dmax)
    h1, h2 = split_stems(counts.shape[1], seed)
    n1, w1 = similarity_matrices(measure, counts, weight, h1, **kw)
    n2, w2 = similarity_matrices(measure, counts, weight, h2, **kw)
    r1, wm1, _ = residuals(n1, w1, dist, mask, dmax, f)
    f2 = smooth_fit(n2, w2, dist, mask, dmax, isotonic=False)
    r2, wm2, _ = residuals(n2, w2, dist, mask, dmax, f2)
    num, den, kind = segmentation_inputs(r1, f[dist] * wm1 if measure == "rare" else wm1, measure)
    segs = segment(num, den, k_max, kind)
    return gain_curve(segs, r2, wm2, k_max)


# --------------------------------------------------------------------------- #
# change points
# --------------------------------------------------------------------------- #


def boundaries_of(labels: Sequence) -> np.ndarray:
    """Positions ``p`` (1..n-1) where ``labels[p] != labels[p-1]``."""
    lab = np.asarray(labels, dtype=object)
    return np.flatnonzero(lab[1:] != lab[:-1]) + 1


def match_boundaries(found: Sequence[int], ref: Sequence[int], tol: int = 1) -> tuple[float, float]:
    """Precision (share of ``found`` within ``tol`` of a reference boundary) and recall (share of ``ref`` hit)."""
    f = np.asarray(sorted(found), dtype=np.int64)
    r = np.asarray(sorted(ref), dtype=np.int64)
    if len(f) == 0 or len(r) == 0:
        return (float("nan") if len(f) == 0 else 0.0), (float("nan") if len(r) == 0 else 0.0)
    close = np.abs(f[:, None] - r[None, :]) <= tol
    return float(close.any(axis=1).mean()), float(close.any(axis=0).mean())


def random_boundary_baseline(n: int, k: int, ref: Sequence[int], tol: int = 1, n_rep: int = 200,
                             seed: int = 0) -> tuple[float, float]:
    """Mean precision and recall of ``k`` random distinct boundaries among positions ``1..n-1``."""
    rng = np.random.default_rng(seed)
    k = min(k, n - 1)
    if k <= 0 or len(ref) == 0:
        return float("nan"), float("nan")
    p, r = [], []
    for _ in range(n_rep):
        found = rng.choice(np.arange(1, n), size=k, replace=False)
        a, b = match_boundaries(found, ref, tol)
        p.append(a)
        r.append(b)
    return float(np.nanmean(p)), float(np.nanmean(r))


# --------------------------------------------------------------------------- #
# physical order
# --------------------------------------------------------------------------- #


def matched_excess(r: np.ndarray, w: np.ndarray, dist: np.ndarray, same: np.ndarray,
                   d_max: int | None = None) -> dict:
    """Pooled residual of the ``same`` pairs minus that of the other pairs at the same distances.

    Only pairs with ``w > 0`` count; ``d_max`` restricts the distances.  The matched value is the
    average over distances of the other pairs' pooled residual, weighted by the ``same`` pairs' weight
    at that distance.  Returns ``n_same``, ``same``, ``matched``, ``excess``.
    """
    iu = np.triu_indices(r.shape[0], 1)
    rr, ww, dd, ss = r[iu], w[iu], dist[iu], same[iu]
    ok = ww > 0
    if d_max is not None:
        ok &= dd <= d_max
    rr, ww, dd, ss = rr[ok], ww[ok], dd[ok], ss[ok]
    n_same = int(ss.sum())
    if n_same == 0 or (~ss).sum() == 0:
        return {"n_same": n_same, "same": np.nan, "matched": np.nan, "excess": np.nan}
    val_same = float(rr[ss].sum() / ww[ss].sum())
    top = int(dd.max()) + 1
    num_o = np.bincount(dd[~ss], weights=rr[~ss], minlength=top)
    den_o = np.bincount(dd[~ss], weights=ww[~ss], minlength=top)
    wt_s = np.bincount(dd[ss], weights=ww[ss], minlength=top)
    use = (den_o > 0) & (wt_s > 0)
    if not use.any():
        return {"n_same": n_same, "same": val_same, "matched": np.nan, "excess": np.nan}
    matched = float((num_o[use] / den_o[use] * wt_s[use]).sum() / wt_s[use].sum())
    return {"n_same": n_same, "same": val_same, "matched": matched, "excess": val_same - matched}


def neighbours(score: np.ndarray, dist: np.ndarray, allowed: np.ndarray, k: int = 3,
               edges: Sequence[int] = NEAR_EDGES, far: int = 8) -> dict:
    """Top-``k`` most similar pages per page (by ``score``, among ``allowed`` pairs) and their distances.

    Returns ``share_le<e>`` for ``e`` in ``edges[:-1]`` (share of neighbour distances <= e),
    ``share_gt<edges[-1]>``, ``n_pages`` (pages with k allowed partners), ``anomalous`` (indices of
    the pages whose ``k`` neighbours are all farther than ``far`` pages) and ``neighbour_dist``
    (``pages x k`` array, -1 where missing).
    """
    n = score.shape[0]
    nd = -np.ones((n, k), dtype=np.int64)
    for i in range(n):
        cand = np.flatnonzero(allowed[i])
        cand = cand[cand != i]
        if len(cand) < k:
            continue
        order = cand[np.argsort(-score[i, cand], kind="stable")[:k]]
        nd[i] = dist[i, order]
    have = nd[:, 0] >= 0
    d = nd[have].ravel()
    out = {"n_pages": int(have.sum())}
    for e in edges[:-1]:
        out[f"share_le{e}"] = float((d <= e).mean()) if d.size else np.nan
    out[f"share_gt{edges[-1]}"] = float((d > edges[-1]).mean()) if d.size else np.nan
    out["anomalous"] = np.flatnonzero(have & (nd.min(axis=1) > far))
    out["n_anomalous"] = int(len(out["anomalous"]))
    out["neighbour_dist"] = nd
    return out


def quire_matrix(r: np.ndarray, w: np.ndarray, quires: Sequence[str]) -> pd.DataFrame:
    """Pooled residual ``sum R / sum W`` for every pair of quires (rows/columns in order of first appearance)."""
    q = np.asarray(quires, dtype=object)
    names = list(dict.fromkeys(q.tolist()))
    idx = {name: i for i, name in enumerate(names)}
    code = np.array([idx[x] for x in q])
    m = len(names)
    num = np.zeros((m, m))
    den = np.zeros((m, m))
    np.add.at(num, (code[:, None].repeat(len(q), 1), code[None, :].repeat(len(q), 0)), r)
    np.add.at(den, (code[:, None].repeat(len(q), 1), code[None, :].repeat(len(q), 0)), w)
    with np.errstate(invalid="ignore", divide="ignore"):
        val = np.where(den > 0, num / den, np.nan)
    return pd.DataFrame(val, index=names, columns=names)


def spectral_order(sim: np.ndarray) -> np.ndarray:
    """Seriation by the Fiedler vector of the Laplacian of ``sim`` (shifted to be non-negative off the diagonal)."""
    a = np.array(sim, dtype=float)
    a = np.where(np.isfinite(a), a, np.nanmin(a[np.isfinite(a)]) if np.isfinite(a).any() else 0.0)
    np.fill_diagonal(a, 0.0)
    off = a[~np.eye(len(a), dtype=bool)]
    a = a - (off.min() if off.size else 0.0)
    np.fill_diagonal(a, 0.0)
    lap = np.diag(a.sum(axis=1)) - a
    vals, vecs = np.linalg.eigh(lap)
    return np.argsort(vecs[:, 1], kind="stable") if len(vals) > 1 else np.arange(len(a))


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rank correlation (average ranks for ties), without scipy."""
    rx = pd.Series(np.asarray(x, dtype=float)).rank().to_numpy()
    ry = pd.Series(np.asarray(y, dtype=float)).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def quire_order_metrics(q: pd.DataFrame) -> dict:
    """``adjacent_share``: share of quires whose most similar quire is a reading-order neighbour; ``rho``:
    |Spearman| of the spectral order with the reading order."""
    m = len(q)
    if m < 3:
        return {"adjacent_share": np.nan, "rho": np.nan}
    v = q.to_numpy(dtype=float).copy()
    np.fill_diagonal(v, -np.inf)
    v = np.where(np.isfinite(v) | np.isneginf(v), v, -np.inf)
    nearest = np.argmax(v, axis=1)
    adj = float(np.mean(np.abs(nearest - np.arange(m)) == 1))
    order = spectral_order(q.to_numpy(dtype=float))
    pos = np.empty(m)
    pos[order] = np.arange(m)
    rho = abs(spearman(pos, np.arange(m)))
    return {"adjacent_share": adj, "rho": rho}
