"""Tests of voytools.locality (experiment 018)."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import pytest

from voytools import locality as L


def _toy(seed: int = 0, n_pages: int = 12, n_stems: int = 60):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(0.3, size=(n_pages, n_stems))
    counts[:, :5] += rng.poisson(3, size=(n_pages, 5))  # frequent stems
    folios = [f"f{i}" for i in range(n_pages)]
    weight = counts.sum(axis=1).astype(float)
    weight = weight / weight.sum()
    dist = L.distance_matrix(np.arange(n_pages))
    return counts, folios, weight, dist


def test_stem_page_counts_rows_follow_folios():
    dfd = pd.DataFrame({"folio": ["b", "a", "a", "c", "b"], "stem": ["x", "x", "y", "", "y"]})
    counts, stems = L.stem_page_counts(dfd, ["a", "b", "zz"])
    assert list(stems) == ["x", "y"]
    assert counts.tolist() == [[1, 1], [1, 1], [0, 0]]  # page c ignored, empty stem dropped


def test_rare_pair_matrices_match_exact_expectation():
    counts, folios, weight, dist = _toy()
    n, w = L.rare_pair_matrices(counts, weight, 2, 5)
    tot = counts.sum(axis=0)
    sel = (tot >= 2) & (tot <= 5)
    pairs = float((tot[sel] * (tot[sel] - 1) / 2).sum())
    # off-diagonal expectation sums to P (1 - sum w^2); observed = all pairs minus same-page pairs
    iu = np.triu_indices(len(folios), 1)
    assert np.isclose(w[iu].sum(), pairs * (1 - (weight ** 2).sum()))
    same_page = float((counts[:, sel] * (counts[:, sel] - 1) / 2).sum())
    assert np.isclose(n[iu].sum(), pairs - same_page)
    assert np.all(np.diag(n) == 0) and np.all(np.diag(w) == 0)
    assert np.allclose(n, n.T) and np.allclose(w, w.T)


def test_cosine_and_jaccard_ranges():
    counts, folios, weight, dist = _toy()
    c = L.cosine_matrix(counts, drop_top=5)
    j = L.jaccard_matrix(counts, drop_top=5)
    for m in (c, j):
        assert m.shape == (12, 12) and np.allclose(m, m.T)
        assert m.min() >= 0 and m.max() <= 1 + 1e-12 and np.all(np.diag(m) == 0)
    # cosine with the top stems kept is not smaller on average (shared function stems)
    assert L.cosine_matrix(counts, drop_top=0).mean() >= c.mean()


def test_isotonic_decreasing_pava():
    y = np.array([3.0, 1.0, 2.0, 2.0, 0.5, 1.0])
    fit = L.isotonic_decreasing(y, np.ones(6))
    assert np.all(np.diff(fit) <= 1e-12)
    assert np.allclose(fit, [3, 5 / 3, 5 / 3, 5 / 3, 0.75, 0.75])
    fit_w = L.isotonic_decreasing(y, np.array([1, 1, 1, 1, 10, 1]))
    assert fit_w[4] == pytest.approx((0.5 * 10 + 1) / 11)


def test_profile_and_smooth_fit():
    counts, folios, weight, dist = _toy()
    n, w = L.similarity_matrices("cos", counts, weight)
    prof = L.pooled_profile(n, w, dist, None, dmax=5)
    assert list(prof["d"]) == [1, 2, 3, 4, 5, 6]
    assert prof["den"].sum() == pytest.approx(12 * 11 / 2)
    f = L.smooth_fit(n, w, dist, None, dmax=5)
    assert len(f) == 12 and np.all(np.diff(f[1:]) <= 1e-12)
    r, wm, f2 = L.residuals(n, w, dist, None, dmax=5)
    assert np.allclose(f, f2) and np.all(np.diag(wm) == 0)
    # the binned non-isotonic model leaves zero residual mass in every bin
    f0 = L.smooth_fit(n, w, dist, None, dmax=5, isotonic=False)
    r0, wm0, _ = L.residuals(n, w, dist, None, dmax=5, f=f0)
    assert abs(r0[dist == 1].sum()) < 1e-9


def test_block_sums_and_segment_against_brute_force():
    rng = np.random.default_rng(1)
    n = 10
    r = rng.normal(size=(n, n))
    r = (r + r.T) / 2
    np.fill_diagonal(r, 0.0)
    s = L._block_sums(r)
    for a, b in ((0, 4), (3, 9), (2, 10), (5, 6)):
        assert s[a, b] == pytest.approx(sum(r[i, j] for i in range(a, b) for j in range(i + 1, b)))
    w = np.ones((n, n))
    np.fill_diagonal(w, 0.0)
    sc = L.block_scores(r, w, "gauss")
    segs = L.segment(r, w, 3, "gauss")
    best = max(sc[0, b1] + sc[b1, b2] + sc[b2, n] for b1, b2 in itertools.combinations(range(1, n), 2))
    b1, b2 = segs[3]
    assert sc[0, b1] + sc[b1, b2] + sc[b2, n] == pytest.approx(best)
    assert segs[1] == [] and len(segs[2]) == 1
    assert L.block_labels([3, 7], n).tolist() == [0, 0, 0, 1, 1, 1, 1, 2, 2, 2]


def test_segment_finds_planted_blocks():
    n = 30
    lab = np.repeat([0, 1, 2], 10)
    r = np.where(lab[:, None] == lab[None, :], 1.0, -0.5)
    np.fill_diagonal(r, 0.0)
    w = np.ones((n, n))
    np.fill_diagonal(w, 0.0)
    assert L.segment(r, w, 3, "gauss")[3] == [10, 20]
    # Poisson version: counts 3x the expectation inside the blocks
    e = np.full((n, n), 0.5)
    np.fill_diagonal(e, 0.0)
    o = np.where(lab[:, None] == lab[None, :], 1.5, 0.5) * (e > 0)
    assert L.segment(o, e, 3, "poisson")[3] == [10, 20]
    assert L.block_gain(r, w, [10, 20]) == pytest.approx(1.0)
    assert L.block_gain(r, w, []) == pytest.approx(r.sum() / w.sum())


def test_crossval_gain_zero_under_page_permutation():
    counts, folios, weight, dist = _toy(seed=3, n_pages=40, n_stems=400)
    g = np.mean([L.crossval_gain("rare", counts, weight, dist, None, s, k_max=6) for s in range(20)], axis=0)
    assert g.shape == (6,) and g[0] == pytest.approx(0.0, abs=1e-9)
    assert abs(g[1:]).max() < 0.5  # independent halves: no systematic gain (noise only)


def test_boundaries_and_matching():
    lab = ["a", "a", "b", "b", "b", "c", "a"]
    assert L.boundaries_of(lab).tolist() == [2, 5, 6]
    prec, rec = L.match_boundaries([1, 5, 10], [2, 5, 6], tol=1)
    assert prec == pytest.approx(2 / 3) and rec == pytest.approx(3 / 3)
    assert L.match_boundaries([], [1, 2]) == (float("nan"), 0.0) or np.isnan(L.match_boundaries([], [1, 2])[0])
    p0, r0 = L.random_boundary_baseline(100, 5, [10, 50, 90], tol=1, n_rep=300, seed=1)
    assert 0.05 < p0 < 0.2 and 0.05 < r0 < 0.3
    assert all(np.isnan(L.random_boundary_baseline(100, 5, [], tol=1)))


def test_matched_excess_zero_without_group_effect():
    n = 40
    rng = np.random.default_rng(0)
    dist = L.distance_matrix(np.arange(n))
    r = rng.normal(size=(n, n))
    r = (r + r.T) / 2
    np.fill_diagonal(r, 0.0)
    w = np.ones((n, n))
    np.fill_diagonal(w, 0.0)
    same = (np.arange(n) // 4)[:, None] == (np.arange(n) // 4)[None, :]
    res = L.matched_excess(r, w, dist, same, d_max=10)
    assert res["n_same"] == 40 * 3 // 2 and abs(res["excess"]) < 0.6
    r2 = r + np.where(same, 2.0, 0.0)
    np.fill_diagonal(r2, 0.0)
    assert L.matched_excess(r2, w, dist, same, d_max=10)["excess"] > 1.5


def test_neighbours_shares_and_anomalies():
    n = 20
    dist = L.distance_matrix(np.arange(n))
    score = -dist.astype(float)  # nearest pages are the most similar
    allowed = np.ones((n, n), dtype=bool)
    nb = L.neighbours(score, dist, allowed, k=3)
    assert nb["n_pages"] == n and nb["share_le2"] >= 0.9 and nb["share_gt16"] == 0 and nb["n_anomalous"] == 0
    score2 = dist.astype(float)  # farthest pages are the most similar
    nb2 = L.neighbours(score2, dist, allowed, k=3)
    assert nb2["share_gt16"] > 0 and nb2["n_anomalous"] == n


def test_quire_matrix_and_order_metrics():
    n = 24
    q = np.repeat(list("ABCDEF"), 4)
    dist = L.distance_matrix(np.arange(n))
    r = -dist.astype(float)
    np.fill_diagonal(r, 0.0)
    w = np.ones((n, n))
    np.fill_diagonal(w, 0.0)
    m = L.quire_matrix(r, w, q)
    assert list(m.index) == list("ABCDEF") and m.loc["A", "B"] > m.loc["A", "F"]
    met = L.quire_order_metrics(m)
    assert met["adjacent_share"] == 1.0 and met["rho"] > 0.9
    assert L.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert L.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
