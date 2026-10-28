"""Register comparison: sharing, affix mapping, matching, page scores (voytools.registers)."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import pytest

from voytools import registers as G


def test_draw_tokens_exact_size_and_exclusion():
    page = np.repeat(np.arange(10), 7)
    rng = np.random.default_rng(0)
    m = G.draw_tokens(page, range(10), 30, rng)
    assert m.sum() == 30
    assert len(np.unique(page[m])) == 5  # 4 whole pages + a cut one
    m2 = G.draw_tokens(page, range(10), 30, rng, exclude=np.unique(page[m]))
    assert not (m & m2).any()
    assert G.draw_tokens(page, [0, 1], 100, rng).sum() == 14
    first, second = G.split_pages(list("abcd"), [10, 10, 10, 10])
    assert first == ["a", "b"] and second == ["c", "d"]


def test_stem_sharing():
    r = G.stem_sharing(list("aaabbbccc"), list("aaabbbddd") + ["e"])
    assert r["jaccard"] == pytest.approx(0.5)
    assert r["cover_xy"] == pytest.approx(6 / 9)
    assert r["cover_yx"] == pytest.approx(6 / 10)
    assert r["n_x"] == 3 and r["n_y"] == 3


def test_hungarian_matches_brute_force():
    rng = np.random.default_rng(1)
    for n, m in ((4, 4), (5, 5), (3, 5), (5, 3)):
        w = rng.random((n, m))
        col = G.max_weight_matching(w)
        got = sum(w[i, j] for i, j in enumerate(col) if j >= 0)
        if n <= m:
            best = max(sum(w[i, p[i]] for i in range(n)) for p in itertools.permutations(range(m), n))
        else:
            best = max(sum(w[p[j], j] for j in range(m)) for p in itertools.permutations(range(n), m))
        assert got == pytest.approx(best)
        used = [j for j in col if j >= 0]
        assert len(used) == len(set(used)) == min(n, m)


def test_transfer_and_mapping_fit():
    rng = np.random.default_rng(2)
    stems = rng.choice(list("abcdefgh"), size=4000)
    sx = rng.choice(["x", "y", "z"], size=4000, p=[0.5, 0.3, 0.2])
    # Y: the same choice, written with a relabelled alphabet x->y, y->z, z->x
    sy = pd.Series(sx).map({"x": "y", "y": "z", "z": "x"}).to_numpy()
    half = slice(0, 2000), slice(2000, 4000)
    M, n = G.transfer_matrix(stems[half[0]], sx[half[0]], stems[half[1]], sy[half[1]], min_count=5)
    assert n == 8
    assert M.to_numpy().sum() == pytest.approx(8)
    cx = pd.Series(sx[half[0]]).value_counts()
    cy = pd.Series(sy[half[1]]).value_counts()
    fit = G.mapping_fit(M, cx, cy)
    assert fit["matching"] == {"x": "y", "y": "z", "z": "x"}
    assert fit["jsd_identity"] > 0.05 and fit["jsd_matched"] < 0.01
    assert fit["matched_mass"] > fit["diag_mass"]
    assert fit["n_moved"] == 3


def test_jsd():
    assert G.jsd([1, 0], [0, 1]) == pytest.approx(1.0)
    assert G.jsd([2, 2], [5, 5]) == pytest.approx(0.0)
    c = G.jsd_contributions([3, 1, 0], [1, 1, 2])
    assert c.sum() == pytest.approx(G.jsd([3, 1, 0], [1, 1, 2]))


def test_page_scores_and_drift():
    rng = np.random.default_rng(3)
    rows = []
    for p in range(40):
        lab = "A" if p < 20 else "B"
        pb = 0.2 + 0.6 * (p >= 20)
        vals = np.where(rng.random(50) < pb, "u", "v")
        rows += [(p, v, lab) for v in vals]
    d = pd.DataFrame(rows, columns=["page", "v", "lab"])
    sc = G.page_scores(d["page"], d["v"], d["lab"])
    assert len(sc) == 40
    assert (sc.loc[sc["label"] == "B", "score"] > 0).mean() > 0.9
    assert (sc.loc[sc["label"] == "A", "score"] < 0).mean() > 0.9
    jumps = G.boundary_jumps(sc["score"], sc["label"])
    assert jumps["n_boundaries"] == 1 and jumps["ratio"] > 2
    a = sc[sc["label"] == "A"]
    dr = G.order_drift(a["score"], a["page"], n_perm=50, seed=0)
    assert dr["p"] > 0.01 and abs(dr["rho"]) < 0.6
    trend = G.order_drift(np.arange(20) + rng.normal(0, 1, 20), np.arange(20), n_perm=50, seed=0)
    assert trend["rho"] > 0.8 and trend["p"] < 0.05
    assert G.spearman([1, 2, 3], [3, 2, 1]) == pytest.approx(-1.0)
