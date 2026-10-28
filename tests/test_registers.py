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


# --------------------------------------------------------------------------- #
# experiment 016: counterparts, edit scripts, inventory, substitution
# --------------------------------------------------------------------------- #


def test_levenshtein_matches_brute_force_and_cap():
    rng = np.random.default_rng(3)
    alphabet = list("abc")

    def slow(a, b):
        if not a:
            return len(b)
        if not b:
            return len(a)
        return min(slow(a[1:], b) + 1, slow(a, b[1:]) + 1, slow(a[1:], b[1:]) + (a[0] != b[0]))

    for _ in range(60):
        a = tuple(rng.choice(alphabet, rng.integers(0, 5)))
        b = tuple(rng.choice(alphabet, rng.integers(0, 5)))
        d = slow(a, b)
        assert G.levenshtein(a, b) == d
        assert G.levenshtein(a, b, 1) == min(d, 2)
    assert G.levenshtein("kitten", "sitting") == 3


def test_edit_script_positions_and_merging():
    assert G.edit_script(("o", "k", "y"), ("q", "o", "k", "e", "d", "y")) == ["ins:q@begin", "ins:ed@mid"]
    assert G.edit_script(("ch", "o", "l"), ("sh", "o", "l")) == ["sub:ch>sh@begin"]
    assert G.edit_script("ol", "old") == ["ins:d@end"]
    assert G.edit_script("chedy", "chy") == ["del:ed@mid"]
    assert G.edit_script("daiin", "aiin") == ["del:d@begin"]
    assert G.edit_script("abc", "abc") == []
    # the script is a valid alignment: number of raw edits == distance
    rng = np.random.default_rng(4)
    for _ in range(40):
        a = tuple(rng.choice(list("xyz"), rng.integers(1, 6)))
        b = tuple(rng.choice(list("xyz"), rng.integers(1, 6)))
        ops = G.edit_script(a, b)
        n = sum(len(op.split(":", 1)[1].rsplit("@", 1)[0]) if not op.startswith("sub") else 1 for op in ops)
        assert n == G.levenshtein(a, b)


def test_nearest_counterparts_exact_within_two_edits():
    rng = np.random.default_rng(5)
    alphabet = list("abcd")
    cands = list({tuple(rng.choice(alphabet, rng.integers(1, 6))) for _ in range(80)})
    queries = [tuple(rng.choice(alphabet, rng.integers(1, 6))) for _ in range(40)]
    got = G.nearest_counterparts(queries, cands, max_dist=2)
    for q in queries:
        best = min(G.levenshtein(q, c) for c in cands)
        if best <= 2:
            assert got[q][1] == best
        else:
            assert got[q] is None
    # ties: the heavier candidate wins
    r = G.nearest_counterparts([("a",)], [("b",), ("c",)], weights={("b",): 1, ("c",): 5})
    assert r[("a",)] == (("c",), 1)


def test_counterpart_rates_inventory_and_substitution():
    cx = {("ch", "o", "l"): 5, ("sh", "o", "l"): 4, ("d", "y"): 3, ("o", "k", "y"): 2, ("x", "x", "x"): 1}
    cy = {("ch", "o", "l"): 6, ("ch", "e", "d", "y"): 3, ("q", "o", "k", "e", "d", "y"): 2, ("z", "z", "z", "z"): 2}
    r = G.counterpart_rates(cx, cy, min_count=2, max_dist=2)
    assert r["n_only"] == 3 and r["n_cand"] == 4
    assert r["share_only"] == pytest.approx(9 / 15)
    # shol -> chol (1), dy -> chedy (2), oky -> qokedy (3: no counterpart)
    assert r["rate_le1"] == pytest.approx(4 / 9)
    assert r["rate_le2"] == pytest.approx(7 / 9)
    assert r["pairs"][("sh", "o", "l")] == (("ch", "o", "l"), 1)
    assert ("o", "k", "y") not in r["pairs"]
    scripts = G.edit_scripts(r["pairs"], "xy")
    assert scripts[("sh", "o", "l")] == ["sub:sh>ch@begin"]
    inv = G.edit_inventory(scripts, cx)
    assert inv["op"].iloc[0] == "sub:sh>ch@begin" and inv["tokens"].iloc[0] == 4
    assert inv["explained"].iloc[0] == pytest.approx(4 / 7)
    assert inv["explained"].iloc[-1] == pytest.approx(1.0)
    assert G.coverage_levels(inv, (0.5, 0.95)) == {0.5: 1, 0.95: 2}
    sub = G.substitution_sharing(cx, cy, r["pairs"], scripts, inv["op"].tolist(), (0, 1, "all"), min_count=2)
    assert sub[0]["jaccard"] == pytest.approx(1 / 7) and sub[0]["n_rewritten"] == 0
    assert sub[1]["jaccard"] == pytest.approx(1 / 6) and sub[1]["n_rewritten"] == 1   # shol merged into chol
    assert sub["all"]["jaccard"] == pytest.approx(2 / 5) and sub["all"]["n_rewritten"] == 2
    assert sub["all"]["cover_xy"] > sub[0]["cover_xy"]


def test_shuffle_units_keeps_lengths_and_multiset():
    stems = [("a", "b"), ("c",), ("d", "e", "f")]
    out = G.shuffle_units(stems, np.random.default_rng(0))
    assert [len(s) for s in out] == [2, 1, 3]
    assert sorted(u for s in out for u in s) == list("abcdef")


def test_profile_test_detects_shared_profiles():
    rng = np.random.default_rng(6)
    prof_x, prof_y, pairs = {}, {}, {}
    for i in range(30):
        base = rng.dirichlet(np.ones(4) * 0.3)
        prof_x[f"x{i}"] = rng.multinomial(40, base)
        prof_y[f"y{i}"] = rng.multinomial(40, base)
        pairs[f"x{i}"] = (f"y{i}", 1)
    r = G.profile_test(pairs, prof_x, prof_y, n_perm=100, seed=0)
    assert r["n_pairs"] == 30 and r["jsd_mean"] < r["perm_mean"] and r["p"] < 0.05
    assert G.profile_test({}, {}, {})["n_pairs"] == 0
