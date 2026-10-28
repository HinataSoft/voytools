"""Line-boundary MI and vertical similarity (voytools.stats, experiment 010)."""

from __future__ import annotations

import numpy as np

from voytools import stats as S
from tests.conftest import make_tokens


def _paragraphs(lines, page_size):
    """Token table with one marked paragraph per page (``<%>`` ... ``<$>``)."""
    marked = []
    for i, toks in enumerate(lines):
        toks = list(toks)
        if i % page_size == 0:
            toks[0] = "<%>" + toks[0]
        if i % page_size == page_size - 1 or i == len(lines) - 1:
            toks[-1] = toks[-1] + "<$>"
        marked.append(toks)
    return make_tokens(marked, page_size=page_size)


def test_pair_kinds_and_boundary_mi():
    # two pages, one paragraph each: inside pairs, line pairs, a page break
    lines = [["ab", "ba", "ab"], ["ba", "ab"], ["ab", "ba"], ["ba", "ab", "ba"]]
    df = _paragraphs(lines, 2)
    kinds = S._pair_kinds(df).tolist()
    assert kinds == ["inside", "inside", "line", "inside", "", "inside", "line", "inside", "inside"]
    r = S.line_boundary_mi(df, shuffles=5).set_index("kind")
    assert r.loc["inside", "n_pairs"] == 6 and r.loc["line", "n_pairs"] == 2 and r.loc["paragraph", "n_pairs"] == 0
    # inside the line 'b' is always followed by 'b...' word start? edges: ab->ba gives b|b, ba->ab gives a|a
    assert r.loc["inside", "mi_obs"] > 0.9
    assert r.loc["inside", "ratio_to_inside"] == 1.0


def test_boundary_mi_detects_line_independence():
    rng = np.random.default_rng(1)
    vocab = ["xa", "ay", "yb", "bx"]
    lines = []
    for _ in range(300):  # inside a line: each word starts with the last letter of the previous one
        w = [vocab[rng.integers(4)]]
        for _ in range(4):
            w.append(next(v for v in rng.permutation(vocab) if v[0] == w[-1][-1]))
        lines.append(w)
    df = _paragraphs(lines, 10)
    r = S.line_boundary_mi(df, shuffles=10, seed=2).set_index("kind")
    assert r.loc["inside", "excess_bits"] > 0.8
    assert abs(r.loc["line", "excess_bits"]) < 0.05
    assert abs(r.loc["line", "ratio_to_inside"]) < 0.1
    # custom labels: the same token class on both sides
    lab = df["token"].str[0].tolist()
    r2 = S.line_boundary_mi(df, left=lab, right=lab, shuffles=3)
    assert list(r2["kind"]) == ["inside", "line", "paragraph"]


def test_vertical_similarity_copy_from_above():
    rng = np.random.default_rng(0)
    vocab = ["chol", "daiin", "qokeedy", "shedy", "otar", "lkam", "sor", "cthy", "pchey", "fal"]
    lines = [rng.choice(vocab, 5).tolist()]
    for _ in range(199):  # each token copies the token above with one edit (one letter appended)
        lines.append([t[:4] + "e" for t in lines[-1]] if len(lines) % 10 else rng.choice(vocab, 5).tolist())
    df = _paragraphs(lines, 10)
    v = S.vertical_similarity(df, k=1, shuffles=20)
    assert v["n_line_pairs"] == 180
    assert v["above"] > 0.9 and v["above"] > v["else_above"] and v["excess_above"] > 0.3
    assert v["p_perm"] < 0.1
    # independent lines: no excess
    ind = _paragraphs([rng.choice(vocab, 5).tolist() for _ in range(200)], 10)
    w = S.vertical_similarity(ind, k=0, shuffles=20, seed=1)
    assert abs(w["excess_above"]) < 0.05
    assert abs(w["perm_mean"] - w["expected_above"]) < 0.02
    assert 0 <= w["left_same"] <= 1


def test_vertical_similarity_empty():
    df = make_tokens([["a", "b"]])
    v = S.vertical_similarity(df, shuffles=3)
    assert v["n_tokens"] == 0 and v["above"] != v["above"]
