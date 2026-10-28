"""voytools.paragraphs: line sharing by distance and paragraph/page boundary (experiment 020)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tests.conftest import make_tokens
from voytools import paragraphs as P


def _table() -> pd.DataFrame:
    """Two pages; page 1 has paragraphs of 2 + 2 lines, page 2 one paragraph of 2 lines."""
    lines = [["a", "b", "c"], ["a", "d"], ["e", "f"], ["e", "g", "h"], ["a", "e"], ["b", "x"]]
    df = make_tokens(lines, page_size=4)
    # paragraphs: page 1 lines 0-1 -> par 0, lines 2-3 -> par 1
    par = np.array([0, 0, 1, 1, 0, 0])
    lid = df.groupby("locus_id", sort=False).ngroup().to_numpy()
    df["par_index"] = par[lid]
    df["line_in_par"] = np.array([0, 1, 0, 1, 0, 1])[lid]
    df["stem"] = df["token"].str.rstrip("y")
    return df


def test_strip_gallows_and_unit_ids():
    assert P.strip_gallows(["kchol", "cthey", "daiin", "pol"]) == ["chol", "ey", "daiin", "ol"]
    assert P.strip_gallows(["gab", "Kab", "ab"], "v101") == ["ab", "ab", "ab"]
    ids, vocab = P.unit_ids(["b", "a", "b"])
    assert ids.tolist() == [0, 1, 0] and vocab == ["b", "a"]


def test_line_table_structure_and_filters():
    L = P.line_table(_table(), "stem", min_count=1)
    assert len(L) == 6
    assert L["par"].tolist() == [0, 0, 1, 1, 0, 0]
    assert L["lip"].tolist() == [0, 1, 0, 1, 0, 1]
    assert L["par_len"].tolist() == [2, 2, 2, 2, 2, 2]
    assert L["n_pars"].tolist() == [2, 2, 2, 2, 1, 1]
    assert L["n"].tolist() == [3, 2, 2, 3, 2, 2]
    L2 = P.line_table(_table(), "stem", min_count=2)  # c, d, f, g, h, x are hapax
    assert L2["n"].tolist() == [2, 1, 1, 1, 2, 1]
    assert L2.attrs["n_units"] == 3  # a, b, e
    L3 = P.line_table(_table(), "stem", min_count=2, drop_top=1)  # a (3 occurrences) dropped
    assert L3.attrs["n_units"] == 2 and L3["n"].tolist() == [1, 0, 1, 1, 1, 1]


def test_pair_content_and_layout():
    L = P.line_table(_table(), "stem", min_count=1)
    pc = P.pair_content(L)
    assert len(pc) == 6 + 1  # C(4,2) on page 1, C(2,2) on page 2
    pairs = P.pair_layout(pc, L)
    r = pairs[(pairs["i"] == 0) & (pairs["j"] == 1)].iloc[0]
    assert r["obs"] == 1 and r["d"] == 1 and r["same"]
    r = pairs[(pairs["i"] == 1) & (pairs["j"] == 2)].iloc[0]
    assert r["obs"] == 0 and r["d"] == 1 and not r["same"]
    # expectation: page 1 has 10 tokens, a:2 b:1 c:1 d:1 e:2 f:1 g:1 h:1; lines 0 (3 tokens) and 1 (2 tokens)
    p = np.array([2, 1, 1, 1, 2, 1, 1, 1]) / 10
    e = np.sum((1 - (1 - p) ** 3) * (1 - (1 - p) ** 2))
    assert pairs[(pairs["i"] == 0) & (pairs["j"] == 1)]["e_page"].iloc[0] == pytest.approx(e)
    assert (pairs["e_glob"] > 0).all()


def test_structure_and_nulls():
    par, lip, plen = P.structure_from_lengths([2, 1])
    assert par.tolist() == [0, 0, 1] and lip.tolist() == [0, 1, 0] and plen.tolist() == [2, 2, 1]
    L = P.line_table(_table(), "stem", min_count=1)
    Lb = P.permute_lines(L, 3)
    assert sorted(Lb[Lb["folio"] == "f1r"]["pos"].tolist()) == [0, 1, 2, 3]
    assert set(Lb["units"].apply(tuple)) == set(L["units"].apply(tuple))  # contents unchanged
    # the structure by position is the same as before
    for f in ("f1r", "f2r"):
        a = L[L["folio"] == f].sort_values("pos")[["par", "lip"]].to_numpy()
        b = Lb[Lb["folio"] == f].sort_values("pos")[["par", "lip"]].to_numpy()
        assert (a == b).all()
    La = P.shuffle_boundaries(L, 0)
    assert P.paragraph_lengths(La)["f1r"] == [2, 2] and La["pos"].tolist() == L["pos"].tolist()


def test_metrics_small():
    L = P.line_table(_table(), "stem", min_count=1)
    pairs = P.pairs_of(L)
    t = P.sharing_by_distance(pairs, max_d=3, min_pars=2)
    assert len(t) == 6
    d1 = t[t["d"] == 1].set_index("same")
    assert d1.loc[True, "n"] == 2 and d1.loc[False, "n"] == 1
    assert d1.loc[True, "obs"] == 2 and d1.loc[False, "obs"] == 0
    step = P.boundary_step(t, 2)
    assert step[1] > 0 and np.isnan(step[2])  # no cross pair at d = 2 on a 2+2 page
    dec, slope = P.within_par_decay(pairs, max_d=2)
    assert dec["n"].tolist()[0] == 3 and np.isnan(slope)  # only d = 1 exists
    fl = P.first_line_sharing(pairs, max_d=2)
    assert set(fl["kind"]) == {"first_own", "first_next", "mid_own", "mid_next"}
    assert fl[(fl["kind"] == "first_own") & (fl["d"] == 1)]["n"].iloc[0] == 3
    assert fl[(fl["kind"] == "first_next") & (fl["d"] == 2)]["n"].iloc[0] == 1
    pb = P.page_boundary(L, pairs, max_d=2, open_end={"f1r": True})
    r = pb[(pb["group"] == "across_all") & (pb["d"] == 1)].iloc[0]
    assert r["n"] == 1 and r["obs"] == 1  # line 3 (e g h) and line 4 (a e) share e
    assert pb[(pb["group"] == "across_open") & (pb["d"] == 1)]["n"].iloc[0] == 1
    fb = P.form_by_boundary(pairs, max_d=1)
    assert fb[(fb["d"] == 1) & fb["same"]]["n_occ"].iloc[0] == 2
    assert fb[(fb["d"] == 0) & fb["same"]]["ident"].iloc[0] == 1.0
    m = P.all_metrics(L, pairs, open_end={"f1r": True})
    assert m["share_p__d1_same"] == pytest.approx(d1.loc[True, "s"], abs=1e-5)
    nf = P.sharing_by_distance(pairs, 3, 2, exclude_first=True)
    assert nf["n"].fillna(0).sum() == 1 and "stepnf_p__d2" in m  # only lines 1 and 3 (both second lines) remain
    assert "page__across_open_d1" in m and "form__excess_same_d0" in m


def test_weights_bootstrap_over_pages():
    L = P.line_table(_table(), "stem", min_count=1)
    pairs = P.pairs_of(L)
    full = P.sharing_by_distance(pairs, 3, 1)
    w = P.sharing_by_distance(pairs, 3, 1, weights={"f1r": 2.0, "f2r": 0.0})
    only1 = P.sharing_by_distance(pairs[pairs["folio"] == "f1r"], 3, 1)
    assert w["s"].fillna(-1).tolist() == pytest.approx(only1["s"].fillna(-1).tolist())
    assert w["n"].sum() == 2 * only1["n"].sum() and full["n"].sum() == 7


def test_permuted_lines_have_no_step_on_average():
    """Sanity check of the expectation: on permuted lines the boundary step averages to ~0."""
    rng = np.random.default_rng(1)
    vocab = [f"w{i}" for i in range(30)]
    lines = [rng.choice(vocab, size=6, p=np.linspace(2, 1, 30) / np.linspace(2, 1, 30).sum()).tolist() for _ in range(80)]
    df = make_tokens(lines, page_size=8)
    lid = df.groupby("locus_id", sort=False).ngroup().to_numpy()
    df["par_index"] = ((lid % 8) // 3)[...]  # paragraphs of 3, 3, 2 lines
    df["stem"] = df["token"]
    L = P.line_table(df, "stem", min_count=1)
    pc = P.pair_content(L)
    steps = []
    for s in range(20):
        Lb = P.permute_lines(L, s)
        steps.append(P.boundary_step(P.sharing_by_distance(P.pair_layout(pc, Lb), 4, 2), 1)[1])
    assert abs(np.mean(steps)) < 0.1
