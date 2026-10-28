"""voytools.recurrence: occurrence pairs, distance bins, exact null expectations (experiment 015)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from voytools import recurrence as R


def _table() -> pd.DataFrame:
    """Four pages, stems a (3 occ.), b (2), c (1), d (2 same page); page sizes 3/2/2/3."""
    rows = [
        ("f1r", "H", "A", "1", "A", "1", "a", "xa"),
        ("f1r", "H", "A", "1", "A", "1", "d", "d"),
        ("f1r", "H", "A", "1", "A", "1", "d", "dy"),
        ("f1v", "H", "A", "1", "A", "1", "a", "xa"),
        ("f1v", "H", "A", "1", "A", "1", "c", "c"),
        ("f2r", "B", "B", "2", "B", "2", "b", "b"),
        ("f2r", "B", "B", "2", "B", "2", "a", "ay"),
        ("f2v", "B", "B", "2", "B", "2", "b", "by"),
        ("f2v", "B", "B", "2", "B", "2", "", "y"),
        ("f2v", "B", "B", "2", "B", "2", "", "o"),
    ]
    df = pd.DataFrame(rows, columns=["folio", "illustration", "language", "hand", "quire", "bifolio", "stem", "token"])
    df["generic_type"] = "P"
    return df


@pytest.fixture
def pages():
    return R.page_table(_table())


def test_page_table_rank_size_weight(pages):
    assert pages.index.tolist() == ["f1r", "f1v", "f2r", "f2v"]
    assert pages["rank"].tolist() == [0, 1, 2, 3]
    assert pages["size"].tolist() == [3.0, 2.0, 2.0, 3.0]
    assert pages["weight"].sum() == pytest.approx(1.0)
    assert pages["bifolio"].tolist() == ["A/1", "A/1", "B/2", "B/2"]
    # explicit file order (pages absent from the table are allowed in the rank dict)
    p2 = R.page_table(_table(), {"f1r": 0, "fX": 1, "f1v": 2, "f2r": 3, "f2v": 4})
    assert p2["rank"].tolist() == [0, 2, 3, 4]


def test_bins():
    assert R.distance_bin([0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 100]).tolist() == [
        "0", "1", "2-3", "2-3", "4-7", "4-7", "8-15", "8-15", "16-31", "16-31", "32+", "32+"]
    assert R.form_bin([0, 1, 3, 4, 15, 16, 40]).tolist() == ["0", "1-3", "1-3", "4-15", "4-15", "16+", "16+"]


def test_rare_stems_and_occurrences(pages):
    df = _table()
    assert R.rare_stems(df, 2, 5).to_dict() == {"a": 3, "b": 2, "d": 2}
    assert R.rare_stems(df, 1, 1).to_dict() == {"c": 1}
    assert R.rare_stems(df, 2, 5, within={"b", "zz"}).to_dict() == {"b": 2}
    occ = R.stem_occurrences(df, pages, 2, 5)
    assert len(occ) == 7 and occ["stem"].is_monotonic_increasing
    assert occ.loc[occ["stem"] == "a", "rank"].tolist() == [0, 1, 2]
    assert occ.loc[occ["stem"] == "a", "section"].tolist() == ["H", "H", "B"]


def test_pairs(pages):
    occ = R.stem_occurrences(_table(), pages, 2, 5)
    pairs = R.occurrence_pairs(occ)
    assert len(pairs) == 3 + 1 + 1
    a = pairs[pairs["stem"] == "a"].sort_values("dist")
    assert a["dist"].tolist() == [1, 1, 2]
    # xa/xa identical, xa/ay and xa/ay different
    assert a["same_form"].tolist() == [True, False, False]
    assert a["same_section"].tolist() == [True, False, False]
    d = pairs[pairs["stem"] == "d"].iloc[0]
    assert d["dist"] == 0 and not d["same_form"] and d["same_bifolio"]
    assert R.occurrence_pairs(occ.iloc[:0]).empty


def test_profiles_sum_to_one(pages):
    occ = R.stem_occurrences(_table(), pages, 2, 5)
    obs = R.observed_profile(R.occurrence_pairs(occ))
    assert sum(obs[b] for b in R.DISTANCE_LABELS) == pytest.approx(1.0)
    exp = R.expected_profile(pages)
    assert sum(exp[b] for b in R.DISTANCE_LABELS) == pytest.approx(1.0)
    w = pages["weight"].to_numpy()
    assert exp["0"] == pytest.approx(float((w ** 2).sum()))
    assert exp["far"] == 0.0  # four pages: no pair at distance >= 16


def test_permutation_matches_exact_expectation():
    """Mean of the sampled permutation = exact expectation (profile, same-attribute shares, clusters)."""
    rng = np.random.default_rng(1)
    n_pages = 40
    df = pd.DataFrame({
        "folio": np.repeat([f"p{i}" for i in range(n_pages)], rng.integers(5, 30, size=n_pages)),
    })
    df["illustration"] = np.where(df["folio"].str[1:].astype(int) < 20, "H", "B")
    df["language"] = df["illustration"]
    df["hand"] = "1"
    df["quire"] = (df["folio"].str[1:].astype(int) // 8).astype(str)
    df["bifolio"] = "1"
    stems = [f"s{i}" for i in range(60)]
    df["stem"] = rng.choice(stems, size=len(df))
    df["token"] = df["stem"] + rng.choice(["", "y", "dy"], size=len(df))
    df["generic_type"] = "P"
    pages = R.page_table(df)
    occ = R.stem_occurrences(df, pages, 2, 8)
    es, ep = R.same_by_distance(pages), R.expected_profile(pages)
    runs = pd.DataFrame([R.recurrence_metrics(R.permute_pages(occ, pages, s), pages, es, ep) for s in range(300)])
    m = runs.mean(numeric_only=True)
    for b in R.DISTANCE_LABELS:
        assert m[f"p_{b}"] == pytest.approx(ep[b], abs=0.02)
    assert m["x_far"] == pytest.approx(0.0, abs=0.02)
    assert m["far_section"] == pytest.approx(m["far_section_exp"], abs=0.03)
    assert m["mid_quire"] == pytest.approx(m["mid_quire_exp"], abs=0.03)
    assert m["clq_n"] == pytest.approx(m["clq_exp"], abs=0.1 * m["clq_exp"] + 1)
    assert m["form_same_16+"] == pytest.approx(m["form_exp_16+"], abs=0.03)


def test_form_expectation_exact():
    df = _table()
    pages = R.page_table(df)
    occ = R.stem_occurrences(df, pages, 2, 5)
    pairs = R.occurrence_pairs(occ)
    f = R.form_by_distance(occ, pairs)
    # stem a: forms xa, xa, ay -> p_same = 2*1 / (3*2) = 1/3; b: 0; d: 0 (distance 0)
    assert f["exp_0"] == pytest.approx(0.0)
    assert f["n_1-3"] == 4
    assert f["exp_1-3"] == pytest.approx((3 * (1 / 3) + 0) / 4)
    assert f["same_1-3"] == pytest.approx(1 / 4)


def test_clusters_exact_probability():
    w = np.array([0.5, 0.3, 0.2])
    # P(3 draws hit 3 distinct) = 3! * 0.5*0.3*0.2 = 0.18
    assert 1 - R._prob_at_most_two(w, 3) == pytest.approx(0.18)
    df = _table()
    pages = R.page_table(df)
    occ = R.stem_occurrences(df, pages, 2, 5)
    c = R.far_clusters(occ, pages, by="quire")
    assert c["n_stems"] == 1 and c["n"] == 0  # stem a: quires A, A, B
    assert 0 <= c["exp"] <= 1


def test_label_stems():
    df = pd.DataFrame({"generic_type": ["L", "L", "L", "P"], "subtype": ["z", "0", "f", "0"],
                       "stem": ["a", "b", "", "c"]})
    assert R.label_stems(df) == {"a"}
    assert R.label_stems(df, subtypes=()) == {"a", "b"}


def test_metrics_keys(pages):
    occ = R.stem_occurrences(_table(), pages, 2, 5)
    m = R.recurrence_metrics(occ, pages)
    assert m["n_stems"] == 3 and m["n_pairs"] == 5
    for k in ("p_0", "e_32+", "x_far", "x_far_nl", "far_section_exp", "mid_bifolio", "form_same_0",
              "clq_share", "cls_exp"):
        assert k in m
