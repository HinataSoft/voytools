"""Within-word variation: entropy decomposition and residual structure (voytools.variation)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from voytools import variation as V


def synthetic(n: int = 6000, seed: int = 0) -> pd.DataFrame:
    """Stems with 4 suffixes: suffix = f(stem, pos) for 'det' rows, random for 'rnd' rows.

    ``pos`` takes 4 values, ``noise`` 50 values unrelated to the suffix, ``page`` 30 values.
    """
    rng = np.random.default_rng(seed)
    stem = rng.choice(["a", "b", "c", "d", "e"], size=n)
    pos = rng.integers(0, 4, size=n)
    det = (np.searchsorted(["a", "b", "c", "d", "e"], stem) + pos) % 4
    rnd = rng.integers(0, 4, size=n)
    return pd.DataFrame({"stem": stem, "pos": pos, "noise": rng.integers(0, 50, size=n),
                         "page": np.sort(rng.integers(0, 30, size=n)),
                         "det": [f"s{v}" for v in det], "rnd": [f"s{v}" for v in rnd]})


def test_codes_and_entropy():
    assert V.codes(["x", "y", "x", None]).tolist() == [0, 1, 0, 2]
    assert V.entropy(np.array([0, 1, 0, 1])) == pytest.approx(1.0)
    assert V.conditional_entropy(np.array([0, 1, 2, 3]), np.array([0, 0, 1, 1])) == pytest.approx(1.0)
    assert V.joint_codes(np.array([0, 1, 0]), np.array([2, 2, 2])).tolist() == [0, 1, 0]


def test_permute_within_keeps_groups():
    g = np.array([0, 0, 1, 1, 1, 2])
    idx = V.permute_within(g, np.random.default_rng(1))
    assert sorted(idx.tolist()) == list(range(6))
    assert (g[idx] == g).all()


def test_decomposition_known_answer():
    df = synthetic()
    out = V.conditional_entropy_decomposition(df, "det", "stem", {"pos": "pos", "noise": "noise"},
                                              n_perm=20, min_count=1, seed=0).set_index("X")
    # suffix fully determined by stem + pos: everything explained by pos, nothing by noise
    assert out.loc["pos", "H_given"] == pytest.approx(2.0, abs=0.01)
    assert out.loc["pos", "explained"] == pytest.approx(1.0, abs=0.03)
    assert abs(out.loc["noise", "explained"]) < 0.03
    # the raw reduction of the 50-valued noise is clearly positive; the permutation removes it
    assert out.loc["noise", "reduction_raw"] > 0.05
    out2 = V.conditional_entropy_decomposition(df, "rnd", "stem", ["pos", "noise", "page"], n_perm=20,
                                               min_count=1, seed=0).set_index("X")
    assert (out2["explained"].abs() < 0.03).all()


def test_crossval_explained():
    df = synthetic(n=20000)
    # seed != data seed: the folds must not reuse the stream that drew the stems
    det = V.crossval_explained(df, "det", "stem", {"pos": "pos", "noise": "noise"}, min_count=1, seed=7)
    rnd = V.crossval_explained(df, "rnd", "stem", ["pos", "noise", "page"], min_count=1, seed=7)
    det = det.set_index("X")
    assert det.loc["pos", "explained_cv"] > 0.9
    assert det.loc["noise", "explained_cv"] < 0.01
    assert (rnd["explained_cv"] < 0.01).all()


def test_decomposition_min_count_and_joint():
    df = synthetic()
    df.loc[:9, "stem"] = "rare"
    out = V.conditional_entropy_decomposition(df, "det", "stem", {"both": ["pos", "noise"]}, n_perm=5,
                                              min_count=20, seed=0)
    assert out.loc[0, "n_stems"] == 5 and out.loc[0, "n_tokens"] == len(df) - 10
    # 1000 (stem, X) cells: the permutation correction is conservative (lower bound), raw is ~1
    assert out.loc[0, "explained_raw"] > 0.95 and out.loc[0, "explained"] > 0.7


def test_neighbour_mi():
    rng = np.random.default_rng(3)
    n = 8000
    stem = rng.choice(["a", "b"], size=n)
    suf = rng.integers(0, 3, size=n)
    copy = suf.copy()
    for i in range(1, n):  # suffix repeats the previous token's suffix with probability 0.9
        if rng.random() < 0.9:
            copy[i] = copy[i - 1]
    line = np.arange(n) // 8
    base = pd.DataFrame({"stem": stem, "cv_line": line, "cv_pos": np.arange(n) % 8})
    agree = V.neighbour_conditional_mi(base.assign(suffix=copy), given=["stem", "cv_pos"], n_perm=20, seed=0)
    indep = V.neighbour_conditional_mi(base.assign(suffix=suf), given=["stem", "cv_pos"], n_perm=20, seed=0)
    assert agree["mi"] > 0.5 and agree["p"] < 0.1
    assert abs(indep["mi"]) < 0.01
    assert agree["n_pairs"] == n - n // 8


def test_mode_consistency():
    rng = np.random.default_rng(4)
    n = 6000
    page = np.repeat(np.arange(200), 30)
    stem = rng.choice(["a", "b", "c"], size=n)
    page_suffix = rng.integers(0, 4, size=200)
    consistent = np.where(rng.random(n) < 0.8, page_suffix[page], rng.integers(0, 4, size=n))
    df = pd.DataFrame({"stem": stem, "cv_page": page})
    c = V.mode_consistency(df.assign(suffix=consistent), n_perm=10, seed=0)
    r = V.mode_consistency(df.assign(suffix=rng.integers(0, 4, size=n)), n_perm=10, seed=0)
    assert c["excess"] > 0.3 and c["z"] > 5
    assert abs(r["excess"]) < 0.03
    assert c["n_stems"] == 3


def test_covariates_and_top_table(zl_tokens_small):
    df = V.add_covariates(zl_tokens_small)
    for col in ("cv_pos", "cv_relpos", "cv_place", "cv_lip", "cv_lang", "cv_hand", "cv_sec", "cv_page",
                "cv_par", "cv_line"):
        assert col in df.columns
    assert set(df["cv_place"]) <= {"first", "mid", "last", "single"}
    dfd = df.assign(stem=df["token"].str[:2], suffix=df["token"].str[-1])
    t = V.affix_by_covariate(dfd, top=3)
    assert t["stem_rank"].max() <= 3
    s = t.groupby(["stem", "cv_place"])["share"].sum()
    assert np.allclose(s.to_numpy(), 1.0)


@pytest.fixture(scope="module")
def zl_tokens_small(synth_tables):
    from voytools import filter_tokens

    return filter_tokens(synth_tables.tokens, generic_type="P").reset_index(drop=True)
