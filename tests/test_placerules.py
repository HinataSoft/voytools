"""Place rules: item decomposition of MI, position chain, seam test, stability, post-processor R1-R3."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from voytools import placerules as PR
from voytools._common import line_ids
from voytools.morphology import decompose_table
from voytools.variation import position_class
from tests.conftest import make_tokens

PRE = ["qo", "o", "ch", "d", "y", "s", "p"]
SUF = ["edy", "dy", "aiin", "ol", "or", "m"]  # no "am": a rewritten "ka"+"m" would re-decompose as k+am


def _corpus(n_lines=120, seed=0, rule=True):
    """Lines of 3..9 tokens; with ``rule`` the first token starts with d-/y-, the last ends in -m."""
    rng = np.random.default_rng(seed)
    stems = ["ke", "ka", "te", "ko", "ai", "ol", "ee", "ka"]
    lines = []
    for _ in range(n_lines):
        n = int(rng.integers(3, 10))
        toks = []
        for i in range(n):
            p = rng.choice(PRE)
            s = rng.choice(SUF[:5])
            if rule and i == 0:
                p = rng.choice(["d", "y"], p=[0.7, 0.3])
            if rule and i == n - 1 and rng.random() < 0.6:
                s = "m"
            toks.append(p + rng.choice(stems) + s)
        lines.append(toks)
    return lines


@pytest.fixture(scope="module")
def df():
    return make_tokens(_corpus(), page_size=10)


@pytest.fixture(scope="module")
def flat():
    return make_tokens(_corpus(seed=5, rule=False), page_size=10)


@pytest.fixture(scope="module")
def labels(df):
    return PR.token_labels(df, PRE, SUF, "basic")


# ------------------------------------------------------------------ decomposition


def test_pointwise_terms_sum_to_mi_and_marginals_are_nonnegative(df, labels):
    pos = position_class(df)
    d = PR.item_decomposition(labels["prefix"], pos, "prefix", "pos")
    assert d["mi"] == pytest.approx(PR.mutual_information(labels["prefix"], pos), abs=1e-9)
    assert d["pairs"]["contribution"].sum() == pytest.approx(d["mi"], abs=1e-9)
    assert (d["by_x"]["contribution"] >= 0).all() and (d["by_y"]["contribution"] >= 0).all()
    assert d["by_x"]["contribution"].sum() == pytest.approx(d["mi"], abs=1e-9)
    assert d["by_x"]["cum_share"].iloc[-1] == pytest.approx(1.0)
    assert list(d["coverage"]["items"]) == ["pair", "prefix", "pos"]
    # the rule d-/y- at line start makes a few items carry most of the information
    cov = d["coverage"].set_index("items").loc["prefix"]
    assert cov["n_50"] <= cov["n_80"] <= cov["n_95"] <= cov["n_items"]
    assert cov["n_80"] <= 5


def test_coverage_counts_items_and_handles_negatives():
    c = PR.coverage([0.5, 0.3, 0.2])
    assert (c["n_50"], c["n_80"], c["n_95"]) == (1, 2, 3)
    c = PR.coverage([0.6, 0.6, -0.2])  # total 1.0; the positive terms reach every level first
    assert c["n_negative"] == 1 and c["n_50"] == 1 and c["n_80"] == 2 and c["n_95"] == 2
    assert np.isnan(PR.coverage([])["n_80"])


def test_mutual_information_matches_sklearn_free_formula():
    x = ["a", "a", "b", "b"]
    y = ["u", "u", "v", "v"]
    assert PR.mutual_information(x, y) == pytest.approx(1.0)
    assert PR.mutual_information(x, ["u", "v", "u", "v"]) == pytest.approx(0.0)
    assert PR.mutual_information(x, y, z=y) == pytest.approx(0.0)


def test_position_chain_terms_add_up(df, labels):
    ch = PR.position_chain(df, labels["prefix"], shuffles=0)
    assert ch["I_pos"] == pytest.approx(ch["I_place"] + ch["I_lip_given_place"] + ch["I_rel_given_place_lip"], abs=1e-9)
    assert ch["I_place"] > 0.1  # the d-/y- rule sits on the first token
    chs = PR.position_chain(df, labels["prefix"], shuffles=5, seed=1)
    assert chs["I_pos"] == pytest.approx(chs["I_pos_raw"] - chs["I_pos_bias"], abs=1e-12)
    assert chs["I_pos"] < chs["I_pos_raw"]


def test_seam_pairs_and_glyph_vs_affix(df, labels):
    L, R, idx = PR.seam_pairs(df, labels["last_glyph"], labels["first_glyph"])
    line = line_ids(df)
    assert len(L) == len(R) == (line[:-1] == line[1:]).sum()
    assert (line[idx] == line[idx + 1]).all()
    r = PR.seam_glyph_vs_affix(df, labels, shuffles=5, seed=0)
    assert r["n_pairs"] == len(L)
    for k in ("I_glyph", "I_affix", "I_affix_given_glyph", "I_glyph_given_affix"):
        assert np.isfinite(r[k]) and np.isfinite(r[f"{k}_raw"])
    # affixes here are random given the neighbour: the corrected terms are near zero
    assert abs(r["I_affix"]) < 0.05


# ------------------------------------------------------------------ stability


def test_conditional_jsd_zero_for_identical_tables():
    x = ["f"] * 50 + ["m"] * 50
    y = ["d", "y"] * 25 + ["o", "ch"] * 25
    r = PR.conditional_jsd(x, y, x, y, min_count=10)
    assert r["jsd"] == pytest.approx(0.0, abs=1e-12) and r["n_contexts"] == 2
    r2 = PR.conditional_jsd(x, y, x, ["o", "ch"] * 25 + ["d", "y"] * 25, min_count=10)
    assert r2["jsd"] > 0.5


def test_stability_detects_shared_and_different_rules(df, labels):
    place = np.where(df["is_first"].to_numpy(dtype=bool), "first", "other")
    groups = np.where(np.arange(len(df)) % 2 == 0, "A", "B")  # arbitrary split: the same rule in both
    st = PR.stability(place, labels["prefix"], groups, line_ids(df), n_rep=20, min_count=10, contexts=["first"])
    assert len(st) == 1
    r = st.iloc[0]
    assert r["n_contexts"] == 1 and r["n_rules"] >= 1
    assert r["ratio"] < 2.0 and r["share_shared"] >= 0.5
    # different rule in B: first prefixes replaced by s-
    pre = labels["prefix"].copy()
    mb = (groups == "B") & (place == "first")
    pre[mb] = "s"
    st2 = PR.stability(place, pre, groups, line_ids(df), n_rep=20, min_count=10, contexts=["first"])
    assert st2.iloc[0]["jsd_between"] > r["jsd_between"] and st2.iloc[0]["ratio"] > 3


# ------------------------------------------------------------------ post-processor


def test_learn_rules_tables(df):
    rules = PR.learn_place_rules(df, PRE, SUF, "basic")
    assert set(rules.r1) <= {"first|pf", "first|other", "inner|pf"}
    t = rules.table("R1")
    assert set(t.columns) == {"rule", "context", "item", "p"}
    assert t.groupby("context")["p"].sum().round(9).eq(1.0).all()
    assert rules.r2.sum() == pytest.approx(1.0) and rules.r2["m"] > 0.4
    assert all(s.sum() == pytest.approx(1.0) for s in rules.r3.values())
    with pytest.raises(ValueError):
        rules.table("R9")


def test_apply_rules_keeps_layout_and_stems_and_moves_distributions(df, flat):
    rules = PR.learn_place_rules(df, PRE, SUF, "basic")
    out = PR.apply_place_rules(flat, rules, ("R1", "R2"), seed=3)
    assert len(out) == len(flat) and list(out.columns) == list(flat.columns)
    positional = [c for c in flat.columns if c not in ("token", "token_raw", "len_chars", "len_glyphs")]
    pd.testing.assert_frame_equal(flat[positional], out[positional])
    assert (out["len_chars"] == out["token"].str.len()).all()
    assert out.attrs["place_rules"] == ["R1", "R2"] and out.attrs["source"].endswith("+R1+R2")
    d0 = decompose_table(flat, prefixes=PRE, suffixes=SUF)
    d1 = decompose_table(out, prefixes=PRE, suffixes=SUF)
    assert (d0["stem"] == d1["stem"]).all()
    first = flat["is_first"].to_numpy(dtype=bool)
    last = flat["is_last"].to_numpy(dtype=bool)
    # untouched: suffixes of non-last tokens and prefixes of inner tokens of non-paragraph-first lines
    inner = ~first & (pd.to_numeric(flat["line_in_par"], errors="coerce").fillna(2).to_numpy() > 0)
    assert (d0.loc[~last, "suffix"] == d1.loc[~last, "suffix"]).all()
    assert (d0.loc[inner, "prefix"] == d1.loc[inner, "prefix"]).all()
    # moved: the line-first prefixes are now mostly d-/y-, the line-last suffixes mostly -m
    share_dy = d1.loc[first, "prefix"].isin(["d", "y"]).mean()
    assert share_dy > 0.8 and d0.loc[first, "prefix"].isin(["d", "y"]).mean() < 0.5
    assert d1.loc[last, "suffix"].eq("m").mean() > 0.4 and d0.loc[last, "suffix"].eq("m").mean() < 0.05
    # minimal rewrite: tokens that already had an allowed prefix mostly kept it
    had = first & d0["prefix"].isin(["d", "y"]).to_numpy()
    assert (d0.loc[had, "prefix"] == d1.loc[had, "prefix"]).mean() > 0.9


def test_apply_rules_r3_conditions_on_previous_last_glyph(df, flat):
    rules = PR.learn_place_rules(df, PRE, SUF, "basic", min_context=5)
    out = PR.apply_place_rules(flat, rules, ("R3",), seed=0)
    lab0 = PR.token_labels(flat, PRE, SUF, "basic")
    lab1 = PR.token_labels(out, PRE, SUF, "basic")
    assert (lab0["last_glyph"] == lab1["last_glyph"]).all()  # stems and suffixes untouched
    L, R, _ = PR.seam_pairs(out, lab1["last_glyph"], lab1["prefix"])
    L0, R0, _ = PR.seam_pairs(flat, lab0["last_glyph"], lab0["prefix"])
    assert PR.mutual_information(L, R) >= PR.mutual_information(L0, R0)
    assert out.attrs["place_rules_changed"]["R3"] > 0
    a = PR.apply_place_rules(flat, rules, ("R3",), seed=0)
    b = PR.apply_place_rules(flat, rules, ("R3",), seed=1)
    assert (a["token"] == out["token"]).all() and (a["token"] != b["token"]).any()
    with pytest.raises(ValueError):
        PR.apply_place_rules(flat, rules, ("R4",))
