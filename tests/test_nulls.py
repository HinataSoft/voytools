"""Null models: shuffles and the auto-copy generator."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
import pytest

from voytools import nulls as N
from voytools import stats as S
from voytools._common import line_ids
from voytools.glyphs import collapse_eva
from voytools.morphology import decompose_table
from voytools.table import TOKEN_COLUMNS
from tests.conftest import make_tokens


def _corpus(n_lines=40, seed=0):
    rng = np.random.default_rng(seed)
    vocab = ["qokeedy", "chedy", "shedy", "daiin", "okaiin", "qokain", "ol", "or", "ar", "chol", "dar", "s"]
    return [rng.choice(vocab, size=int(rng.integers(3, 12))).tolist() for _ in range(n_lines)]


@pytest.fixture(scope="module")
def df():
    return make_tokens(_corpus(), page_size=8)


# ------------------------------------------------------------------ shuffles


def test_shuffle_within_lines_preserves_line_multisets_and_positions(df):
    sh = N.shuffle_within_lines(df, seed=1)
    assert len(sh) == len(df) and list(sh.columns) == list(df.columns)
    for _, a in df.groupby("locus_id", sort=False):
        b = sh.loc[a.index]
        assert Counter(a["token"]) == Counter(b["token"])
    positional = [c for c in df.columns if c not in N._TOKEN_PROPERTY_COLUMNS]
    pd.testing.assert_frame_equal(df[positional], sh[positional])
    assert not (sh["token"].to_numpy() == df["token"].to_numpy()).all()  # something actually moved


def test_shuffle_is_deterministic_and_seed_dependent(df):
    a = N.shuffle_within_lines(df, seed=3)
    b = N.shuffle_within_lines(df, seed=3)
    c = N.shuffle_within_lines(df, seed=4)
    pd.testing.assert_frame_equal(a, b)
    assert (a["token"].to_numpy() != c["token"].to_numpy()).any()


def test_shuffle_keeps_token_properties_consistent(df):
    dfd = decompose_table(df, prefixes=["qo", "o", "ch"], suffixes=["edy", "dy", "aiin"])
    sh = N.shuffle_within_lines(dfd, seed=5)
    assert (sh["len_chars"] == sh["token"].str.len()).all()
    assert (sh["len_glyphs"] == [len(collapse_eva(t, "basic")) for t in sh["token"]]).all()
    assert ((sh["prefix"] + sh["stem"] + sh["suffix"]) == sh["token"]).all()
    assert (sh["token_raw"] == sh["token"]).all()
    assert sh.attrs["prefixes"] == dfd.attrs["prefixes"] and sh.attrs["suffixes"] == dfd.attrs["suffixes"]


def test_shuffle_does_not_modify_its_input(df):
    before = df.copy()
    N.shuffle_within_lines(df, seed=0)
    N.shuffle_within_pages(df, seed=0)
    pd.testing.assert_frame_equal(df, before)


def test_shuffle_within_pages_moves_tokens_across_lines_only_inside_a_page(df):
    sh = N.shuffle_within_pages(df, seed=2)
    for _, a in df.groupby("folio", sort=False):
        assert Counter(a["token"]) == Counter(sh.loc[a.index, "token"])
    assert (sh["n_tokens"] == df["n_tokens"]).all()  # line lengths unchanged
    moved = 0
    for _, a in df.groupby("locus_id", sort=False):
        if Counter(a["token"]) != Counter(sh.loc[a.index, "token"]):
            moved += 1
    assert moved > 0


def test_shuffle_destroys_order_dependence_but_not_bag_of_words(df):
    big = make_tokens([["ab", "cd", "ab", "cd", "ab", "cd", "ab", "cd"]] * 40)
    sh = N.shuffle_within_lines(big, seed=0)
    assert S.vocabulary_summary(big).iloc[0].tolist() == S.vocabulary_summary(sh).iloc[0].tolist()
    mi = S.adjacent_token_mi(big, shuffles=0)
    mi_sh = S.adjacent_token_mi(sh, shuffles=0)
    assert mi.iloc[0]["mi_obs"] == pytest.approx(1.0, abs=0.05) and mi_sh.iloc[0]["mi_obs"] < 0.2


def test_shuffle_of_empty_table_is_empty(df):
    empty = df.iloc[:0]
    assert len(N.shuffle_within_lines(empty)) == 0 and len(N.shuffle_within_pages(empty)) == 0


# ------------------------------------------------------------------ autocopy_generate


SEED = ["qokeedy", "chedy", "shedy", "daiin", "okaiin", "qokain", "ol", "or", "ar", "chol", "dar", "s"]


def _flat(lines):
    return [t for l in lines for t in l]


def test_autocopy_generate_is_deterministic_and_has_the_requested_size():
    a = N.autocopy_generate(SEED, 500, seed=7)
    b = N.autocopy_generate(SEED, 500, seed=7)
    c = N.autocopy_generate(SEED, 500, seed=8)
    assert a == b and a != c
    assert len(_flat(a)) == 500 and all(len(l) >= 1 for l in a) and all(isinstance(t, str) and t for t in _flat(a))


def test_autocopy_exact_layout_and_sampled_layout():
    lengths = [3, 5, 8]
    exact = N.autocopy_generate(SEED, 64, line_lengths=lengths, layout="exact", seed=0)
    assert [len(l) for l in exact[:6]] == [3, 5, 8, 3, 5, 8]
    assert len(_flat(exact)) == 64
    sampled = N.autocopy_generate(SEED, 200, line_lengths=lengths, layout="sample", seed=0)
    assert {len(l) for l in sampled[:-1]} <= set(lengths)
    with pytest.raises(ValueError):
        N.autocopy_generate(SEED, 10, layout="diagonal")


def test_autocopy_pure_copying_only_repeats_what_it_has_seen():
    out = _flat(N.autocopy_generate(SEED, 300, p_copy_exact=1.0, p_new=0.0, seed=1))
    assert set(out) <= set(SEED)
    last = _flat(N.autocopy_generate(["x", "y", "ar"], 50, window=1, p_copy_exact=1.0, p_new=0.0, seed=1))
    assert set(last) == {"ar"}  # window=1: always a copy of the previous token


def test_autocopy_edits_respect_the_bigram_model_and_the_length_limit():
    ref = SEED
    pairs = set()
    for t in ref:
        g = collapse_eva(t, "basic")
        pairs.update(zip(("^",) + g, g + ("$",)))
    longest = max(len(collapse_eva(t, "basic")) for t in ref)
    out = _flat(N.autocopy_generate(ref, 2000, p_copy_exact=0.1, p_two_edits=0.5, seed=3))
    assert len(set(out) - set(ref)) > 20  # it does produce new word forms
    for t in out:
        g = collapse_eva(t, "basic")
        assert len(g) <= longest
        assert set(zip(("^",) + g, g + ("$",))) <= pairs  # every glyph pair is attested


def test_autocopy_allowed_from_extends_the_model_and_max_len_caps_words():
    out = _flat(N.autocopy_generate(["ar", "or"], 400, p_copy_exact=0.0, allowed_from=SEED, max_len=4, seed=2))
    assert max(len(collapse_eva(t, "basic")) for t in out) <= 4
    # with only the 2-glyph seed words as the model nothing can grow beyond 2 glyphs
    small = _flat(N.autocopy_generate(["ar", "or"], 400, p_copy_exact=0.0, seed=2))
    assert max(len(t) for t in small) <= 2


def test_autocopy_input_validation():
    with pytest.raises(ValueError):
        N.autocopy_generate([], 10)
    with pytest.raises(ValueError):
        N.autocopy_generate(SEED, 10, edit_type_probs=(0, 0, 0))
    with pytest.raises(ValueError):
        N.autocopy_generate(SEED, 10, edit_type_probs=(1, 2))
    assert N.autocopy_generate(SEED, 0) == []


def test_autocopy_unconstrained_and_recency_variants_run():
    a = _flat(N.autocopy_generate(SEED, 200, constrain_bigrams=False, recency=0.5, seed=4))
    assert len(a) == 200
    b = _flat(N.autocopy_generate(list("abcab"), 100, preset=None, seed=4))
    assert len(b) == 100


# ------------------------------------------------------------------ as_table / autocopy_like


def test_as_table_has_the_token_table_columns_and_works_with_stats():
    lines = N.autocopy_generate(SEED, 600, seed=0)
    t = N.as_table(lines, lines_per_page=10)
    assert list(t.columns) == list(TOKEN_COLUMNS)
    assert len(t) == 600 and t.attrs["alphabet"] == "Eva-" and t.attrs["glyph_preset"] == "basic"
    assert t["token_index"].max() < t["n_tokens"].max() and (t["pos_from_end"] >= 0).all()
    assert t.groupby("locus_id")["n_tokens"].first().sum() == 600
    assert (t["folio"].nunique() == int(np.ceil(len(lines) / 10)))
    assert (t["generic_type"] == "P").all() and (t["locator"].isin(["@", "+"])).all()
    # downstream functions accept it
    assert S.vocabulary_summary(t).iloc[0]["n_tokens"] == 600
    assert len(S.glyph_entropy_report(t, "basic")) == 1
    assert len(S.positional_profile(t, which="first").n) == 6
    assert S.adjacent_token_mi(t, shuffles=2).iloc[0]["n_pairs"] > 0
    assert S.neighbor_similarity(t, shuffles=2) is not None


def test_as_table_with_template_reuses_metadata_and_requires_equal_size(df):
    lengths = np.bincount(line_ids(df)).tolist()
    lines = N.autocopy_generate(SEED, len(df), line_lengths=lengths, layout="exact", seed=0)
    t = N.as_table(lines, df)
    meta = [c for c in df.columns if c not in N._TOKEN_PROPERTY_COLUMNS]
    pd.testing.assert_frame_equal(t[meta], df[meta])
    assert (t["len_chars"] == t["token"].str.len()).all() and t.attrs["source"] == "autocopy"
    with pytest.raises(ValueError):
        N.as_table(lines[:-1], df)


def test_autocopy_like_matches_the_template_layout(df):
    t = N.autocopy_like(df, seed_tokens=SEED, window=6, seed=1)
    assert len(t) == len(df)
    assert (t["n_tokens"] == df["n_tokens"]).all() and (t["locus_id"] == df["locus_id"]).all()
    t2 = N.autocopy_like(df, seed_tokens=SEED, window=6, seed=1)
    pd.testing.assert_frame_equal(t, t2)
    t3 = N.autocopy_like(df, seed_tokens=SEED, window=6, seed=2)
    assert (t3["token"].to_numpy() != t["token"].to_numpy()).any()
    # default seed: the first tokens of the template
    t4 = N.autocopy_like(df, seed=0)
    assert len(t4) == len(df)


# ------------------------------------------------------------------ shuffle_within_position_class


def _par_corpus(n_pars=12, seed=0):
    """Lines with paragraph markup: paragraphs of 5 lines (two per page), 2-10 tokens per line."""
    rng = np.random.default_rng(seed)
    vocab = ["qokeedy", "chedy", "shedy", "daiin", "okaiin", "qokain", "ol", "or", "ar", "chol", "dar", "s"]
    lines = []
    for _ in range(n_pars):
        n_lines = 5
        for k in range(n_lines):
            toks = rng.choice(vocab, size=int(rng.integers(2, 11))).tolist()
            if k == 0:
                toks[0] = "<%>" + toks[0]
            if k == n_lines - 1:
                toks[-1] = toks[-1] + "<$>"
            lines.append(toks)
    return lines


def _default_class(t):
    place = np.where(t["is_first"], "first", np.where(t["is_last"], "last", "mid"))
    rel = np.clip((t["rel_pos"].to_numpy(float) * 5).astype(int), 0, 4)
    lip = np.minimum(t["line_in_par"].to_numpy(float), 2)
    return [f"{a}{b}{c}" for a, b, c in zip(place, rel, lip)]


def test_shuffle_within_position_class_keeps_class_multisets_per_paragraph():
    t = make_tokens(_par_corpus(), page_size=10)
    assert t["par_index"].notna().all()
    sh = N.shuffle_within_position_class(t, seed=1)
    positional = [c for c in t.columns if c not in N._TOKEN_PROPERTY_COLUMNS]
    pd.testing.assert_frame_equal(t[positional], sh[positional])
    key = pd.Series(_default_class(t), index=t.index)
    for _, a in t.groupby([t["folio"], t["par_index"], key], sort=False):
        assert Counter(a["token"]) == Counter(sh.loc[a.index, "token"])
    # line-first and line-last tokens stay first/last; tokens move across lines of a paragraph
    assert set(sh.loc[t["is_first"], "token"]) <= set(t.loc[t["is_first"], "token"])
    moved_lines = sum(Counter(a["token"]) != Counter(sh.loc[a.index, "token"]) for _, a in t.groupby("locus_id"))
    assert moved_lines > 0
    assert (sh["token"].to_numpy() != t["token"].to_numpy()).any()


def test_shuffle_within_position_class_seed_scope_and_custom_classes():
    t = make_tokens(_par_corpus(seed=2), page_size=10)
    a = N.shuffle_within_position_class(t, seed=3)
    pd.testing.assert_frame_equal(a, N.shuffle_within_position_class(t, seed=3))
    assert (a["token"].to_numpy() != N.shuffle_within_position_class(t, seed=4)["token"].to_numpy()).any()
    # one class per row: nothing can move
    same = N.shuffle_within_position_class(t, seed=0, classes=np.arange(len(t)))
    assert (same["token"].to_numpy() == t["token"].to_numpy()).all()
    # a single class with scope 'page' is the page shuffle's multiset per page
    pg = N.shuffle_within_position_class(t, seed=0, scope="page", classes=["x"] * len(t))
    for _, g in t.groupby("folio"):
        assert Counter(g["token"]) == Counter(pg.loc[g.index, "token"])
    # without paragraph markup every line is its own unit
    nopar = make_tokens(_corpus(), page_size=8)
    s2 = N.shuffle_within_position_class(nopar, seed=0, rel_bins=1)
    for _, g in nopar.groupby("locus_id"):
        assert Counter(g["token"]) == Counter(s2.loc[g.index, "token"])
    with pytest.raises(ValueError):
        N.shuffle_within_position_class(t, scope="quire")
    with pytest.raises(ValueError):
        N.shuffle_within_position_class(t, classes=["x"])
    assert len(N.shuffle_within_position_class(t.iloc[:0])) == 0
