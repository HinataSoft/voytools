"""Structured-template null (voytools.templates)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from voytools._common import line_ids
from voytools.morphology import decompose_table, paradigm_metrics, slot_independent_sample
from voytools.templates import ParadigmTemplate
from voytools.wordgrammar import learn_word_grammar
from tests.conftest import make_tokens

PRE = ["qo", "o", "ch", "sh", "d"]
SUF = ["aiin", "edy", "dy", "y", "ol", "ar"]
STEMS = [c + v + e for c in ("k", "t", "p", "f", "ck", "ch") for v in ("e", "ee", "o", "a") for e in ("", "d", "l", "r")]


def _corpus(n_lines=120, seed=0):
    rng = np.random.default_rng(seed)
    lines = []
    for _ in range(n_lines):
        lines.append([str(rng.choice(["", *PRE])) + str(rng.choice(STEMS)) + str(rng.choice(["", *SUF]))
                      for _ in range(int(rng.integers(4, 10)))])
    return lines


@pytest.fixture(scope="module")
def real():
    df = make_tokens(_corpus(), page_size=10)
    # two artificial sections and two languages so that topic / suffix_by have something to work on
    page = df["folio"].astype(str)
    first = sorted(page.unique())[: page.nunique() // 2]
    df["illustration"] = np.where(page.isin(first), "H", "B")
    df["language"] = np.where(page.isin(first), "A", "B")
    return df


@pytest.fixture(scope="module")
def dfd(real):
    return decompose_table(real, prefixes=PRE, suffixes=SUF)


def test_sample_is_deterministic_and_seed_dependent(real, dfd):
    t = ParadigmTemplate.learn(dfd, n_stems=60, topic_concentration=1.0, page_locality=0.3,
                               prefix_dependence="stem")
    a, b, c = t.sample(real, seed=3), t.sample(real, seed=3), t.sample(real, seed=4)
    pd.testing.assert_frame_equal(a, b)
    assert (a["token"].to_numpy() != c["token"].to_numpy()).any()
    # a fresh instance gives the same text (no hidden state besides the stem-pool cache)
    t2 = ParadigmTemplate.learn(dfd, n_stems=60, topic_concentration=1.0, page_locality=0.3,
                                prefix_dependence="stem")
    pd.testing.assert_frame_equal(a, t2.sample(real, seed=3))
    # pool_seed selects the stem inventory only
    d = t.sample(real, seed=4, pool_seed=3)
    assert t.stem_pool(3) == t2.stem_pool(3)
    assert (d["token"].to_numpy() != a["token"].to_numpy()).any()


def test_sample_keeps_layout_and_metadata(real, dfd):
    out = ParadigmTemplate.learn(dfd, n_stems=50).sample(real, seed=1)
    assert len(out) == len(real)
    assert (np.bincount(line_ids(out)) == np.bincount(line_ids(real))).all()
    meta = [c for c in real.columns if c not in ("token", "token_raw", "len_chars", "len_glyphs",
                                                    "has_unreadable", "has_alternative", "has_ligature",
                                                    "has_high_ascii")]
    pd.testing.assert_frame_equal(out[meta], real[meta])
    assert (out["len_chars"] == out["token"].str.len()).all()
    assert out.attrs["source"] == "template"
    info = out.attrs["template_info"]
    assert 0 < info["n_stems_used"] <= 50 and 0.0 <= info["share_stems_real"] <= 1.0
    assert "prefix" not in out.columns  # goes through decompose_table like the real text


def test_stems_are_grammatical_and_distinct(dfd):
    t = ParadigmTemplate.learn(dfd, n_stems=40)
    pool = t.stem_pool(0)
    assert len(pool) == 40 and len(set(pool)) == 40
    g = learn_word_grammar(dfd["stem"].tolist(), "basic")
    assert all(g.accepts_token(s, "slot") and g.accepts_token(s, "trigram") for s in pool)
    assert t.stem_pool(0, 10) == pool[:10]  # rank order does not depend on n_stems


def test_tokens_are_prefix_stem_suffix(real, dfd):
    t = ParadigmTemplate.learn(dfd, n_stems=30, prefix_dependence="stem_suffix")
    out = t.sample(real, seed=2)
    stems = set(t.stem_pool(2))
    ok = 0
    for tok in out["token"]:
        ok += any(tok.startswith(p) and tok.endswith(s) and tok[len(p): len(tok) - len(s)] in stems
                  for p in ["", *PRE] for s in ["", *SUF] if len(tok) > len(p) + len(s))
    assert ok == len(out)


def test_fixed_paradigm_restricts_suffixes_and_free_does_not(real, dfd):
    fixed = ParadigmTemplate.learn(dfd, n_stems=30, zipf_s=0.5).sample(real, seed=5)
    free = ParadigmTemplate.learn(dfd, n_stems=30, zipf_s=0.5, paradigm="free").sample(real, seed=5)
    mf = paradigm_metrics(decompose_table(fixed, prefixes=PRE, suffixes=SUF)).summary
    mr = paradigm_metrics(decompose_table(free, prefixes=PRE, suffixes=SUF)).summary
    assert mf["mean_suffixes_per_stem"] < mr["mean_suffixes_per_stem"]


def test_prefix_dependence_creates_prefix_suffix_mi(real, dfd):
    def mi(**kw):
        out = ParadigmTemplate.learn(dfd, n_stems=40, **kw).sample(real, seed=7)
        return paradigm_metrics(decompose_table(out, prefixes=PRE, suffixes=SUF)).summary
    none, suff = mi(), mi(prefix_dependence="suffix")
    assert suff["mi_prefix_suffix"] - suff["mi_prefix_suffix_shuffled"] >= none["mi_prefix_suffix"] - none["mi_prefix_suffix_shuffled"]


def test_topic_separates_sections(real, dfd):
    from voytools.morphology import inventory_jaccard

    def jsd(alpha):
        out = ParadigmTemplate.learn(dfd, n_stems=40, topic_concentration=alpha).sample(real, seed=9)
        d = decompose_table(out, prefixes=PRE, suffixes=SUF)
        return float(inventory_jaccard(d, "illustration", column="stem")["jsd"].mean())

    assert jsd(0.05) > jsd(math.inf)


def test_parameter_validation(dfd):
    with pytest.raises(TypeError):
        ParadigmTemplate.learn(dfd, n_stem=10)
    with pytest.raises(ValueError):
        ParadigmTemplate.learn(dfd, paradigm="open")
    with pytest.raises(ValueError):
        ParadigmTemplate.learn(dfd, prefix_dependence="page")
    with pytest.raises(ValueError):
        ParadigmTemplate.learn(dfd, page_locality=1.5)
    with pytest.raises(ValueError):
        ParadigmTemplate.learn(dfd, topic_concentration=0)


def test_free_paradigm_cmi_excess_like_slot_null(zl_tables):
    """paradigm='free' (slot independence of the suffix) gives cmi_excess near the slot null of 002."""
    from voytools.table import filter_tokens

    real = filter_tokens(zl_tables.tokens, generic_type="P", clean=True)
    dfd = decompose_table(real)
    pre, suf = dfd.attrs["prefixes"], dfd.attrs["suffixes"]
    t = ParadigmTemplate.learn(dfd, n_stems=2000, paradigm="free")
    cm_t = paradigm_metrics(decompose_table(t.sample(real, seed=0), prefixes=pre, suffixes=suf), seed=0).summary
    cm_s = paradigm_metrics(decompose_table(slot_independent_sample(dfd, seed=0), prefixes=pre, suffixes=suf),
                            seed=0).summary
    cm_r = paradigm_metrics(dfd, seed=0).summary
    assert abs(cm_t["cmi_excess"] - cm_s["cmi_excess"]) < 0.025
    assert cm_t["cmi_excess"] < cm_r["cmi_excess"] - 0.03
