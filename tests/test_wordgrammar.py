"""Word grammars and the grammar-bound auto-copy generator (experiment 003)."""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from voytools import nulls as N
from voytools.glyphs import collapse_eva
from voytools.table import filter_tokens, load_tables
from voytools.wordgrammar import BEGIN, END, WordGrammar, grammaticality, learn_word_grammar

SEED = ["qokeedy", "chedy", "shedy", "daiin", "okaiin", "qokain", "ol", "or", "ar", "chol", "dar", "s"]
REF = SEED + ["qokedy", "okedy", "cheol", "shol", "dain", "aiin", "qol", "chor", "sheedy", "otedy", "otaiin",
              "dal", "sar", "kar", "okal", "chey", "shey", "qotedy", "cheody", "lchedy"] * 3


def _flat(lines):
    return [t for line in lines for t in line]


def _digest(lines) -> str:
    return hashlib.sha256(json.dumps(lines).encode()).hexdigest()[:16]


# ------------------------------------------------------------------ regression of the default generator


def test_default_autocopy_output_is_unchanged():
    """Fingerprints recorded with the generator before edit_model/new_token_source existed."""
    assert _digest(N.autocopy_generate(SEED, 2000, seed=11)) == "a14bc553fde69ca0"
    variant = dict(window=16, p_copy_exact=0.2, p_two_edits=0.5, recency=0.3, edit_type_probs=(0.4, 0.3, 0.3),
                   p_new=0.05, seed=12)
    assert _digest(N.autocopy_generate(SEED, 2000, **variant)) == "b953e888e0a4f533"
    explicit = N.autocopy_generate(SEED, 2000, edit_model="bigram", new_token_source="seed", info={}, **variant)
    assert _digest(explicit) == "b953e888e0a4f533"


def test_default_autocopy_like_on_real_files_is_unchanged(zl_path, gc_path):
    zl = filter_tokens(load_tables(zl_path, space_policy="split").tokens, generic_type="P", clean=True)
    df = N.autocopy_like(zl, allowed_from=zl["token"].tolist(), seed=5, window=32, p_new=0.15, p_copy_exact=0.5,
                         p_two_edits=0.5)
    assert hashlib.sha256("\n".join(df["token"]).encode()).hexdigest()[:16] == "b31daa2140b4f9bd"
    gc = filter_tokens(load_tables(gc_path, space_policy="split").tokens, generic_type="P", clean=True)
    dg = N.autocopy_like(gc, allowed_from=gc["token"].tolist(), seed=6)
    assert hashlib.sha256("\n".join(dg["token"]).encode()).hexdigest()[:16] == "b54d0db605ca69cb"


# ------------------------------------------------------------------ grammar learning and checks


def test_ngram_models_admit_exactly_the_attested_ngrams():
    g = learn_word_grammar(["abc", "bcd"], preset=None)
    assert g.accepts(("a", "b", "c"), "trigram") and g.accepts(("b", "c", "d"), "trigram")
    assert g.accepts(("a", "b", "c", "d"), "trigram")  # ^ab, abc, bcd, cd$ are all attested
    assert not g.accepts(("c", "d"), "trigram")  # cd$ is attested, ^cd is not
    assert g.accepts(("c", "d"), "bigram") is False  # ^c is not attested either
    assert g.accepts(("b", "c"), "bigram")
    assert not g.accepts(("c", "b"), "bigram")
    assert not g.accepts((), "slot")
    with pytest.raises(ValueError):
        g.accepts(("a",), "fourgram")


def test_trigram_boundaries_use_sentinels_not_glyphs():
    g = learn_word_grammar(["a$b", "$"], preset=None)  # '$' is a glyph in v101
    assert g.accepts(("a", "$", "b"), "trigram") and g.accepts(("$",), "trigram")
    assert (BEGIN, "$", END) in g.trigrams and not g.accepts(("a", "$"), "trigram")


def test_slot_grammar_orders_by_mean_position_and_allows_repeat_runs():
    g = learn_word_grammar(["qok", "qoeek", "ok", "qoeeek", "oek"] * 10, preset=None, target_coverage=1.0)
    assert g.slots.index("q") < g.slots.index("o") < g.slots.index("k")
    assert "e" in g.repeatable and g.max_run["e"] == 3
    assert g.accepts(tuple("qoeek"), "slot") and g.coverage == 1.0
    assert not g.accepts(tuple("koq"), "slot")  # wrong order (and unattested bigrams)
    assert not g.accepts(tuple("qoeeeek"), "slot")  # run longer than attested


def test_slot_grammar_adds_slots_for_glyphs_in_two_positions():
    words = ["oda", "da", "do", "odo", "ad"] * 20
    g = learn_word_grammar(words, preset=None, target_coverage=1.0, min_gain=0.0)
    assert g.coverage == 1.0 and all(g.accepts(tuple(w), "slot") for w in set(words))
    assert len(g.slots) > len(set("".join(words)))


def test_grammaticality_on_real_text_is_close_to_one(zl_path, gc_path):
    for path in (zl_path, gc_path):
        real = filter_tokens(load_tables(path, space_policy="split").tokens, generic_type="P", clean=True)
        g = learn_word_grammar(real["token"], real.attrs.get("glyph_preset"))
        assert grammaticality(real, "trigram", grammar=g) == 1.0
        assert grammaticality(real, "bigram", grammar=g) == 1.0
        assert 0.98 <= grammaticality(real, "slot", grammar=g) <= 1.0
        assert grammaticality(real, "slot", grammar=g, weight="type") >= 0.9


def test_grammaticality_input_forms():
    toks = ["qokeedy", "chedy", "ydeh"]
    g = learn_word_grammar(SEED + ["qokeedy", "chedy"], "basic")
    assert grammaticality(toks, "bigram", grammar=g) == pytest.approx(2 / 3)
    assert grammaticality(toks, "bigram", grammar=g, weight="type") == pytest.approx(2 / 3)
    assert grammaticality(SEED, "trigram") == 1.0  # learned on itself
    assert np.isnan(grammaticality([], "slot", grammar=g))
    with pytest.raises(ValueError):
        grammaticality(toks, weight="page")


# ------------------------------------------------------------------ grammar-bound generator


def test_trigram_edits_never_create_unattested_trigrams():
    g = learn_word_grammar(REF, "basic")
    out = _flat(N.autocopy_generate(SEED, 3000, p_copy_exact=0.0, p_two_edits=0.5, allowed_from=REF,
                                    edit_model="trigram", seed=3))
    assert all(g.accepts(collapse_eva(t, "basic"), "trigram") for t in out)
    bigram = _flat(N.autocopy_generate(SEED, 3000, p_copy_exact=0.0, p_two_edits=0.5, allowed_from=REF, seed=3))
    assert not all(g.accepts(collapse_eva(t, "basic"), "trigram") for t in bigram)  # the constraint binds
    assert len(set(out)) > len(set(SEED))  # and edits still happen


def test_slot_edits_stay_grammatical():
    g = learn_word_grammar(REF, "basic")
    seeds = [t for t in SEED if g.accepts_token(t, "slot")]
    out = _flat(N.autocopy_generate(seeds, 2000, p_copy_exact=0.0, allowed_from=REF, edit_model="slot",
                                    word_grammar=g, seed=4))
    assert all(g.accepts_token(t, "slot") for t in out)


def test_grammar_sampler_returns_only_grammatical_words():
    g = learn_word_grammar(REF, "basic")
    rng = np.random.default_rng(0)
    words = [g.sample(rng) for _ in range(500)]
    assert all(w is not None for w in words)
    assert all(g.accepts(w, "slot") and g.accepts(w, "trigram") for w in words)
    assert len(set(words)) > 5
    assert set(len(w) for w in words) <= set(g.lengths.tolist())


def test_new_tokens_from_the_grammar():
    g = learn_word_grammar(REF, "basic")
    info: dict = {}
    out = _flat(N.autocopy_generate(SEED[:2], 1000, p_new=1.0, p_copy_exact=1.0, allowed_from=REF,
                                    new_token_source="grammar", info=info, seed=5))
    assert all(g.accepts_token(t, "slot") for t in out)
    assert len(set(out) - set(SEED)) > 10
    assert info["n_new"] == 1000 and info["n_grammar_fallback"] == 0


def test_info_counts_unchanged_copies():
    info: dict = {}
    N.autocopy_generate(SEED, 500, p_copy_exact=1.0, info=info, seed=1)
    assert info["p_unchanged"] == 1.0 and info["n_edit_attempted"] == 0
    info = {}
    N.autocopy_generate(SEED, 500, p_copy_exact=0.0, allowed_from=REF, edit_model="trigram", info=info, seed=1)
    assert info["n_edit_attempted"] == 500 and 0 <= info["n_edit_failed"] == info["n_unchanged"] < 500
    assert info["p_unchanged"] == info["n_unchanged"] / 500


def test_autocopy_like_records_info(zl_path):
    zl = filter_tokens(load_tables(zl_path, space_policy="split").tokens, generic_type="P", clean=True).iloc[:3000]
    df = N.autocopy_like(zl, allowed_from=zl["token"].tolist(), edit_model="slot", new_token_source="grammar", seed=2)
    assert len(df) == len(zl)
    assert df.attrs["autocopy_info"]["n_tokens"] == len(zl)
    a = N.autocopy_like(zl, allowed_from=zl["token"].tolist(), edit_model="trigram", seed=2)
    b = N.autocopy_like(zl, allowed_from=zl["token"].tolist(), edit_model="trigram", seed=2)
    assert a["token"].tolist() == b["token"].tolist()


def test_new_parameters_are_validated():
    with pytest.raises(ValueError):
        N.autocopy_generate(SEED, 10, edit_model="fourgram")
    with pytest.raises(ValueError):
        N.autocopy_generate(SEED, 10, new_token_source="dictionary")
    with pytest.raises(ValueError):
        N.autocopy_generate(SEED, 10, edit_model="slot", constrain_bigrams=False)


def test_word_grammar_requires_tokens():
    with pytest.raises(ValueError):
        WordGrammar.fit([], "basic")
