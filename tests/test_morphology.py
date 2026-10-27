"""Morphology: affix inventories, decomposition, paradigms and within-token BPE."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pytest

from voytools import morphology as M
from voytools import stats as S
from tests.conftest import make_tokens

TOY = ["qokeedy", "qokeedy", "qokedy", "okeedy", "daiin", "chedy", "shedy", "chedy"]


# ------------------------------------------------------------------ affix inventory


def test_affix_inventory_prefix_hand_computed():
    inv = M.affix_inventory(make_tokens([TOY]), "prefix", max_len=2, min_count=1).set_index("affix")
    q = inv.loc["q"]
    assert q["count"] == 3 and q["share"] == pytest.approx(3 / 8) and q["n_types"] == 2
    assert q["n_remainders"] == 2
    assert q["n_free"] == 1  # 'okeedy' is itself a token, 'okedy' is not
    assert q["exact"] == 0 and q["ext_share"] == pytest.approx(1.0)  # always continues as 'qo'
    assert inv.loc["qo", "count"] == 3 and inv.loc["qo", "n_glyphs"] == 2
    assert inv.loc["ch", "n_glyphs"] == 1  # 'ch' is a single glyph in the basic preset
    assert inv.loc["ch", "count"] == 2 and inv.loc["ch", "n_types"] == 1
    assert inv.loc["che", "n_glyphs"] == 2
    # sorted by count (descending), ties alphabetically
    assert inv["count"].is_monotonic_decreasing


def test_affix_inventory_suffix_and_exact():
    inv = M.affix_inventory(make_tokens([TOY + ["y", "dy", "d"]]), "suffix", max_len=2, min_count=1).set_index("affix")
    dy = inv.loc["dy"]
    assert dy["count"] == 8 and dy["n_types"] == 6  # 5 words + the bare token 'dy'
    assert dy["exact"] == 1
    assert dy["n_free"] == 0  # 'qokee', 'qoke', 'okee', 'che', 'she' are no tokens
    assert inv.loc["y", "exact"] == 1
    assert inv.loc["y", "ext_share"] == pytest.approx(8 / 9)
    assert inv.loc["y", "n_free"] == 1  # only the remainder 'd' (from 'dy') is a token itself


def test_affix_inventory_min_count_by_and_errors():
    df = make_tokens([TOY]).assign(g=["a"] * 4 + ["b"] * 4)
    inv = M.affix_inventory(df, "prefix", max_len=1, min_count=2, by="g")
    assert set(inv["g"]) <= {"a", "b"} and (inv["count"] >= 2).all()
    assert inv.columns[0] == "g"
    assert inv.loc[inv["g"] == "a", "affix"].tolist()[0] == "q"
    with pytest.raises(ValueError):
        M.affix_inventory(df, "middle")
    empty = M.affix_inventory(make_tokens([TOY]), "prefix", min_count=1000)
    assert len(empty) == 0 and "affix" in empty.columns


def test_affix_inventory_char_unit_treats_ch_as_two_chars():
    inv = M.affix_inventory(make_tokens([TOY]), "prefix", max_len=1, min_count=1, unit="char").set_index("affix")
    assert inv.loc["c", "count"] == 2 and "ch" not in inv.index
    with pytest.raises(ValueError):
        M.affix_inventory(make_tokens([TOY]), "prefix", unit="syllable")


def test_learn_affixes_skips_dominated_affixes():
    df = make_tokens([TOY * 3])
    pre = M.learn_affixes(df, "prefix", top_n=10, max_len=2, min_count=2)
    assert "qo" in pre and "q" not in pre  # q is (almost) always qo
    lens = [len(M.collapse_eva(a, "basic")) for a in pre]
    assert lens == sorted(lens, reverse=True)  # longest first, ready for greedy matching
    assert M.learn_affixes(df, "suffix", top_n=3, max_len=2, min_count=10_000) == []


# ------------------------------------------------------------------ decomposition


def test_decompose_known_splits():
    assert M.decompose("qokeedy", ["qo", "o"], ["edy", "dy", "y"]) == ("qo", "ke", "edy")
    assert M.decompose("chedy", ["ch", "c"], ["edy"]) == ("ch", "edy", "")  # stem must keep >= 1 glyph
    assert M.decompose("daiin", ["qo"], ["aiin", "iin"]) == ("", "d", "aiin")
    assert M.decompose("or", ["o"], ["or", "r"]) == ("o", "r", "")
    assert M.decompose("x", ["x"], ["x"]) == ("", "x", "")
    assert M.decompose("qo", ["qo"], []) == ("", "qo", "")  # nothing left for the stem
    assert M.decompose("qo", ["qo"], [], min_stem=0) == ("qo", "", "")


def test_decompose_longest_prefix_and_glyph_awareness():
    # 'ch' is one glyph under 'basic': the prefix 'c' can never match inside it
    assert M.decompose("chol", ["c", "ch"], ["l"], preset="basic") == ("ch", "o", "l")
    assert M.decompose("chol", ["c"], ["l"], preset=None) == ("c", "ho", "l")
    # glyph sequences are accepted too
    assert M.decompose(("q", "o", "k", "y"), ["qo"], ["y"]) == ("qo", "k", "y")


def test_decompose_table_reassembles_tokens_and_stores_lists():
    df = make_tokens([TOY * 2, TOY])
    out = M.decompose_table(df, prefixes=["qo", "o", "ch", "sh"], suffixes=["edy", "dy", "aiin"])
    assert len(out) == len(df) and "prefix" in out.columns
    assert ((out["prefix"] + out["stem"] + out["suffix"]) == out["token"]).all()
    assert (out["stem"].str.len() >= 1).all()
    assert out.attrs["prefixes"] == ["qo", "o", "ch", "sh"] and out.attrs["suffixes"] == ["edy", "dy", "aiin"]
    assert out.attrs.get("alphabet") == df.attrs.get("alphabet")
    assert "prefix" not in df.columns  # input untouched
    row = out[out["token"] == "qokeedy"].iloc[0]
    assert (row["prefix"], row["stem"], row["suffix"]) == ("qo", "ke", "edy")
    row = out[out["token"] == "daiin"].iloc[0]
    assert (row["prefix"], row["stem"], row["suffix"]) == ("", "d", "aiin")


def test_decompose_suffix_first_order():
    # default unchanged; mirrored order takes the longest suffix first
    assert M.decompose("olor", ["ol"], ["lor", "r"]) == ("ol", "o", "r")
    assert M.decompose("olor", ["ol"], ["lor", "r"], order="suffix_first") == ("", "o", "lor")
    assert M.decompose("qokeedy", ["qo", "o"], ["edy", "dy", "y"], order="suffix_first") == ("qo", "ke", "edy")
    assert M.decompose("dy", ["d"], ["dy", "y"], order="suffix_first") == ("", "d", "y")  # stem keeps >= 1 glyph
    with pytest.raises(ValueError):
        M.decompose("dy", ["d"], ["y"], order="middle")
    df = make_tokens([TOY * 2, ["olor", "olor", "chedy"]])
    a = M.decompose_table(df, prefixes=["qo", "o", "ol"], suffixes=["edy", "dy", "lor", "r"])
    b = M.decompose_table(df, prefixes=["qo", "o", "ol"], suffixes=["edy", "dy", "lor", "r"], order="prefix_first")
    assert a[["prefix", "stem", "suffix"]].equals(b[["prefix", "stem", "suffix"]])
    c = M.decompose_table(df, prefixes=["qo", "o", "ol"], suffixes=["edy", "dy", "lor", "r"], order="suffix_first")
    assert ((c["prefix"] + c["stem"] + c["suffix"]) == c["token"]).all()
    row = c[c["token"] == "olor"].iloc[0]
    assert (row["prefix"], row["stem"], row["suffix"]) == ("", "o", "lor")


def test_decompose_table_learns_and_reference_lists():
    df = make_tokens([TOY * 20])
    learned = M.decompose_table(df, min_count=5, max_len=2, n_prefixes=3, n_suffixes=3)
    assert learned.attrs["prefixes"] and learned.attrs["suffixes"]
    assert ((learned["prefix"] + learned["stem"] + learned["suffix"]) == learned["token"]).all()
    ref = M.decompose_table(df, prefixes="reference", suffixes="reference")
    assert ref.attrs["prefixes"] == list(M.REFERENCE_PREFIXES)
    assert ref.attrs["suffixes"] == list(M.REFERENCE_SUFFIXES)
    assert ((ref["prefix"] + ref["stem"] + ref["suffix"]) == ref["token"]).all()
    # learning on another table
    other = M.decompose_table(make_tokens([["qokedy", "qoky"]]), learn_on=df, min_count=5, max_len=2, n_prefixes=3)
    assert other.attrs["prefixes"] == learned.attrs["prefixes"] == M.learn_affixes(df, "prefix", 3, 2, 5)
    with pytest.raises(ValueError):
        M.decompose_table(df, prefixes="nonsense")


# ------------------------------------------------------------------ paradigms


def _decomposed(tokens, prefixes, suffixes):
    return M.decompose_table(make_tokens([tokens]), prefixes=prefixes, suffixes=suffixes)


def test_paradigm_table_counts():
    dfd = _decomposed(["qokedy", "kedy", "qoky", "ky", "qoky"], ["qo"], ["dy", "y"])
    t = M.paradigm_table(dfd)
    assert t.loc["ke", "dy"] == 2 and t.loc["k", "y"] == 3
    assert t.index.name == "stem" and t.loc["k", "dy"] == 0 and t.loc["ke", "y"] == 0
    assert t.sum().sum() == 5
    ty = M.paradigm_table(dfd, counts="types")
    assert ty.loc["k", "y"] == 2 and ty.loc["ke", "dy"] == 2  # qoky/ky and qokedy/kedy
    assert M.paradigm_table(dfd, min_stem_count=3).index.tolist() == ["k"]
    pre = M.paradigm_table(dfd, axis="prefix")
    assert M.NULL_LABEL in pre.columns and pre.loc["ke", "qo"] == 1 and pre.loc["k", "qo"] == 2
    assert M.paradigm_table(dfd, top_stems=1).shape[0] == 1
    assert M.paradigm_table(dfd, top_suffixes=1).shape[1] == 1
    with pytest.raises(ValueError):
        M.paradigm_table(dfd, axis="stem")
    with pytest.raises(ValueError):
        M.paradigm_table(dfd, counts="chars")


def test_paradigm_metrics_detects_prefix_suffix_dependence():
    # within the single stem 'k' the prefix determines the suffix: I(P;S|stem) = 1 bit
    dfd = _decomposed(["qoky"] * 10 + ["okdy"] * 10, ["qo", "o"], ["dy", "y"])
    pm = M.paradigm_metrics(dfd, min_stem_count=5, shuffles=30, seed=1)
    s = pm.summary
    assert s["n_tokens"] == 20 and s["n_stems"] == 1 and s["n_stems_k"] == 1
    assert s["cmi_prefix_suffix_given_stem"] == pytest.approx(1.0)
    assert s["mi_prefix_suffix"] == pytest.approx(1.0)
    assert s["cmi_excess"] > 0.8 and s["cmi_shuffled"] < 0.2
    assert s["H_suffix_given_stem"] == pytest.approx(1.0)
    assert s["mean_suffixes_per_stem"] == 2 and s["frac_stems_multi_suffix"] == 1.0
    assert pm.per_stem.loc[0, "count"] == 20 and pm.per_stem.loc[0, "n_types"] == 2


def test_paradigm_metrics_independent_affixes_have_no_excess():
    rng = np.random.default_rng(0)
    pre = rng.choice(["qo", "o"], size=600)
    suf = rng.choice(["dy", "y"], size=600)
    toks = [p + "k" + s for p, s in zip(pre, suf)]
    dfd = _decomposed(toks, ["qo", "o"], ["dy", "y"])
    s = M.paradigm_metrics(dfd, min_stem_count=5, shuffles=30, seed=2).summary
    assert abs(s["cmi_excess"]) < 0.02 and s["cmi_prefix_suffix_given_stem"] < 0.03


def test_paradigm_metrics_ignores_rare_stems_and_is_reproducible():
    dfd = _decomposed(["qoky"] * 6 + ["odz"], ["qo", "o"], ["y", "z"])
    a = M.paradigm_metrics(dfd, min_stem_count=5, shuffles=5, seed=3)
    b = M.paradigm_metrics(dfd, min_stem_count=5, shuffles=5, seed=3)
    assert a.summary == b.summary
    assert a.summary["n_stems"] == 2 and a.summary["n_stems_k"] == 1
    assert a.summary["token_coverage_k"] == pytest.approx(6 / 7)
    empty = M.paradigm_metrics(dfd, min_stem_count=100, shuffles=2)
    assert np.isnan(empty.summary["cmi_prefix_suffix_given_stem"]) and len(empty.per_stem) == 0


def test_paradigm_metrics_is_invariant_under_token_reordering():
    from voytools.nulls import shuffle_within_lines

    df = make_tokens([TOY * 5 for _ in range(6)])
    dfd = M.decompose_table(df, prefixes=["qo", "o"], suffixes=["edy", "dy"])
    sh = M.decompose_table(shuffle_within_lines(df, seed=4), prefixes=["qo", "o"], suffixes=["edy", "dy"])
    a = M.paradigm_metrics(dfd, min_stem_count=2, shuffles=0).summary
    b = M.paradigm_metrics(sh, min_stem_count=2, shuffles=0).summary
    for k in ("n_stems", "n_stems_k", "H_suffix_given_stem", "mi_prefix_suffix", "cmi_prefix_suffix_given_stem"):
        assert a[k] == pytest.approx(b[k])


# ------------------------------------------------------------------ inventories, productivity, dispersion


def _with_groups(tokens, groups, prefixes, suffixes):
    dfd = _decomposed(tokens, prefixes, suffixes)
    dfd["illustration"] = groups
    return dfd


def test_inventory_jaccard_hand_computed():
    toks = ["kedy"] * 3 + ["ky"] * 3 + ["tedy"] * 3 + ["kedy"] * 3 + ["py"] * 2
    grp = ["H"] * 9 + ["B"] * 5
    dfd = _with_groups(toks, grp, [], ["edy", "y"])
    st = M.inventory_jaccard(dfd, "illustration", column="stem", min_count=3)
    assert st.shape[0] == 1 and (st.loc[0, "group_a"], st.loc[0, "group_b"]) == ("B", "H")
    # H stems {k, t}; B stems with >= 3: {k} ('p' occurs twice)
    assert st.loc[0, "n_inv_a"] == 1 and st.loc[0, "n_inv_b"] == 2 and st.loc[0, "jaccard"] == pytest.approx(0.5)
    su = M.inventory_jaccard(dfd, "illustration", column="suffix", min_count=2)
    assert su.loc[0, "jaccard"] == pytest.approx(1.0)  # {edy, y} in both
    assert 0 < su.loc[0, "jsd"] < 1
    same = M.inventory_jaccard(_with_groups(["ky"] * 6, ["H"] * 3 + ["B"] * 3, [], ["y"]), "illustration", "stem", 3)
    assert same.loc[0, "jaccard"] == 1.0 and same.loc[0, "jsd"] == pytest.approx(0.0)
    empty = M.inventory_jaccard(dfd, "illustration", column="stem", min_count=50)
    assert np.isnan(empty.loc[0, "jaccard"])


def test_inventory_jaccard_drops_missing_groups_and_respects_order():
    dfd = _with_groups(["ky", "ky", "ty"], ["H", None, "B"], [], ["y"])
    out = M.inventory_jaccard(dfd, "illustration", "stem", min_count=1, groups=["H", "B", "Z"])
    assert out[["group_a", "group_b"]].values.tolist() == [["H", "B"]]
    assert out.loc[0, "n_tokens_a"] == 1 and out.loc[0, "jaccard"] == 0.0


def test_productivity_slope():
    # stem with frequency f takes exactly f distinct suffixes -> slope 1, intercept 0
    sufs = ["a", "b", "c", "d", "e", "f", "g", "h"]
    toks = ["k" + s for s in sufs[:5]] + ["t" + s for s in sufs[:8]]
    dfd = M.decompose_table(make_tokens([toks]), prefixes=[], suffixes=sufs)
    r = M.productivity_slope(dfd, min_stem_count=5)
    assert r["n_stems"] == 2 and r["slope"] == pytest.approx(1.0) and r["intercept"] == pytest.approx(0.0, abs=1e-9)
    # saturated paradigm: both stems use the same 2 suffixes -> slope 0
    toks = ["ka", "kb"] * 3 + ["ta", "tb"] * 10
    r = M.productivity_slope(M.decompose_table(make_tokens([toks]), prefixes=[], suffixes=["a", "b"]), 5)
    assert r["slope"] == pytest.approx(0.0, abs=1e-12)
    assert np.isnan(M.productivity_slope(dfd, min_stem_count=100)["slope"])
    with pytest.raises(ValueError):
        M.productivity_slope(dfd, axis="stem")


def test_stem_dispersion_even_random_and_clustered():
    # 'k' every 4th token in a 40-token text: ratio exactly 1
    toks = (["kdy"] + ["tdy", "pdy", "fdy"]) * 10
    dfd = M.decompose_table(make_tokens([toks]), prefixes=[], suffixes=["dy"])
    per = M.stem_dispersion(dfd, min_stem_count=5, per_stem=True)
    assert set(per["stem"]) == {"k", "t", "p", "f"} and (per["ratio"] == 1.0).all()
    assert M.stem_dispersion(dfd, min_stem_count=5) == 1.0
    # all occurrences in one burst at the start -> ratio 1/(N/k) = 10/40
    toks = ["kdy"] * 10 + ["tdy", "pdy", "fdy"] * 10
    dfd = M.decompose_table(make_tokens([toks]), prefixes=[], suffixes=["dy"])
    per = M.stem_dispersion(dfd, per_stem=True).set_index("stem")
    assert per.loc["k", "ratio"] == pytest.approx(1 / (40 / 10))
    # independent random positions: about ln 2
    rng = np.random.default_rng(0)
    toks = [f"{c}dy" for c in rng.choice(list("ktpf"), size=4000)]
    dfd = M.decompose_table(make_tokens([toks]), prefixes=[], suffixes=["dy"])
    assert abs(M.stem_dispersion(dfd) - np.log(2)) < 0.1
    assert np.isnan(M.stem_dispersion(dfd, min_stem_count=10**6))


def test_slot_independent_sample_keeps_marginals_and_layout():
    from voytools.nulls import shuffle_within_lines  # noqa: F401  (same table contract)

    rng = np.random.default_rng(1)
    pre = rng.choice(["qo", "o"], size=300)
    suf = np.where(pre == "qo", "dy", "y")  # suffix determined by prefix
    toks = [p + "k" + s for p, s in zip(pre, suf)]
    df = make_tokens([toks[i : i + 10] for i in range(0, 300, 10)])
    dfd = M.decompose_table(df, prefixes=["qo", "o"], suffixes=["dy", "y"])
    a = M.slot_independent_sample(dfd, seed=3)
    b = M.slot_independent_sample(dfd, seed=3)
    assert a["token"].tolist() == b["token"].tolist()
    assert len(a) == len(dfd) and a["locus_id"].tolist() == dfd["locus_id"].tolist()
    assert "stem" not in a.columns and a.attrs["source"] == "slot_independent"
    assert (a["len_glyphs"] == a["token"].map(len)).all()
    re = M.decompose_table(a, prefixes=["qo", "o"], suffixes=["dy", "y"])
    assert abs((re["prefix"] == "qo").mean() - (dfd["prefix"] == "qo").mean()) < 0.12
    real = M.paradigm_metrics(dfd, min_stem_count=5, shuffles=20, seed=0).summary
    null = M.paradigm_metrics(re, min_stem_count=5, shuffles=20, seed=0).summary
    assert real["cmi_excess"] > 0.5 and abs(null["cmi_excess"]) < 0.05


# ------------------------------------------------------------------ BPE


CLASSIC = ["low"] * 5 + ["lower"] * 2 + ["newest"] * 6 + ["widest"] * 3


def test_bpe_classic_example():
    merges = M.learn_bpe(CLASSIC, merges=4, preset=None)
    assert merges == [("s", "t"), ("e", "st"), ("o", "w"), ("l", "ow")]
    assert M.apply_bpe(["lowest", "newest", "low"], merges) == [("low", "est"), ("n", "e", "w", "est"), ("low",)]
    assert M.apply_bpe(["lowest"], merges, n_merges=2) == [("l", "o", "w", "est")]
    assert M.apply_bpe(["lowest"], merges, n_merges=0) == [tuple("lowest")]


def test_bpe_never_merges_across_tokens_and_respects_min_pair_count():
    assert M.learn_bpe(["a", "a", "a", "a"], merges=5) == []
    assert M.learn_bpe(["ab", "cd"], merges=5, min_pair_count=2) == []
    assert M.learn_bpe(["ab", "cd"], merges=5, min_pair_count=1) == [("c", "d"), ("a", "b")]
    # nested input (lines of tokens) gives the same merges as the flat list
    assert M.learn_bpe([CLASSIC[:7], CLASSIC[7:]], merges=4) == M.learn_bpe(CLASSIC, merges=4)


def test_bpe_uses_eva_glyphs_as_initial_units():
    merges = M.learn_bpe(["chedy"] * 5 + ["shedy"] * 5, merges=1, preset="basic")
    assert merges == [("e", "d")]  # ties (e,d)=(d,y)=10 go to the lexicographically larger pair
    units = M.apply_bpe(["chedy"], [], preset="basic")[0]
    assert units == ("ch", "e", "d", "y")


def test_is_nested_tokens():
    assert M.is_nested_tokens([["a", "b"], ["c"]]) and not M.is_nested_tokens(["a", "b"])
    assert not M.is_nested_tokens([("a", "b"), ("c",)]) and not M.is_nested_tokens([])
    assert M.is_nested_tokens([[("a", "b")], [("c",)]])


def test_bpe_dependence_curve_columns_and_monotonicity():
    rng = np.random.default_rng(3)
    words = ["qokeedy", "chedy", "shedy", "daiin", "okaiin", "qokain", "ol", "or", "ar", "chol"]
    lines = [rng.choice(words, size=8).tolist() for _ in range(80)]
    curve = M.bpe_dependence_curve(lines, (0, 4, 16, 64), preset="basic")
    assert list(curve.columns) == [
        "merges", "n_units", "n_symbols", "mean_glyphs_per_unit", "H1", "H2", "gap", "bits_per_glyph",
    ]
    assert curve["merges"].is_monotonic_increasing
    assert curve["n_units"].is_monotonic_decreasing and curve["mean_glyphs_per_unit"].is_monotonic_increasing
    assert (curve["gap"] == curve["H1"] - curve["H2"]).all()
    # at 0 merges the units are glyphs and the numbers equal the direct computation
    glyph_lines = [[g for t in l for g in M.collapse_eva(t, "basic")] for l in lines]
    flat = [g for l in glyph_lines for g in l]
    row = curve.iloc[0]
    assert row["merges"] == 0 and row["n_units"] == len(flat)
    assert row["H1"] == pytest.approx(S.entropy(Counter(flat)))
    assert row["H2"] == pytest.approx(S.conditional_entropy(glyph_lines, 1))
    assert row["gap"] == pytest.approx(row["H1"] - row["H2"])
    assert row["mean_glyphs_per_unit"] == pytest.approx(1.0)


def test_bpe_curve_flat_input_and_separate_training():
    toks = CLASSIC * 3
    a = M.bpe_dependence_curve(toks, (0, 2, 4))
    assert len(a) == 3 and a["merges"].tolist() == [0, 2, 4]
    b = M.bpe_dependence_curve(toks, (0, 2, 4), train=CLASSIC)
    assert b["n_units"].tolist() == a["n_units"].tolist()
    # asking for more merges than exist is clipped, not an error
    c = M.bpe_dependence_curve(["ab", "ab", "ab"], (0, 50))
    assert c["merges"].tolist() == [0, 1]
    with pytest.raises(ValueError):
        M.bpe_dependence_curve(CLASSIC, ())
