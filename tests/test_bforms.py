"""Tests of voytools.bforms (experiment 019)."""

import numpy as np
import pandas as pd
import pytest

from voytools import bforms as B
from voytools.registers import levenshtein


def _t(s):
    return tuple(s)


def test_neighbours_within_matches_brute_force():
    rng = np.random.default_rng(3)
    alph = list("abcd")
    forms = list({tuple(rng.choice(alph, size=rng.integers(1, 5)).tolist()) for _ in range(80)})
    cands = forms[:50]
    got = B.neighbours_within(forms[50:], cands, 2)
    for q in forms[50:]:
        exp = sorted(((c, levenshtein(q, c)) for c in cands if levenshtein(q, c) <= 2), key=lambda t: (t[1], t[0]))
        assert got[q] == exp


def test_is_d_op_schemes():
    assert B.is_d_op("ins:d@mid")
    assert B.is_d_op("ins:ee@begin")
    assert B.is_d_op("sub:sh>ch@begin")
    assert not B.is_d_op("ins:l@begin")
    assert not B.is_d_op("sub:o>l@begin")
    assert B.is_d_op("ins:8@mid", "v101") and B.is_d_op("sub:1>2@begin", "v101")
    assert not B.is_d_op("ins:o@mid", "v101")
    assert not B.is_d_op("ins:d@mid", None)


def test_classify_forms_classes_and_new_stems():
    cx = {_t("chey"): 3, _t("chedy"): 5, _t("okol"): 2, _t("zzz"): 4, _t("qqq"): 1, _t("dal"): 2}
    cy = {_t("chey"): 4, _t("okal"): 3, _t("dar"): 2}
    stem = {f: "".join(f)[:2] for f in cx}
    stem.update({f: "".join(f)[:2] for f in cy})
    c = B.classify_forms(cx, cy, stem, "extended")
    c.index = ["".join(f) for f in c["form"]]
    assert c.loc["chey", "cls"] == "shared"
    assert c.loc["chedy", "cls"] == "d_variant" and c.loc["chedy", "op"] == "ins:d@mid"
    assert c.loc["okol", "cls"] == "other_le1" and c.loc["okol", "counterpart"] == _t("okal")
    assert c.loc["zzz", "cls"] == "no_counterpart"
    assert c.loc["qqq", "cls"] == "hapax"
    assert c.loc["dal", "cls"] == "other_le1"
    assert bool(c.loc["zzz", "stem_new"]) and not bool(c.loc["chedy", "stem_new"])
    assert not bool(c.loc["dal", "stem_new"])  # stem "da" has tokens in Y (dar)
    comp = B.composition(c)
    assert comp["n_tokens"] == 17
    assert comp["share_only"] == pytest.approx((5 + 2 + 4 + 2) / 17)
    assert comp["share_hapax"] == pytest.approx(1 / 17)
    assert comp["only_d_variant"] == pytest.approx(5 / 13)
    assert comp["only_new_stem"] == pytest.approx(4 / 13)
    assert comp["only_new_stem"] + comp["only_new_form"] == pytest.approx(1)


def test_match_by_count_keeps_counts():
    rows = [{"form": (f"s{i}",), "n_x": c, "cls": "shared"} for i, c in enumerate([2, 2, 3, 5, 5, 9, 1, 1])]
    rows += [{"form": (f"n{i}",), "n_x": c, "cls": "no_counterpart"} for i, c in enumerate([2, 5, 9, 4])]
    df = pd.DataFrame(rows)
    forms, matched = B.match_by_count(df, "no_counterpart", np.random.default_rng(0))
    assert len(forms) == len(matched) == 4
    assert len(set(matched)) == 4
    counts = dict(zip(df["form"], df["n_x"]))
    got = sorted(counts[m] for m in matched)
    assert got in ([2, 3, 5, 9], [2, 5, 5, 9])  # 4 -> nearest 3 or 5
    assert B.match_by_count(df, "hapax", np.random.default_rng(0)) == ([], [])


def test_gries_dp_and_entropy():
    assert B.gries_dp({"a": 5, "b": 5}, {"a": 50, "b": 50}) == pytest.approx(0.0)
    assert B.gries_dp({"a": 10}, {"a": 50, "b": 50}) == pytest.approx(0.5)
    assert B.entropy_bits([1, 1]) == pytest.approx(1.0)
    assert B.entropy_bits([]) == 0.0


def test_bigram_model_scores_seen_forms_higher():
    forms = {_t("chey"): 10, _t("chedy"): 10, _t("okey"): 5}
    vocab = ["c", "h", "e", "y", "d", "o", "k", "z"]
    m = B.bigram_model(forms, vocab)
    assert B.bigram_loglik(_t("chey"), m) > B.bigram_loglik(_t("zzzz"), m)


def _sample():
    toks = ["chey", "chedy", "okal", "chey", "zzz", "zzz", "chedy", "okal"]
    pages = ["p1", "p1", "p1", "p2", "p2", "p2", "p3", "p3"]
    lines = [0, 0, 0, 1, 1, 1, 2, 2]
    form = np.empty(len(toks), dtype=object)
    for i, t in enumerate(toks):
        form[i] = tuple(t)
    nf = np.array([tuple(toks[i + 1])[0] if i + 1 < len(toks) and lines[i + 1] == lines[i] else "" for i in range(len(toks))],
                  dtype=object)
    return {"form": form, "stem": np.array([t[:2] for t in toks], dtype=object),
            "suffix": np.array([t[2:] for t in toks], dtype=object),
            "page": np.array(pages, dtype=object), "section": np.array(["H"] * 6 + ["B"] * 2, dtype=object),
            "line": np.array(lines), "is_first": np.array([l == 0 for l in [0, 1, 2, 0, 1, 2, 0, 1]]),
            "is_last": np.array([False] * 8), "is_par_first": np.array([True] * 3 + [False] * 5),
            "next_first": nf}


def test_form_profile_metrics():
    tok = _sample()
    rng = np.random.default_rng(0)
    shared = {_t("chey"): 2, _t("chedy"): 2, _t("okal"): 2}
    model = B.bigram_model(shared, ["c", "h", "e", "y", "d", "o", "k", "a", "l", "z"])
    ref = np.fromiter((f in shared for f in tok["form"]), dtype=bool, count=8)
    p = B.form_profile([_t("zzz")], tok, {_t("zzz"): (1, 2)}, model, rng, seam_ref=ref)
    assert p["n_forms"] == 1 and p["n_tokens"] == 2
    assert p["paradigm_forms"] == 1
    assert p["pages_per_token"] == pytest.approx(0.5)
    assert p["top_page_share"] == pytest.approx(1.0)
    assert p["section_entropy"] == 0.0
    assert p["label_share"] == pytest.approx(1 / 3) and p["label_any"] == 1.0
    assert p["first_share"] == pytest.approx(0.0)
    assert p["par_first_share"] == pytest.approx(0.0)
    assert p["cooc_ratio"] > 1  # both zzz in one line
    assert np.isfinite(p["bigram_ll"]) and np.isfinite(p["seam_jsd"])
    q = B.form_profile([_t("chey")], tok, None, model, rng)
    assert q["paradigm_forms"] == 2  # chey, chedy share the stem "ch"
    assert q["H_suffix_given_stem"] == pytest.approx(1.0)
    assert B.form_profile([], tok)["n_forms"] == 0


def test_suffix_pairs_and_pair_choice():
    sp = B.suffix_pairs(["chey", "chedy", "dy", "ol", "okey", "ey"], "ey", "edy")
    assert sp["choice"].tolist() == [0, 1, 0] and sp["base"].tolist() == ["ch", "ch", "ok"]
    sp2 = B.suffix_pairs(["chey", "chedy", "oky", "okdy", "dy"], "y", "dy", exclude_base_end="e")
    assert sp2["base"].tolist() == ["ok", "ok"] and sp2["choice"].tolist() == [0, 1]
    ch = np.array([0, 0, 0, 0, 1, 1, 1, 1, 0, 1])
    ctx = np.array(["a"] * 4 + ["b"] * 4 + ["c"] * 2, dtype=object)
    r = B.pair_choice(ch, ctx, min_n=2, n_perm=50, seed=1)
    assert r["H"] == pytest.approx(1.0)
    assert r["H_given"] == pytest.approx(0.2)
    assert r["share_determined"] == pytest.approx(0.8) and r["share_middle"] == pytest.approx(0.2)
    assert r["z"] > 1
    e = B.pair_choice(np.array([]), np.array([]))
    assert e["n"] == 0 and np.isnan(e["H"])


def test_both_forms():
    r = B.both_forms(["a", "a", "b", "b", "c"], [0, 1, 0, 0, 1], ["p1", "p2", "p1", "p2", "p1"])
    assert r["n_bases"] == 2 and r["share_both"] == 0.5
    assert r["share_both_same_page"] == 0.0
    r2 = B.both_forms(["a", "a"], [0, 1], ["p1", "p1"])
    assert r2["share_both_same_page"] == 1.0
