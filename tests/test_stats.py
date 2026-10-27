"""Statistics: known values on toy data plus behavioural checks."""

from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from voytools import stats as S
from tests.conftest import make_tokens


# ------------------------------------------------------------------ entropy


def test_entropy_known_values():
    assert S.entropy([1, 1, 1, 1]) == pytest.approx(2.0)
    assert S.entropy({"a": 3, "b": 1}) == pytest.approx(0.8112781244591328)
    assert S.entropy(Counter("aabb")) == pytest.approx(1.0)
    assert S.entropy([5]) == 0.0 and S.entropy([]) == 0.0 and S.entropy([0, 0, 4]) == 0.0
    assert S.entropy(pd.Series([2, 2])) == pytest.approx(1.0)
    assert S.perplexity(2.0) == pytest.approx(4.0) and S.perplexity(0.0) == 1.0


def test_conditional_entropy_known_values():
    assert S.conditional_entropy("aaab") == pytest.approx(0.9182958340544896)
    assert S.conditional_entropy("abababab") == 0.0
    assert S.conditional_entropy([]) == 0.0 and S.conditional_entropy("a") == 0.0
    # flat sequence crosses the boundary, nested does not
    assert S.conditional_entropy("aabb") == pytest.approx(2 / 3, abs=1e-12)
    assert S.conditional_entropy([list("aa"), list("bb")]) == 0.0
    # tokens (strings) as symbols
    assert S.conditional_entropy(["x", "y", "x", "y", "x"]) == 0.0
    with pytest.raises(ValueError):
        S.conditional_entropy("abc", order=0)


def _reference_cond_entropy(seq, order):
    joint = Counter(tuple(seq[i : i + order + 1]) for i in range(len(seq) - order))
    ctx = Counter()
    for g, c in joint.items():
        ctx[g[:-1]] += c
    n = sum(joint.values())
    return -sum(c / n * math.log2(c / ctx[g[:-1]]) for g, c in joint.items())


@pytest.mark.parametrize("order", [1, 2, 3])
def test_conditional_entropy_matches_reference_formula(order):
    rng = np.random.default_rng(order)
    seq = rng.choice(list("abcde"), size=400, p=[0.4, 0.3, 0.15, 0.1, 0.05]).tolist()
    assert S.conditional_entropy(seq, order) == pytest.approx(_reference_cond_entropy(seq, order), abs=1e-9)


def test_glyph_entropy_report_toy():
    df = make_tokens([["ch", "ch"], ["ch", "ch"]])
    basic = S.glyph_entropy_report(df, "basic").iloc[0]
    assert basic["n_glyphs"] == 4 and basic["n_symbols"] == 1 and basic["H1"] == 0.0 and basic["H2"] == 0.0
    chars = S.glyph_entropy_report(df, None).iloc[0]
    assert chars["n_symbols"] == 2 and chars["H1"] == pytest.approx(1.0) and chars["H2"] == pytest.approx(0.0)
    assert chars["perplexity"] == pytest.approx(1.0) and chars["mutual_info"] == pytest.approx(1.0)
    with_b = S.glyph_entropy_report(df, "basic", boundary="_").iloc[0]
    assert with_b["n_symbols"] == 2 and with_b["n_glyphs"] == 6
    assert with_b["H1"] == pytest.approx(S.entropy([4, 2]))


def test_glyph_entropy_report_by_group_and_total():
    df = make_tokens([["ab", "ab"], ["cd", "cd"]], page_size=1)
    df = df.assign(g=["x", "x", "y", "y"])
    rep = S.glyph_entropy_report(df, None, by="g")
    assert rep["g"].tolist() == ["x", "y", "all"]
    assert rep["n_glyphs"].tolist() == [4, 4, 8]
    assert rep.iloc[:2]["H1"].tolist() == [1.0, 1.0] and rep.iloc[2]["H1"] == pytest.approx(2.0)


# ------------------------------------------------------------------ divergence, distance


def test_jsd_properties_and_values():
    assert S.jsd([1, 2, 3], [1, 2, 3]) == 0.0
    assert S.jsd([1, 0], [0, 1]) == pytest.approx(1.0)
    assert S.jsd([1, 0], [0.5, 0.5]) == pytest.approx(0.31127812445913283)
    p, q = [3, 1, 0, 2], [1, 1, 1, 1]
    assert S.jsd(p, q) == pytest.approx(S.jsd(q, p))
    assert S.jsd([10, 20], [1, 2]) == pytest.approx(0.0, abs=1e-12)  # normalised
    assert math.isnan(S.jsd([0, 0], [1, 1]))
    with pytest.raises(ValueError):
        S.jsd([1, 2], [1, 2, 3])


def test_levenshtein_known_values():
    assert S.levenshtein("kitten", "sitting") == 3
    assert S.levenshtein("", "abc") == 3 and S.levenshtein("abc", "") == 3 and S.levenshtein("", "") == 0
    assert S.levenshtein("daiin", "daiin") == 0
    assert S.levenshtein("chedy", "shedy") == 1
    assert S.levenshtein("ab", "ba") == 2
    assert S.levenshtein("flaw", "lawn") == S.levenshtein("lawn", "flaw") == 2


def test_within_edit_distance_agrees_with_levenshtein():
    rng = np.random.default_rng(1)
    words = ["".join(rng.choice(list("abc"), size=int(rng.integers(0, 6)))) for _ in range(150)]
    for a in words[:40]:
        for b in words:
            d = S.levenshtein(a, b)
            for k in (0, 1, 2):
                assert S.within_edit_distance(a, b, k) == (d <= k), (a, b, k, d)


# ------------------------------------------------------------------ frequency / vocabulary


def test_token_frequency_ranks_and_ties():
    df = make_tokens([["b", "a", "a", "c"], ["b", "d"]])
    f = S.token_frequency(df)
    assert f["token"].tolist() == ["a", "b", "c", "d"]  # ties alphabetical
    assert f["count"].tolist() == [2, 2, 1, 1] and f["rank"].tolist() == [1, 2, 3, 4]
    assert f["freq"].sum() == pytest.approx(1.0)
    g = S.token_frequency(df.assign(grp=["x"] * 4 + ["y"] * 2), by="grp")
    assert g["grp"].tolist() == ["x", "x", "x", "y", "y"]
    assert g[g["grp"] == "y"]["rank"].tolist() == [1, 2]


def test_vocabulary_summary():
    df = make_tokens([["a", "a", "a", "b", "cc"]])
    v = S.vocabulary_summary(df).iloc[0]
    assert v["n_tokens"] == 5 and v["n_types"] == 3 and v["n_hapax"] == 2
    assert v["hapax_share"] == pytest.approx(2 / 3) and v["ttr"] == pytest.approx(0.6)
    assert v["top50_coverage"] == pytest.approx(1.0) and v["mean_len"] == pytest.approx(1.2)
    by = S.vocabulary_summary(df.assign(g=["u", "u", "v", "v", "v"]), by="g")
    assert by["n_tokens"].tolist() == [2, 3] and by["n_types"].tolist() == [1, 3]


def test_vocabulary_summary_handles_missing_group_values():
    df = make_tokens([["a", "b"]]).assign(g=[None, "x"])
    by = S.vocabulary_summary(df, by="g")
    assert by["g"].tolist() == ["x", "NA"]


def test_zipf_slope_on_perfect_zipf():
    counts = [1000.0 / r for r in range(1, 201)]
    assert S.zipf_slope(counts) == pytest.approx(-1.0, abs=1e-9)
    fit = S.zipf_fit(counts)
    assert fit["r2"] == pytest.approx(1.0) and fit["n_ranks"] == 200
    assert S.zipf_slope([1000.0 / r**2 for r in range(1, 51)]) == pytest.approx(-2.0, abs=1e-9)
    df = make_tokens([["a"] * 8 + ["b"] * 4 + ["c"] * 2 + ["d"]])
    expected = np.polyfit(np.log([1, 2, 3, 4]), np.log([8, 4, 2, 1]), 1)[0]
    assert S.zipf_slope(S.token_frequency(df)) == pytest.approx(expected, abs=1e-9)
    assert math.isnan(S.zipf_slope([5]))


def test_length_distribution():
    df = make_tokens([["a", "bb", "cc", "chd"]])
    chars = S.length_distribution(df, unit="chars")
    assert chars["length"].tolist() == [1, 2, 3] and chars["count"].tolist() == [1, 2, 1]
    glyphs = S.length_distribution(df, unit="glyphs")  # 'chd' = ch + d
    assert glyphs["length"].tolist() == [1, 2] and glyphs["count"].tolist() == [1, 3]
    with pytest.raises(ValueError):
        S.length_distribution(df, unit="bytes")


def test_ngram_counts():
    c = S.ngram_counts("abab", 2)
    assert dict(zip(c["ngram"], c["count"])) == {("a", "b"): 2, ("b", "a"): 1}
    nested = S.ngram_counts([["a", "b"], ["a", "b"]], 2)
    assert dict(zip(nested["ngram"], nested["count"])) == {("a", "b"): 2}
    toks = S.ngram_counts(["daiin", "ol", "daiin", "ol"], 1)
    assert toks.iloc[0]["count"] == 2 and toks["share"].sum() == pytest.approx(1.0)
    assert len(S.ngram_counts("ab", 3)) == 0


# ------------------------------------------------------------------ mutual information


def test_adjacent_token_mi_detects_order_and_ignores_independence():
    cycle = [["x", "y", "z"] * 4 for _ in range(40)]
    dep = S.adjacent_token_mi(make_tokens(cycle), shuffles=20).iloc[0]
    assert dep["mi_obs"] == pytest.approx(math.log2(3), abs=0.05)
    assert dep["excess_bits"] > 1.2 and dep["excess_share"] > 0.7 and dep["n_pairs"] == 40 * 11
    rng = np.random.default_rng(0)
    indep = [rng.choice(list("pqrst"), size=12).tolist() for _ in range(80)]
    ind = S.adjacent_token_mi(make_tokens(indep), shuffles=50, seed=3).iloc[0]
    assert abs(ind["excess_bits"]) < 0.03
    again = S.adjacent_token_mi(make_tokens(indep), shuffles=50, seed=3).iloc[0]
    assert again["excess_bits"] == ind["excess_bits"]  # deterministic


def test_adjacent_token_mi_cap_and_groups():
    rng = np.random.default_rng(2)
    lines = [[f"t{rng.integers(0, 50)}" for _ in range(8)] for _ in range(40)]
    df = make_tokens(lines)
    r = S.adjacent_token_mi(df, cap=10, shuffles=5).iloc[0]
    assert r["n_types"] == 11  # 10 kept + '<other>'
    by = S.adjacent_token_mi(df.assign(g=["a"] * (len(df) // 2) + ["b"] * (len(df) - len(df) // 2)), shuffles=5, by="g")
    assert by["g"].tolist() == ["a", "b"]


def test_edge_glyph_mi():
    lines = [["ab", "ca"] * 6 for _ in range(30)]
    r = S.edge_glyph_mi(make_tokens(lines), shuffles=30).iloc[0]
    assert r["mi_obs"] == pytest.approx(1.0, abs=0.02) and r["excess_bits"] > 0.8
    rng = np.random.default_rng(5)
    rnd = [["".join(rng.choice(list("abc"), size=3)) for _ in range(10)] for _ in range(60)]
    r2 = S.edge_glyph_mi(make_tokens(rnd), shuffles=40).iloc[0]
    assert abs(r2["excess_bits"]) < 0.03


# ------------------------------------------------------------------ positional profile


def test_positional_profile_classes_and_shares():
    text = (
        "#=IVTFF Eva- 2.0 M\n<f1r>\n"
        "<f1r.1,@P0>   <%>pa.oa.oa.oa\n"
        "<f1r.2,+P0>   ta.oa.oa.sa<$>\n"
        "<f1r.3,*P0>   <%>x<$>\n"
    )
    from voytools.ivtff import parse_ivtff_text
    from voytools.table import build_tables

    df = build_tables(parse_ivtff_text(text)).tokens
    prof = S.positional_profile(df, which="first")
    assert prof.n.to_dict() == {
        "line_first": 1, "mid": 2, "line_last": 1,
        "parfirst_line_first": 1, "parfirst_mid": 2, "parfirst_line_last": 1,
    }  # the one-token line is ignored
    assert prof.shares.loc["t", "line_first"] == 1.0 and prof.shares.loc["p", "parfirst_line_first"] == 1.0
    assert prof.shares.loc["o", "mid"] == 1.0 and prof.shares.loc["s", "line_last"] == 1.0
    assert prof.jsd.loc["line_first", "line_first"] == 0.0
    assert prof.jsd.loc["line_first", "mid"] == pytest.approx(1.0)
    assert prof.jsd.loc["mid", "parfirst_mid"] == pytest.approx(0.0)
    assert prof.top("mid", 1).index.tolist() == ["o"]
    last = S.positional_profile(df, which="last")
    assert last.shares.loc["a"].dropna().eq(1.0).all() and last.jsd.fillna(0).to_numpy().max() == 0.0
    with pytest.raises(ValueError):
        S.positional_profile(df, which="middle")


def test_positional_profile_permutation_test():
    rng = np.random.default_rng(7)
    lines = []
    for _ in range(200):
        n = int(rng.integers(4, 9))
        lines.append(["p" + "a" * 2] + [rng.choice(["oa", "ob", "oc"]) for _ in range(n - 1)])
    df = make_tokens(lines)
    prof = S.positional_profile(df, which="first", permutations=50, seed=1)
    assert prof.pvalues is not None
    assert prof.pvalues.loc["line_first", "mid"] <= 0.05  # clearly different distributions
    # a genuinely null comparison: every token drawn from the same distribution
    null_lines = [[rng.choice(["oa", "ob", "oc"]) for _ in range(int(rng.integers(4, 9)))] for _ in range(200)]
    null = S.positional_profile(make_tokens(null_lines), which="last", permutations=50, seed=1)
    assert null.pvalues.loc["line_first", "mid"] > 0.05


# ------------------------------------------------------------------ neighbour similarity


def _random_words(rng, n, length=8):
    seen, out = set(), []
    while len(out) < n:
        w = "".join(rng.choice(list("abcdefghijklmnopqrstuvwxyz"), size=length))
        if w not in seen:
            seen.add(w)
            out.append(w)
    return out


def test_neighbor_similarity_repeats_and_edits():
    rng = np.random.default_rng(0)
    words = _random_words(rng, 120)
    doubled = [w for w in words for _ in range(2)]  # each word twice in a row
    df = make_tokens([doubled[i : i + 8] for i in range(0, len(doubled), 8)], page_size=30)
    r = S.neighbor_similarity(df, window=5, max_edit=1, shuffles=10).iloc[0]
    assert r["share_ident"] == pytest.approx(0.5, abs=0.02)
    assert r["share_ident_shuffled"] < 0.06 and r["ratio_ident"] > 8 and r["z"] > 10
    # near-copies: same word with the last letter changed
    mutated = []
    for w in words:
        mutated += [w, w[:-1] + ("z" if w[-1] != "z" else "y")]
    df2 = make_tokens([mutated[i : i + 8] for i in range(0, len(mutated), 8)], page_size=30)
    r2 = S.neighbor_similarity(df2, window=5, max_edit=1, shuffles=10).iloc[0]
    assert r2["share_ident"] < 0.05 and r2["share_edit_only"] == pytest.approx(0.5, abs=0.03)
    r3 = S.neighbor_similarity(df2, window=5, max_edit=0, shuffles=3).iloc[0]
    assert r3["share_near"] < 0.05  # with max_edit=0 only identical tokens count


def test_neighbor_similarity_respects_pages():
    # the repeated word is at the start of the next page: not within the window of the same page
    df = make_tokens([["abcdefgh", "ijklmnop"], ["abcdefgh", "qrstuvwx"]], page_size=1)
    r = S.neighbor_similarity(df, window=5, shuffles=2).iloc[0]
    assert r["share_ident"] == 0.0 and r["n_tokens"] == 2


# ------------------------------------------------------------------ overlap and keyness


def test_vocabulary_overlap():
    df = pd.DataFrame({"token": list("aabc") + list("bcd"), "grp": ["A"] * 4 + ["B"] * 3})
    o = S.vocabulary_overlap(df, "grp", "A", "B")
    assert o["shared_types"] == 2 and o["jaccard"] == pytest.approx(2 / 4)
    assert o["a_token_in_b"] == pytest.approx(0.5) and o["a_type_in_b"] == pytest.approx(2 / 3)
    assert o["b_token_in_a"] == pytest.approx(2 / 3) and o["b_type_in_a"] == pytest.approx(2 / 3)
    assert o["n_tokens_a"] == 4 and o["n_types_b"] == 3


def test_labels_vs_text(synth_tables):
    o = S.labels_vs_text(synth_tables.tokens)
    assert o["n_tokens_a"] == 4 and o["n_tokens_b"] == 24 and o["shared_types"] == 0
    df = pd.DataFrame({"token": ["a", "b", "a", "c"], "generic_type": ["L", "L", "P", "P"]})
    o2 = S.labels_vs_text(df)
    assert o2["shared_types"] == 1 and o2["a_token_in_b"] == pytest.approx(0.5)


def test_keyness_sign_and_value():
    tokens = ["x"] * 30 + ["y"] * 5 + ["z"] * 10 + ["x"] * 5 + ["y"] * 30 + ["z"] * 10
    groups = ["A"] * 45 + ["B"] * 45
    df = pd.DataFrame({"token": tokens, "lang": groups})
    k = S.keyness(df, "lang", min_count=1)
    ax = k[(k["group"] == "A") & (k["item"] == "x")].iloc[0]
    bx = k[(k["group"] == "B") & (k["item"] == "x")].iloc[0]
    az = k[(k["group"] == "A") & (k["item"] == "z")].iloc[0]
    assert ax["sign"] == 1 and bx["sign"] == -1 and ax["signed_G2"] > 0 > bx["signed_G2"]
    expected = 2 * (30 * math.log(30 / 17.5) + 5 * math.log(5 / 17.5))
    assert ax["G2"] == pytest.approx(expected, rel=1e-9) == pytest.approx(19.81, abs=0.01)
    assert ax["G2"] == pytest.approx(bx["G2"])  # two equal-sized groups mirror each other
    assert az["G2"] == pytest.approx(0.0, abs=1e-12)
    assert ax["count_in"] == 30 and ax["count_out"] == 5 and ax["log2_ratio"] == pytest.approx(math.log2(6))
    top = S.top_keyness(k, k=1)
    assert top[top["group"] == "A"]["item"].tolist() == ["x"] and top[top["group"] == "B"]["item"].tolist() == ["y"]
    both = S.top_keyness(k, k=1, direction="both")
    assert len(both) == 4
    assert len(S.keyness(df, "lang", min_count=1000)) == 0


def test_keyness_three_groups_and_stem_requirement():
    df = pd.DataFrame({"token": list("aaabbbccc") * 3, "g": list("xxxyyyzzz") * 3})
    k = S.keyness(df, "g", min_count=1)
    assert set(k["group"]) == {"x", "y", "z"}
    assert k[(k["group"] == "x") & (k["item"] == "a")]["sign"].iloc[0] == 1
    with pytest.raises(ValueError):
        S.keyness(df, "g", unit="stem")
