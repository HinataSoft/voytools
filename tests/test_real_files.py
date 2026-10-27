"""Smoke and regression tests on the real transliteration files (skipped when absent)."""

from __future__ import annotations

import pytest

from voytools import stats as S
from voytools.glyphs import EVA_EXTENDED, collapse_eva
from voytools.ivtff import parse_ivtff
from voytools.morphology import decompose_table, paradigm_metrics
from voytools.table import build_tables, filter_tokens, load_tables, most_complete_transcriber


@pytest.fixture(scope="module")
def zl_clean(zl_tables):
    return filter_tokens(zl_tables.tokens, generic_type="P")  # paragraph text, clean tokens


# ------------------------------------------------------------------ ZL3b


def test_zl_structure(zl_path, zl_tables):
    t = zl_tables
    assert (t.header.alphabet, t.header.version) == ("Eva-", "2.0")
    assert len(t.pages) == 227 and len(t.loci) == 5385
    assert t.warnings == []
    assert t.tokens.attrs["glyph_preset"] == "basic"
    n = len(t.tokens)
    assert 34_000 <= n <= 40_000
    n_p = int((t.tokens["generic_type"] == "P").sum())
    assert 34_000 <= n_p <= 36_500
    p0 = filter_tokens(t.tokens, generic_type="P", locus_type="P0")
    assert 31_000 <= len(p0) <= 34_000
    assert set(t.tokens["generic_type"]) == {"P", "C", "L", "R"}
    assert set(t.tokens["language"].dropna()) == {"A", "B"}
    assert t.tokens["token"].map(len).min() >= 1 and not t.tokens["token"].str.contains(r"[\s.]", regex=True).any()


def test_zl_exact_counts_for_version_3b(zl_path, zl_tables):
    """Regression guard: exact counts of ZL3b-n.txt (Version 3b of 13/05/2025)."""
    translit = parse_ivtff(zl_path)
    if not any("Version 3b of 13/05/2025" in text for _, text in translit.comment_lines):
        pytest.skip("different edition of the ZL file")
    assert len(zl_tables.tokens) == 39_026
    assert int((zl_tables.tokens["generic_type"] == "P").sum()) == 35_075
    assert len(filter_tokens(zl_tables.tokens, generic_type="P")) == 34_268


def test_zl_loads_deterministically(zl_path, zl_tables):
    again = load_tables(zl_path)
    assert again.tokens.equals(zl_tables.tokens) and again.loci.equals(zl_tables.loci)


def test_zl_token_counts_by_language_add_up(zl_tables):
    clean = filter_tokens(zl_tables.tokens, generic_type="P")
    by = S.vocabulary_summary(clean, by="language").set_index("language")
    assert by["n_tokens"].sum() == len(clean) and {"A", "B"} <= set(by.index)
    assert by.loc["B", "n_types"] > 1000 and 0.3 < by.loc["A", "hapax_share"] < 0.8


def test_collapse_eva_equals_sequential_replacement_on_all_zl_tokens(zl_tables):
    marks = {pat: chr(0x2460 + i) for i, pat in enumerate(EVA_EXTENDED)}
    back = {m: p for p, m in marks.items()}
    order = ["cth", "ckh", "cph", "cfh", "ch", "sh", "iin", "in", "ee"]
    assert tuple(order) == EVA_EXTENDED
    for tok in zl_tables.tokens["token"].unique():
        s = tok
        for pat in order:
            s = s.replace(pat, marks[pat])
        expected = tuple(back.get(ch, ch) for ch in s)
        assert collapse_eva(tok, "extended") == expected, tok


def test_published_glyph_entropy_values(zl_clean):
    """Rozanova & Temerev (2026): glyph-level conditional entropy about 2.69 bits."""
    ext = S.glyph_entropy_report(zl_clean, "extended", scope="all").iloc[0]
    assert 2.59 <= ext["H2"] <= 2.79
    assert 3.8 <= ext["H1"] <= 4.1
    basic = S.glyph_entropy_report(zl_clean, "basic", scope="all").iloc[0]
    assert basic["H2"] < ext["H2"] and basic["H1"] < ext["H1"]
    by_line = S.glyph_entropy_report(zl_clean, "extended", scope="line").iloc[0]
    assert by_line["H2"] <= ext["H2"] + 1e-9


def test_published_mi_values(zl_clean):
    """Adjacent-token MI excess share below 1-2 percent; edge-glyph MI about 0.2 bits."""
    tok = S.adjacent_token_mi(zl_clean, shuffles=20, seed=0).iloc[0]
    assert 0 < tok["excess_share"] < 0.02
    edge = S.edge_glyph_mi(zl_clean, shuffles=20, seed=0).iloc[0]
    assert 0.1 <= edge["excess_bits"] <= 0.3 and edge["z"] > 10


def test_positional_asymmetry_of_gallows_at_paragraph_starts(zl_clean):
    prof = S.positional_profile(zl_clean, which="first")
    gallows_first = prof.shares.loc[["p", "t", "k", "f"], "parfirst_line_first"].sum()
    gallows_mid = prof.shares.loc[["p", "t", "k", "f"], "mid"].sum()
    assert gallows_first > 3 * gallows_mid
    assert prof.jsd.loc["parfirst_line_first", "mid"] > 0.1


def test_morphology_pipeline_on_zl(zl_clean):
    dfd = decompose_table(zl_clean, prefixes="reference", suffixes="reference")
    assert ((dfd["prefix"] + dfd["stem"] + dfd["suffix"]) == dfd["token"]).all()
    pm = paradigm_metrics(dfd, min_stem_count=5, shuffles=3, seed=0)
    assert pm.summary["n_stems_k"] > 50 and 0 <= pm.summary["cmi_prefix_suffix_given_stem"] < 1.5


# ------------------------------------------------------------------ other files


def test_gc2a_v101_alphabet(gc_path):
    t = load_tables(gc_path)
    assert t.header.alphabet == "v101" and t.tokens.attrs["glyph_preset"] is None
    assert (t.tokens["len_glyphs"] == t.tokens["len_chars"]).all()
    assert t.tokens["token"].str.contains("!", regex=False).any()  # '!' is a letter in v101
    assert t.tokens["has_high_ascii"].any()
    assert t.warnings == []
    clean = filter_tokens(t.tokens, generic_type="P")
    assert not clean["has_high_ascii"].any() and 30_000 < len(clean) < 40_000
    ent = S.glyph_entropy_report(clean, None).iloc[0]
    assert ent["H1"] > 3.5 and ent["H2"] < ent["H1"]


def test_rf1b_and_it2a_parse(it_path):
    from tests.conftest import data_dir

    rf = data_dir() / "RF1b-er.txt"
    if not rf.exists():
        pytest.skip("RF1b-er.txt not found")
    t = load_tables(rf)
    assert t.header.alphabet == "Eva-" and t.warnings == [] and len(t.pages) == 227
    assert 36_000 <= len(t.tokens) <= 40_000
    t = load_tables(it_path)
    assert t.header.alphabet == "EvaT" and t.warnings == [] and len(t.pages) == 225
    assert len(t.tokens) == 37_919


def test_lsi_interlinear(lsi_path, it_path):
    translit = parse_ivtff(lsi_path)
    assert translit.is_interlinear and len(translit.transcribers) == 18
    assert translit.header.version == "1.5"
    best = most_complete_transcriber(translit)
    assert best == "H"
    h = build_tables(translit, transcriber="H")
    it = load_tables(it_path)
    # the Takeshi Takahashi transcription is the H transcriber of the interlinear file
    assert len(h.tokens) == len(it.tokens) == 37_919
    assert (h.tokens.groupby("folio", sort=False).size().sort_index().to_dict()
            == it.tokens.groupby("folio", sort=False).size().sort_index().to_dict())
    assert not h.tokens["token"].str.contains("[!%]", regex=True).any()  # old-style fillers removed
    allt = build_tables(translit)
    assert set(allt.tokens["transcriber"].dropna()) >= {"H", "C", "F"}
    assert len(allt.tokens) > 2 * len(h.tokens)  # 18 partial transcriptions of the same text


def test_space_policies_on_zl(zl_path):
    split = load_tables(zl_path, space_policy="split").tokens
    join = load_tables(zl_path, space_policy="join").tokens
    strict = load_tables(zl_path, space_policy="strict").tokens
    assert len(split) > len(join) and len(strict) == len(join)  # strict keeps the join segmentation ...
    assert strict["token"].str.contains(",", regex=False).any()  # ... but remembers the uncertain space
    assert not split["token"].str.contains(",", regex=False).any()
    assert not join["token"].str.contains(",", regex=False).any()
    assert (strict["len_chars"] == join["len_chars"]).all()  # the marker is not counted as a character
