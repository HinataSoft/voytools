"""Token / locus / page table tests (synthetic fixture)."""

from __future__ import annotations

import pandas as pd
import pytest

from voytools.ivtff import HIGH_ASCII_BASE, parse_ivtff_text
from voytools.table import TOKEN_COLUMNS, build_tables, export, filter_tokens, most_complete_transcriber
from tests.conftest import SYNTH_INTERLINEAR

HI150 = chr(HIGH_ASCII_BASE + 150)

EXPECTED_COLUMNS = [
    "idx", "folio", "quire", "page_in_quire", "folio_in_quire", "bifolio", "illustration", "language",
    "hand", "currier_hand", "extraneous", "locus_id", "locus_num", "locator", "locus_type",
    "generic_type", "subtype", "transcriber", "par_index", "line_in_par", "is_par_first_line",
    "is_par_last_line", "line_index", "token_index", "n_tokens", "pos_from_end", "rel_pos",
    "is_first", "is_last", "token", "token_raw", "len_chars", "len_glyphs", "has_unreadable",
    "has_alternative", "has_ligature", "has_high_ascii", "after_uncertain_space",
    "before_uncertain_space", "after_drawing_break", "before_drawing_break",
]


def col(df, name):
    return df[name].tolist()


def test_exact_column_names(synth_tables):
    assert list(synth_tables.tokens.columns) == EXPECTED_COLUMNS == TOKEN_COLUMNS


def test_token_counts_and_order(synth_tables):
    t = synth_tables.tokens
    assert len(t) == 34
    assert col(t, "idx") == list(range(34))
    assert t.groupby("locus_id", sort=False)["n_tokens"].first().tolist() == [7, 9, 5, 1, 2, 2, 1, 1, 1, 2, 2, 1]
    assert len(synth_tables.loci) == 12 and len(synth_tables.pages) == 2


def test_page_variables_and_text_tags(synth_tables):
    t = synth_tables.tokens
    first = t.iloc[0]
    assert (first["folio"], first["quire"], first["page_in_quire"], first["folio_in_quire"]) == ("f1r", "A", "A", "a")
    assert (first["bifolio"], first["illustration"], first["hand"], first["currier_hand"], first["extraneous"]) == (
        "1", "T", "1", "1", "V")
    per_locus = t.groupby("locus_id", sort=False)["language"].first().tolist()
    assert per_locus[:2] == ["A", "A"]
    assert per_locus[2:9] == ["B"] * 7  # text tag <@L=B> on f1r.3 applies from there on
    assert pd.isna(per_locus[9])  # f1v: $L=@ and no tag yet
    assert per_locus[10:] == ["A", "A"]
    v = t[t["folio"] == "f1v"].iloc[0]
    assert pd.isna(v["folio_in_quire"]) and pd.isna(v["bifolio"]) and pd.isna(v["currier_hand"]) and pd.isna(v["extraneous"])
    assert (v["quire"], v["page_in_quire"], v["illustration"]) == ("A", "B", "H")


def test_position_columns(synth_tables):
    t = synth_tables.tokens
    l1 = t[t["locus_id"] == "f1r.1"]
    assert col(l1, "token_index") == list(range(7))
    assert col(l1, "pos_from_end") == list(range(6, -1, -1))
    assert col(l1, "rel_pos") == pytest.approx([i / 6 for i in range(7)])
    assert col(l1, "is_first") == [True] + [False] * 6
    assert col(l1, "is_last") == [False] * 6 + [True]
    single = t[t["locus_id"] == "f1r.4"]
    assert col(single, "rel_pos") == [0.0] and col(single, "is_first") == [True] and col(single, "is_last") == [True]
    assert col(t[t["folio"] == "f1v"].drop_duplicates("locus_id"), "line_index") == [0, 1, 2]
    assert col(t[t["folio"] == "f1r"].drop_duplicates("locus_id"), "line_index") == list(range(9))


def test_paragraph_columns(synth_tables):
    loci = synth_tables.loci.set_index("locus_id")
    assert loci.loc["f1r.1", ["par_index", "line_in_par"]].tolist() == [0, 0]
    assert loci.loc["f1r.2", ["par_index", "line_in_par"]].tolist() == [0, 1]
    assert loci.loc["f1r.3", ["par_index", "line_in_par"]].tolist() == [0, 2]
    assert loci.loc["f1r.4", ["par_index", "line_in_par"]].tolist() == [0, 3]
    assert loci.loc["f1r.5", ["par_index", "line_in_par"]].tolist() == [1, 0]
    assert loci.loc["f1r.1", "is_par_first_line"] and not loci.loc["f1r.2", "is_par_first_line"]
    assert loci.loc["f1r.4", "is_par_last_line"] and not loci.loc["f1r.3", "is_par_last_line"]
    assert loci.loc["f1r.5", "is_par_first_line"] and loci.loc["f1r.5", "is_par_last_line"]
    for lid in ("f1r.6", "f1r.7", "f1r.8", "f1r.9", "f1v.1", "f1v.3"):
        assert pd.isna(loci.loc[lid, "par_index"]) and pd.isna(loci.loc[lid, "line_in_par"])
        assert not loci.loc[lid, "is_par_first_line"] and not loci.loc[lid, "is_par_last_line"]


def test_token_strings_lengths_and_flags(synth_tables):
    t = synth_tables.tokens
    r = t[t["locus_id"] == "f1r.1"].set_index("token_index")
    assert r.loc[3, ["token", "token_raw", "len_chars", "len_glyphs"]].tolist() == ["cthres", "[cth:oto]res", 6, 4]
    assert r.loc[0, ["token", "len_chars", "len_glyphs"]].tolist() == ["fachys", 6, 5]  # f a ch y s
    assert r.loc[6, "len_glyphs"] == 5  # sh o l d y
    assert bool(r.loc[3, "has_alternative"]) and not bool(r.loc[0, "has_alternative"])
    assert bool(r.loc[4, "before_uncertain_space"]) and bool(r.loc[5, "after_uncertain_space"])
    r2 = t[t["locus_id"] == "f1r.2"].set_index("token_index")
    assert r2.loc[2, "token"] == HI150 + "ar" and bool(r2.loc[2, "has_high_ascii"]) and r2.loc[2, "len_chars"] == 3
    assert bool(r2.loc[1, "has_ligature"]) and r2.loc[1, "token"] == "ckhy" and r2.loc[1, "len_glyphs"] == 2
    assert bool(r2.loc[3, "has_unreadable"]) and bool(r2.loc[4, "has_unreadable"])
    assert bool(r2.loc[5, "before_drawing_break"]) and bool(r2.loc[6, "after_drawing_break"])
    assert bool(r2.loc[7, "before_drawing_break"]) and bool(r2.loc[8, "after_drawing_break"])
    assert not r2.loc[5, ["after_drawing_break", "before_uncertain_space"]].any()


def test_space_policies_in_table(synth):
    join = build_tables(synth, space_policy="join").tokens
    split = build_tables(synth, space_policy="split").tokens
    strict = build_tables(synth, space_policy="strict").tokens
    assert len(split) == 34 and len(join) == 31 and len(strict) == 31
    assert "ykor" in join["token"].tolist() and "y,kor" in strict["token"].tolist()
    # strict: the comma marker does not count as a character
    row = strict[strict["token"] == "y,kor"].iloc[0]
    assert row["len_chars"] == 4 and row["len_glyphs"] == 4
    assert not strict["after_uncertain_space"].any() and not strict["before_uncertain_space"].any()
    with pytest.raises(TypeError):
        build_tables(synth, bogus=1)


def test_alphabet_dependent_glyph_length():
    t = parse_ivtff_text("#=IVTFF v101 2.0 M\n<f1r>\n<f1r.1,@P0> ch9.sh9\n")
    df = build_tables(t).tokens
    assert df["len_glyphs"].tolist() == df["len_chars"].tolist() == [3, 3]
    assert df.attrs["alphabet"] == "v101" and df.attrs["glyph_preset"] is None
    t2 = parse_ivtff_text("#=IVTFF Eva- 2.0 M\n<f1r>\n<f1r.1,@P0> chsh.cthol\n")
    assert build_tables(t2).tokens["len_glyphs"].tolist() == [2, 3]
    assert build_tables(t2, glyph_preset=None).tokens["len_glyphs"].tolist() == [4, 5]
    assert build_tables(parse_ivtff_text("#=IVTFF Eva- 2.0 M\n<f1r>\n<f1r.1,@P0> daiin\n"),
                        glyph_preset="extended").tokens["len_glyphs"].tolist() == [3]


def test_infer_paragraphs():
    text = ("#=IVTFF Eva- 2.0 M\n<f1r>\n<f1r.1,@P0> a.b\n<f1r.2,+P0> c.d\n"
            "<f1r.3,*P0> e.f\n<f1r.4,+P0> g\n<f1r.5,@L0> lab\n")
    t = parse_ivtff_text(text)
    plain = build_tables(t).loci
    assert plain["par_index"].isna().all()
    inf = build_tables(t, infer_paragraphs=True).loci
    assert inf["par_index"].tolist()[:4] == [0, 0, 1, 1] and pd.isna(inf["par_index"].iloc[4])
    assert inf["line_in_par"].tolist()[:4] == [0, 1, 0, 1]
    assert inf["is_par_last_line"].tolist() == [False, True, False, True, False]
    assert inf["is_par_first_line"].tolist() == [True, False, True, False, False]


def test_interlinear_tables():
    t = parse_ivtff_text(SYNTH_INTERLINEAR)
    allt = build_tables(t)
    assert allt.transcribers == ["C", "H", "U", "Z"]
    assert len(allt.tokens) == 20  # H: 4+3+2, C: 4+4, Z: 3, U: empty line
    h = build_tables(t, transcriber="H")
    assert len(h.tokens) == 9 and set(h.tokens["transcriber"]) == {"H"}
    assert h.tokens["locus_id"].drop_duplicates().tolist() == ["f1r.1", "f1r.2", "f2r.44a"]
    assert h.tokens[h.tokens["locus_id"] == "f2r.44a"]["locus_num"].iloc[0] == 44
    c = build_tables(t, transcriber="C").loci
    assert c["line_index"].tolist() == [0, 1]
    assert most_complete_transcriber(allt.loci) == "H"
    # all transcribers at once: own line numbering per transcriber
    li = allt.loci.set_index(["locus_id", "transcriber"])["line_index"]
    assert li[("f1r.2", "C")] == 1 and li[("f1r.2", "H")] == 1  # both also have f1r.1
    assert li[("f1r.2", "U")] == 0 and li[("f1r.2", "Z")] == 0  # these transcribers only have f1r.2
    assert allt.loci.set_index(["locus_id", "transcriber"]).loc[("f1r.2", "U"), "n_tokens"] == 0
    assert len(h.pages) == 2


def test_filter_tokens(synth_tables):
    t = synth_tables.tokens
    p = filter_tokens(t, generic_type="P", clean=False)
    assert len(p) == 24 and p.index.tolist() == list(range(24)) and p["idx"].iloc[-1] == 23
    assert len(filter_tokens(t, generic_type="P")) == 19  # 2 unreadable, 2 alternative, 1 high-ASCII dropped
    assert len(filter_tokens(t, generic_type="P", clean=("ligature",))) == 23
    assert len(filter_tokens(t, generic_type="P", clean=("unreadable",))) == 22
    assert len(filter_tokens(t, generic_type="P,L", clean=False)) == 24 + 4
    assert len(filter_tokens(t, generic_type=["C", "R"], clean=False)) == 6
    assert len(filter_tokens(t, language="B", clean=False)) == 13
    assert len(filter_tokens(t, language="A", clean=False)) == 19
    assert len(filter_tokens(t, language=["A", "B"], clean=False)) == 32  # NA never matches
    assert len(filter_tokens(t, illustration="H", clean=False)) == 5
    assert len(filter_tokens(t, hand=1, clean=False)) == 34 and len(filter_tokens(t, hand="2", clean=False)) == 0
    assert len(filter_tokens(t, generic_type="all", clean=False)) == 34
    assert len(filter_tokens(t, locus_type="Pt", clean=False)) == 1
    assert len(filter_tokens(t, folio="f1v", language="A", clean=False)) == 3
    assert filter_tokens(t, generic_type="P").attrs["alphabet"] == "Eva-"
    with pytest.raises(ValueError):
        filter_tokens(t, clean=("nonsense",))


def test_filter_does_not_modify_input(synth_tables):
    before = synth_tables.tokens.copy()
    out = filter_tokens(synth_tables.tokens, generic_type="P")
    out["token"] = "zzz"
    pd.testing.assert_frame_equal(before, synth_tables.tokens)


def test_loci_and_pages_tables(synth_tables):
    loci, pages = synth_tables.loci, synth_tables.pages
    assert loci["comments"].iloc[0] == "strange r" and loci["inline_tags"].iloc[2] == "L=B"
    assert bool(loci["continued"].iloc[2]) and loci["n_tokens"].sum() == 34
    assert pages["folio"].tolist() == ["f1r", "f1v"]
    assert pages["n_loci"].tolist() == [9, 3] and pages["n_loci_P"].tolist() == [5, 0]
    assert pages["n_tokens"].tolist() == [29, 5]
    assert pages.loc[1, "language"] == "@"  # page variables are kept as written


def test_export_csv_and_parquet_fallback(synth_tables, tmp_path):
    p = export(synth_tables.tokens, tmp_path / "t.csv")
    back = pd.read_csv(p, dtype=str, keep_default_na=False)
    assert list(back.columns) == EXPECTED_COLUMNS
    assert "@150;ar" in back["token"].tolist() and len(back) == 34
    raw = export(synth_tables.tokens, tmp_path / "raw.csv", render_high_ascii_codes=False)
    assert HI150 + "ar" in raw.read_text(encoding="utf-8")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        with pytest.warns(RuntimeWarning):
            q = export(synth_tables.tokens, tmp_path / "t.parquet")
        assert q.suffix == ".csv" and q.exists()
    else:  # pragma: no cover
        q = export(synth_tables.tokens, tmp_path / "t.parquet")
        assert q.suffix == ".parquet"
    paths = synth_tables.export(tmp_path / "all")
    assert [x.name for x in paths] == ["tokens.csv", "loci.csv", "pages.csv"]
