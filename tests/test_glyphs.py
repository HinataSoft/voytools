"""Glyph collapsing and symbol helpers."""

from __future__ import annotations

import pytest

from voytools.glyphs import (
    EVA_BASIC,
    EVA_EXTENDED,
    collapse_eva,
    default_preset,
    glyph_sequence,
    is_eva_alphabet,
    resolve_patterns,
    to_symbol_string,
    to_symbol_strings,
)
from tests.conftest import make_tokens


def test_presets_content():
    assert EVA_BASIC == ("cth", "ckh", "cph", "cfh", "ch", "sh")
    assert EVA_EXTENDED == EVA_BASIC + ("iin", "in", "ee")
    assert resolve_patterns(None) == () and resolve_patterns("none") == ()
    with pytest.raises(ValueError):
        resolve_patterns("nonsense")


@pytest.mark.parametrize(
    "token, preset, expected",
    [
        ("chedy", "basic", ("ch", "e", "d", "y")),
        ("cthol", "basic", ("cth", "o", "l")),
        ("ckhhy", "basic", ("ckh", "h", "y")),
        ("cphy", "basic", ("cph", "y")),
        ("cfhaiin", "basic", ("cfh", "a", "i", "i", "n")),
        ("shsh", "basic", ("sh", "sh")),
        ("chch", "basic", ("ch", "ch")),
        ("sch", "basic", ("s", "ch")),
        ("daiin", "extended", ("d", "a", "iin")),
        ("daiiin", "extended", ("d", "a", "i", "iin")),
        ("ain", "extended", ("a", "in")),
        ("qokeedy", "extended", ("q", "o", "k", "ee", "d", "y")),
        ("qokeeedy", "extended", ("q", "o", "k", "ee", "e", "d", "y")),
        ("cheey", "extended", ("ch", "ee", "y")),
        ("daiin", "basic", ("d", "a", "i", "i", "n")),
        ("daiin", None, tuple("daiin")),
        ("", "basic", ()),
        ("o?x", "basic", ("o", "?", "x")),
        ("or,y", "basic", ("o", "r", "y")),  # uncertain-space marker ignored
        ("oteo[s:r]res", "basic", ("o", "t", "e", "o", "[s:r]", "r", "e", "s")),
    ],
)
def test_collapse_eva(token, preset, expected):
    assert collapse_eva(token, preset) == expected


def test_collapse_custom_patterns():
    assert collapse_eva("qokal", ["qo", "ka"]) == ("qo", "ka", "l")


def test_collapse_keeps_high_ascii_as_one_glyph():
    token = "ar"
    assert collapse_eva(token) == ("", "a", "r")


def test_alphabet_helpers():
    assert is_eva_alphabet("Eva-") and is_eva_alphabet("EvaT") and not is_eva_alphabet("v101")
    assert not is_eva_alphabet(None)
    assert default_preset("Eva-") == "basic" and default_preset("v101") is None


def test_glyph_sequence_per_line_and_boundary():
    df = make_tokens([["chedy", "qol"], ["daiin"]])
    seq = glyph_sequence(df, "basic")
    assert seq == [["ch", "e", "d", "y", "q", "o", "l"], ["d", "a", "i", "i", "n"]]
    seq = glyph_sequence(df, "extended", separator="_")
    assert seq == [["ch", "e", "d", "y", "_", "q", "o", "l"], ["d", "a", "iin"]]
    flat = glyph_sequence(df, "basic", separator="_", scope="all")
    assert len(flat) == 1 and flat[0].count("_") == 2 and flat[0][:5] == ["ch", "e", "d", "y", "_"]
    assert glyph_sequence(df, None)[1] == list("daiin")
    assert glyph_sequence(df)[0][0] == "ch"  # 'auto' follows the table (Eva -> basic)
    with pytest.raises(ValueError):
        glyph_sequence(df, "basic", scope="nonsense")


def test_to_symbol_string():
    s, mapping = to_symbol_string(["ch", "e", "d", "ch", "e"])
    assert len(s) == 5 and len(set(s)) == 3
    assert mapping["e"] == "e" and mapping["d"] == "d"  # single characters map to themselves
    assert len(mapping["ch"]) == 1 and mapping["ch"] not in "ed"
    s2, mapping2 = to_symbol_string(["e", "ch", "x"], mapping)
    assert mapping2["ch"] == mapping["ch"] and s2[0] == "e" and s2[1] == mapping["ch"]
    many, m = to_symbol_strings([["ch", "a"], ["a", "sh", "ch"]])
    assert [len(x) for x in many] == [2, 3] and len(set("".join(many))) == 3
    assert to_symbol_string([]) == ("", {})
