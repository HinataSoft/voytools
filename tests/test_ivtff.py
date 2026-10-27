"""Parser and tokenizer tests on the synthetic IVTFF fixture."""

from __future__ import annotations

import pytest

from voytools.ivtff import (
    HIGH_ASCII_BASE,
    IVTFFError,
    decode_high_ascii,
    encode_high_ascii,
    is_high_ascii_char,
    parse_ivtff,
    parse_ivtff_text,
    render_high_ascii,
    tokenize_locus,
    tokenize_text,
)
from tests.conftest import SYNTH_INTERLINEAR, SYNTH_V101

HI150 = chr(HIGH_ASCII_BASE + 150)


def toks(locus, **kw):
    return [t.text for t in tokenize_locus(locus, **kw)]


def locus(synth, lid):
    return next(lc for lc in synth.loci if lc.id == lid)


# ------------------------------------------------------------------ structure


def test_header(synth):
    h = synth.header
    assert (h.alphabet, h.version, h.origin, h.extra) == ("Eva-", "2.0", "M", "1")
    assert h.version_tuple == (2, 0)
    assert h.is_eva_like


def test_pages_and_variables(synth):
    assert [p.name for p in synth.pages] == ["f1r", "f1v"]
    assert synth.pages[0].variables == dict(Q="A", P="A", F="a", B="1", I="T", L="A", H="1", C="1", X="V")
    assert synth.pages[1].variables["L"] == "@"
    assert "F" not in synth.pages[1].variables
    assert synth.pages[0].comments == ["first page"]


def test_locus_count_and_fields(synth):
    assert len(synth.loci) == 12
    lc = locus(synth, "f1r.1")
    assert (lc.folio, lc.num, lc.locator, lc.locus_type) == ("f1r", 1, "@", "P0")
    assert (lc.generic_type, lc.subtype, lc.transcriber) == ("P", "0", None)
    assert lc.paragraph_start and not lc.paragraph_end
    assert lc.comments == ["strange r"]
    assert lc.text == "fachys.ykal.ar.[cth:oto]res.y,kor.sholdy"
    assert lc.raw_text.startswith("<%>fachys")


def test_all_locators_and_types(synth):
    got = {lc.id: (lc.locator, lc.locus_type) for lc in synth.loci}
    assert got["f1r.1"] == ("@", "P0")
    assert got["f1r.2"] == ("+", "P0")
    assert got["f1r.3"] == ("*", "P0")
    assert got["f1r.4"] == ("=", "Pt")
    assert got["f1r.6"] == ("&", "Cc")
    assert got["f1r.7"] == ("~", "Lz")
    assert got["f1r.8"] == ("/", "Ls")
    assert got["f1r.9"] == ("!", "L0")
    assert got["f1v.1"] == ("@", "Ri")
    assert got["f1v.2"] == ("+", "Ro")
    assert {lc.generic_type for lc in synth.loci} == {"P", "C", "L", "R"}


def test_continuation_lines_are_joined(synth):
    lc = locus(synth, "f1r.3")
    assert lc.continued
    assert lc.text == "qokeedy.dai[n:r].ol,chedy.shy"
    assert lc.line_no == 8  # first physical line of the locus


def test_paragraph_end_and_comments(synth):
    assert locus(synth, "f1r.4").paragraph_end
    assert locus(synth, "f1r.5").paragraph_start and locus(synth, "f1r.5").paragraph_end
    assert locus(synth, "f1r.3").comments == ["x"]
    assert locus(synth, "f1r.8").comments == ["label"]
    assert locus(synth, "f1v.1").comments == ["10:00"]


def test_inline_text_tags_apply_to_following_loci(synth):
    assert locus(synth, "f1r.3").inline_tags == [("L", "B")]
    assert locus(synth, "f1r.2").variables["L"] == "A"
    assert locus(synth, "f1r.3").variables["L"] == "B"
    assert locus(synth, "f1r.9").variables["L"] == "B"  # persists until the end of the page
    # new page resets; L=@ means 'unset' until a tag sets it
    assert locus(synth, "f1v.1").variables["L"] is None
    assert locus(synth, "f1v.2").variables["L"] == "A"
    assert locus(synth, "f1v.3").variables["L"] == "A"


def test_no_warnings_for_clean_file(synth):
    assert synth.warnings == []


# ------------------------------------------------------------------ tokens


def test_tokens_split_policy(synth):
    lc = locus(synth, "f1r.1")
    ts = tokenize_locus(lc)
    assert [t.text for t in ts] == ["fachys", "ykal", "ar", "cthres", "y", "kor", "sholdy"]
    alt = ts[3]
    assert alt.raw == "[cth:oto]res" and alt.has_alternative
    assert ts[4].before_uncertain_space and not ts[4].after_uncertain_space
    assert ts[5].after_uncertain_space and not ts[5].before_uncertain_space
    assert not any(t.has_unreadable or t.has_ligature or t.has_high_ascii for t in ts)
    assert [t.index for t in ts] == list(range(7))


def test_tokens_join_policy(synth):
    ts = tokenize_locus(locus(synth, "f1r.1"), space_policy="join")
    assert [t.text for t in ts] == ["fachys", "ykal", "ar", "cthres", "ykor", "sholdy"]
    assert not any(t.after_uncertain_space or t.before_uncertain_space for t in ts)


def test_tokens_strict_policy(synth):
    ts = tokenize_locus(locus(synth, "f1r.1"), space_policy="strict")
    assert [t.text for t in ts] == ["fachys", "ykal", "ar", "cthres", "y,kor", "sholdy"]
    assert ts[4].raw == "y,kor" and ts[4].n_uncertain_inside == 1
    assert not any(t.after_uncertain_space or t.before_uncertain_space for t in ts)


def test_special_characters_flags_and_drawing_breaks(synth):
    ts = tokenize_locus(locus(synth, "f1r.2"))
    assert [t.text for t in ts] == ["sory", "ckhy", HI150 + "ar", "o?x", "???", "d", "kor", "a", "b"]
    lig, hi, unr, unk = ts[1], ts[2], ts[3], ts[4]
    assert lig.has_ligature and lig.raw == "{ck}hy" and not lig.has_high_ascii
    assert hi.has_high_ascii and hi.raw == "@150;ar" and len(hi.text) == 3
    assert unr.has_unreadable and unk.has_unreadable
    assert ts[5].before_drawing_break and not ts[5].after_drawing_break
    assert ts[6].after_drawing_break
    assert ts[7].before_drawing_break and ts[8].after_drawing_break  # <~>
    # drawing breaks split in every policy
    for policy in ("split", "join", "strict"):
        assert toks(locus(synth, "f1r.2"), space_policy=policy)[5:] == ["d", "kor", "a", "b"]


def test_alternatives_mark_policy(synth):
    ts = tokenize_locus(locus(synth, "f1r.1"), alternatives="mark")
    assert ts[3].text == "[cth:oto]res" and ts[3].has_alternative
    ts3 = tokenize_locus(locus(synth, "f1r.3"), alternatives="mark")
    assert [t.text for t in ts3][1] == "dai[n:r]"
    # default keeps the first option
    assert [t.text for t in tokenize_locus(locus(synth, "f1r.3"))][1] == "dain"


def test_continued_locus_tokens(synth):
    ts = tokenize_locus(locus(synth, "f1r.3"))
    assert [t.text for t in ts] == ["qokeedy", "dain", "ol", "chedy", "shy"]
    assert ts[2].before_uncertain_space and ts[3].after_uncertain_space


def test_separator_runs_and_edges():
    # leading/trailing separators ignored; strongest separator of a run wins
    ts = tokenize_text(".a..b,.c,d,", "split")
    assert [t.text for t in ts] == ["a", "b", "c", "d"]
    assert not ts[1].after_uncertain_space  # '..' = definite
    assert not ts[2].after_uncertain_space  # ',.' = definite
    assert ts[3].after_uncertain_space and ts[2].before_uncertain_space
    assert not ts[3].before_uncertain_space  # trailing comma has no following token
    ts = tokenize_text("a.<->b", "split")
    assert ts[1].after_drawing_break and not ts[1].after_uncertain_space


def test_alternative_with_ligature_and_high_ascii():
    ts = tokenize_text("dai[{cto}:@194;]y")
    assert len(ts) == 1
    assert ts[0].text == "daictoy" and ts[0].has_ligature and ts[0].has_alternative
    ts = tokenize_text("dai[@194;:{cto}]y")
    assert ts[0].has_high_ascii and ts[0].text == "dai" + chr(HIGH_ASCII_BASE + 194) + "y"


def test_bad_policy_raises():
    with pytest.raises(ValueError):
        tokenize_text("a.b", space_policy="nonsense")
    with pytest.raises(ValueError):
        tokenize_text("a.b", alternatives="nonsense")


# ------------------------------------------------------------------ high ASCII helpers


def test_high_ascii_helpers():
    ch = encode_high_ascii(150)
    assert ch == HI150 and is_high_ascii_char(ch) and not is_high_ascii_char("a")
    assert render_high_ascii("x" + ch + "y") == "x@150;y"
    assert decode_high_ascii("x@150;y") == "x" + ch + "y"
    assert decode_high_ascii("@099;") == "@099;"  # out of range stays literal
    with pytest.raises(ValueError):
        encode_high_ascii(100)


# ------------------------------------------------------------------ robustness


def test_malformed_locus_identifier_raises():
    base = "#=IVTFF Eva- 2.0 M\n<f1r>\n"
    for bad in ["<f1r.x,+P0> abc", "<f1r.1,+p0> abc", "<f1r.1+P0> abc", "<f1r.1,+P> abc", "<f1r.1,+P0 abc"]:
        with pytest.raises(IVTFFError):
            parse_ivtff_text(base + bad + "\n")
    with pytest.raises(IVTFFError):
        parse_ivtff_text(base + "<f1r.1, +P0> abc\n")  # whitespace inside identifier


def test_unknown_constructs_collect_warnings_not_errors():
    text = (
        "#=IVTFF Eva- 2.0 M\n<f1r>\n"
        "<f1r.1,@P0>   ab<?foo>c.d[e.f{g.h\n"
        "<f1r.2,+Q0>   x\n"
        "<f2r.1,+P0>   y\n"
        "stray line\n"
        "/ orphan continuation\n"
        "<AA001>   alt\n"
    )
    t = parse_ivtff_text(text)
    cats = {w.category for w in t.warnings}
    assert {"unknown-markup", "unbalanced-bracket", "unknown-locus-type", "page-mismatch",
            "unexpected-line", "orphan-continuation", "alternative-locus"} <= cats
    assert len(t.loci) == 3
    assert [x.text for x in tokenize_locus(t.loci[0])][0] == "abc"


def test_missing_header_is_a_warning():
    t = parse_ivtff_text("<f1r>\n<f1r.1,@P0> ab.cd\n")
    assert t.header.alphabet is None
    assert any(w.category == "missing-header" for w in t.warnings)
    assert len(t.loci) == 1


def test_unclosed_identifier_raises():
    with pytest.raises(IVTFFError):
        parse_ivtff_text("#=IVTFF Eva- 2.0 M\n<f1r.1,@P0 abc\n")


def test_parse_file_roundtrip(synth_path, synth):
    t = parse_ivtff(synth_path)
    assert len(t.loci) == len(synth.loci) and t.path == str(synth_path)


# ------------------------------------------------------------------ interlinear / old files


def test_interlinear_transcribers_and_fillers():
    t = parse_ivtff_text(SYNTH_INTERLINEAR)
    assert t.transcribers == ["C", "H", "U", "Z"]
    assert t.is_interlinear
    h1 = [lc for lc in t.loci if lc.id == "f1r.1" and lc.transcriber == "H"][0]
    assert h1.raw_text == "fachys.ykal.ar!.cth%res"
    assert h1.text == "fachys.ykal.ar.cthres"
    c1 = [lc for lc in t.loci if lc.id == "f1r.1" and lc.transcriber == "C"][0]
    assert toks(c1) == ["fachys", "ykal", "ar", "cthorys"]
    h2 = [lc for lc in t.loci if lc.id == "f1r.2" and lc.transcriber == "H"][0]
    assert toks(h2) == ["sory", "ckhar", "or"]
    u2 = [lc for lc in t.loci if lc.transcriber == "U"][0]
    assert toks(u2) == []
    z2 = [lc for lc in t.loci if lc.transcriber == "Z"][0]
    assert z2.paragraph_end
    odd = t.loci[-1]
    assert (odd.num, odd.num_suffix, odd.id) == (44, "a", "f2r.44a")
    assert [lc.id for lc in t.iter_loci(transcriber="H")] == ["f1r.1", "f1r.2", "f2r.44a"]


def test_fillers_kept_when_disabled_or_v101():
    t = parse_ivtff_text(SYNTH_INTERLINEAR, strip_fillers=False)
    h1 = [lc for lc in t.loci if lc.transcriber == "H"][0]
    assert toks(h1) == ["fachys", "ykal", "ar!", "cth%res"]
    t = parse_ivtff_text(SYNTH_V101)
    assert t.header.alphabet == "v101"
    assert toks(t.loci[0]) == ["fa!s", "%a9", "hae", "ch9"]  # ! and % are v101 letters
    assert t.loci[0].paragraph_start and t.loci[0].paragraph_end
