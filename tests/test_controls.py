"""Control texts: parsing, division selection and layout in a template (voytools.controls)."""

from __future__ import annotations

import numpy as np
import pytest

from voytools import controls as C
from voytools.stats import vocabulary_summary
from tests.conftest import make_tokens

GUT = """header junk
*** START OF THIS PROJECT GUTENBERG EBOOK X ***
body one
*** END OF THIS PROJECT GUTENBERG EBOOK X ***
licence
"""

LINN = """*** START OF THE PROJECT GUTENBERG EBOOK SPECIES PLANTARUM ***
preface words

_Classis I._

MONANDRIA

_MONOGYNIA._

CANNA.

[_Indica._]

1. Canna foliis ovatis. _Roy. lugdb. 11._

  Habitat in _India_. ♃

_DIGYNIA._

BLITUM.

1. Blitum capitellis.

_Classis II._

DIANDRIA

_MONOGYNIA._

OLEA.

1. Olea foliis lanceolatis.

       *       *       *       *       *

AUTHORS AND TEXTS
Alpini, Prospero.
"""

CULP = """intro

    AMARA DULCIS.

CONSIDERING divers shires. It grows. Well.

    ASH TREE.

_Government._] It is governed by the Sun; good.

    WILD BETONY, CALLED
    WOOD BETONY.

_Descript._] Leaves. Stalks.

DIRECTIONS FOR MAKING SYRUPS
more text
"""


def _template(n_lines=40, page_size=8, seed=0):
    rng = np.random.default_rng(seed)
    vocab = ["qokeedy", "chedy", "daiin", "ol", "or", "ar", "chol", "dar"]
    lines = [rng.choice(vocab, size=int(rng.integers(3, 9))).tolist() for _ in range(n_lines)]
    return make_tokens(lines, page_size=page_size)


def _divs(sizes, item_len=5):
    out = []
    for d, n in enumerate(sizes):
        toks = [f"w{d}x{i % 37}" for i in range(n)]
        out.append((f"d{d}", [toks[i: i + item_len] for i in range(0, n, item_len)]))
    return out


def test_words_and_paragraphs():
    assert C.words("Haec, nusquam 12 est_Ae!") == ["haec", "nusquam", "est", "ae"]
    assert C.words("Æthiopiæ dell'acqua") == ["æthiopiæ", "dell", "acqua"]
    assert C.paragraphs("a b\nc\n\n\nd\n") == ["a b c", "d"]
    assert C.sentences("One. Two; three: four") == ["One.", "Two;", "three:", "four"]


def test_gutenberg_body_and_read_text(tmp_path):
    assert C.gutenberg_body(GUT).strip() == "body one"
    assert C.gutenberg_body("plain") == "plain"
    p = tmp_path / "latin1.txt"
    p.write_bytes("caf\xe9\r\nx".encode("latin-1"))
    assert C.read_text(p) == "café\nx"
    q = tmp_path / "utf8.txt"
    q.write_bytes("﻿é\r\n".encode("utf-8"))
    assert C.read_text(q) == "é\n"


def test_parse_tess_levels():
    text = "<cels. med. 1.pr> Ut alimenta.\n<cels. med. 1.1> Sanis.\nno tag\n<cels. med. 2.1> Corporibus agri."
    d = C.parse_tess(text)
    assert d == [("1", [["ut", "alimenta"], ["sanis"]]), ("2", [["corporibus", "agri"]])]
    assert [x[0] for x in C.parse_tess(text, level=2)] == ["1.pr", "1.1", "2.1"]


def test_parse_headed():
    text = "front\nLiber V\nA b.\n\nC d e.\nLiber VI\nF."
    assert C.parse_headed(text, r"Liber ([IVX]+)", "{1}") == [("V", [["a", "b"], ["c", "d", "e"]]), ("VI", [["f"]])]


def test_parse_linnaeus():
    d = C.parse_linnaeus(LINN)
    assert [x[0] for x in d] == ["I MONOGYNIA", "I DIGYNIA", "II MONOGYNIA"]
    assert d[0][1] == [["canna"], ["indica"], ["canna", "foliis", "ovatis", "roy", "lugdb"], ["habitat", "in", "india"]]
    flat = [w for _, its in d for it in its for w in it]
    assert "monandria" not in flat and "alpini" not in flat and "preface" not in flat


def test_parse_culpeper():
    d = C.parse_culpeper(CULP, window=1)
    assert [x[0] for x in d] == ["A", "W"]
    assert d[0][1][0] == ["amara", "dulcis"]
    assert ["wild", "betony", "called", "wood", "betony"] in d[1][1]  # two-line heading
    assert "syrups" not in [w for _, its in d for it in its for w in it]
    s = C.parse_culpeper(CULP, window=1, item="sentence")
    assert ["considering", "divers", "shires"] in s[0][1] and ["well"] in s[0][1]
    # running lower median over 3 entries (A, A, W): the lone W at the edge joins the A run
    assert [x[0] for x in C.parse_culpeper(CULP, window=3)] == ["A"]
    with pytest.raises(ValueError):
        C.parse_culpeper(CULP, item="word")


def test_blocks_and_counts():
    b = C.blocks([str(i) for i in range(100)], 4, limit=80, item_len=7)
    assert [x[0] for x in b] == ["b1", "b2", "b3", "b4"]
    assert C.n_tokens(b) == 80
    assert b[1][1][0][0] == "20"


def test_section_quotas_equal_mass():
    sec = C.section_quotas([2] * 10, 5)
    assert sec.tolist() == [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]
    assert C.section_quotas([5, 1, 1, 1, 1, 1], 2).tolist() == [0, 1, 1, 1, 1, 1]


def test_fill_prose_and_items():
    items = [["a", "b", "c", "d"], ["e"], ["f", "g"], ["h", "i", "j"]]
    assert C._fill(items, [3, 3], "prose") == [["a", "b", "c"], ["d", "e", "f"]]
    # items: every line starts a new item; long items are cut, short ones followed by the next
    assert C._fill(items, [3, 3], "items") == [["a", "b", "c"], ["e", "f", "g"]]
    assert C._fill(items, [3, 3, 3, 3], "items") is None
    assert C._fill(items, [20], "prose") is None
    with pytest.raises(ValueError):
        C._fill(items, [1], "rows")


def test_merge_small_and_select():
    divs = _divs([30, 5, 5, 40, 12, 3])
    m = C.merge_small(divs, lambda items: sum(map(len, items)) >= 10)
    assert [x[0] for x in m] == ["d0", "d1+d2", "d3", "d4+d5"]
    assert C.n_tokens(m) == C.n_tokens(divs)
    sel = C.select_divisions(_divs([50, 50, 50, 50, 50, 50]), [5] * 20, 3)
    assert [x[0] for x in sel] == ["d0", "d2", "d5"]
    # too few divisions: k is lowered
    sel = C.select_divisions(_divs([60, 60]), [5] * 20, 4)
    assert len(sel) == 2


def test_layout_and_control_table():
    tpl = _template()
    lens = C.template_lines(tpl)
    divs = _divs([200, 200, 200, 200])
    lines, labels = C.layout_lines(divs[:3], lens, "prose")
    assert [len(x) for x in lines] == lens
    assert labels[0] == "d0" and labels[-1] == "d2"
    assert lines[0][0] == "w0x0"  # each block starts at the beginning of its division
    t = C.control_table(lines, labels, tpl, name="toy")
    assert len(t) == len(tpl)
    assert (t["folio"].to_numpy() == tpl["folio"].to_numpy()).all()
    assert (t["hand"].fillna("?").to_numpy() == tpl["hand"].fillna("?").to_numpy()).all()
    assert t.attrs["glyph_preset"] is None and t.attrs["source"] == "toy"
    assert (t["len_glyphs"] == t["token"].str.len()).all()
    # one section per page (majority of the page's lines)
    assert t.groupby("folio")["illustration"].nunique().max() == 1
    assert set(t["illustration"]) <= {"d0", "d1", "d2"}
    assert vocabulary_summary(t).iloc[0]["n_tokens"] == len(tpl)
    with pytest.raises(ValueError):
        C.layout_lines(_divs([10, 10]), lens, "prose")


def test_layout_sequential_and_build_control():
    tpl = _template()
    lens = C.template_lines(tpl)
    divs = _divs([40, 40, 40, 40, 40, 40, 40, 40])
    lines, labels = C.layout_sequential(divs, lens, "prose", k=4)
    flat = [w for ln in lines for w in ln]
    assert flat == [w for _, its in divs for it in its for w in it][: sum(lens)]
    assert labels[0] == "d0+d1" and len(set(labels)) <= 4
    lines_i, labels_i = C.layout_sequential(divs, lens, "items", k=4)
    assert all(int(ln[0].split("x")[1]) % 5 == 0 for ln in lines_i)  # every line starts an item
    assert [len(x) for x in lines_i] == lens
    t = C.build_control(_divs([300] * 6), tpl, k=3, mode="prose", name="bal", min_k=3)
    assert t.attrs["control"]["strategy"] == "balanced" and len(t.attrs["control"]["sections"]) == 3
    t2 = C.build_control(_divs([150, 150]), tpl, k=7, mode="prose", name="seq", min_k=4)
    assert t2.attrs["control"]["strategy"] == "sequential"
    t3 = C.build_control(_divs([300] * 6), tpl, k=3, mode="items", strategy="sequential")
    assert t3.attrs["control"]["mode"] == "items" and len(t3) == len(tpl)
    with pytest.raises(ValueError):
        C.build_control(divs, tpl, k=3, strategy="random")
    with pytest.raises(ValueError):
        C.layout_sequential(_divs([10]), lens, "prose")
