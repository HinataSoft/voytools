"""Port of the Timm & Schinner self-citation generator (voytools.selfcitation)."""

from __future__ import annotations

import hashlib
import re

import pytest

from voytools import filter_tokens, load_tables
from voytools import nulls as N
from voytools import selfcitation as SC
from voytools._common import line_ids
from voytools.glyphs import collapse_eva
from tests.conftest import make_tokens
from tests.test_wordgrammar import SEED, _digest

EVA_LETTERS = set("abcdefghijklmnopqrstuvxyz")


def _template(n_pages=6, lines=12, n=9):
    rows = []
    for p in range(n_pages):
        for i in range(lines):
            rows.append((f"g{p}", n - (i % 4), i % 6 == 0, i % 6 == 5))
    return rows


def _sha(lines) -> str:
    return hashlib.sha256("\n".join(" ".join(line) for line in lines).encode()).hexdigest()[:16]


def test_generate_is_deterministic_and_seed_dependent():
    t = _template()
    a = SC.selfcitation_generate(t, seed=3)
    assert a == SC.selfcitation_generate(t, seed=3)
    assert a != SC.selfcitation_generate(t, seed=4)
    b = SC.selfcitation_generate(t, seed=3, p_copy=0.2, replace_weights=(0, 0, 100))
    assert b == SC.selfcitation_generate(t, seed=3, p_copy=0.2, replace_weights=(0, 0, 100))
    assert b != a


def test_generate_matches_template_line_counts():
    t = _template()
    out = SC.selfcitation_generate(t, seed=1)
    assert [len(line) for line in out] == [r[1] for r in t]
    # 2- and 3-element rows are accepted
    out2 = SC.selfcitation_generate([(r[1], r[2]) for r in t], seed=1)
    assert [len(line) for line in out2] == [r[1] for r in t]
    with pytest.raises(ValueError):
        SC.selfcitation_generate([(1,)], seed=0)


def test_generate_fingerprint_is_stable():
    """Regression: recorded when the module was ported (identical to the R&T class output)."""
    t = _template()
    assert _sha(SC.selfcitation_generate(t, seed=1)) == "88b8b7d8a4fd1cf6"
    assert _sha(SC.selfcitation_generate(t, seed=2, replace_weights=(0, 0, 100))) == "9ee8511835892127"


def test_output_is_basic_eva_on_2000_tokens():
    t = [(f"g{i // 20}", 10, i % 8 == 0) for i in range(200)]  # 2 000 tokens
    out = SC.selfcitation_generate(t, seed=7)
    toks = [w for line in out for w in line]
    assert len(toks) == 2000
    assert all(re.fullmatch(r"[a-z]+", w) for w in toks)
    assert set("".join(toks)) <= EVA_LETTERS
    for w in set(toks):
        assert "".join(collapse_eva(w, "extended")) == w


def test_tokenize_longest_ligature_match():
    assert SC.tokenize("qokeedy") == ["qo", "k", "ee", "dy"]
    assert SC.tokenize("daiin") == ["d", "a", "iin"]
    assert SC.G("chol").toks == ["ch", "ol"]


def test_template_and_like_keep_layout():
    df = make_tokens([["qokeedy", "chedy", "ol"], ["daiin", "shedy"], ["okaiin", "dar", "s", "or"]] * 5, page_size=4)
    t = SC.selfcitation_template(df)
    assert [r[1] for r in t] == [3, 2, 4] * 5
    s = SC.selfcitation_like(df, seed=2)
    assert len(s) == len(df)
    assert (line_ids(s) == line_ids(df)).all()
    for c in ("folio", "locus_id", "token_index", "n_tokens", "is_par_first_line"):
        assert s[c].tolist() == df[c].tolist()
    assert s.attrs["source"] == "selfcitation"
    assert s["len_chars"].tolist() == [len(w) for w in s["token"]]
    assert s.equals(SC.selfcitation_like(df, seed=2))
    with pytest.raises(TypeError):
        SC.selfcitation_like(df)


def test_selfcitation_like_on_real_file(zl_path):
    zl = filter_tokens(load_tables(zl_path, space_policy="split").tokens, generic_type="P", clean=True)
    sub = zl[zl["folio"].isin(zl["folio"].drop_duplicates().head(8))].reset_index(drop=True)
    s = SC.selfcitation_like(sub, seed=1)
    assert len(s) == len(sub)
    assert s["locus_id"].tolist() == sub["locus_id"].tolist()
    assert set("".join(s["token"])) <= EVA_LETTERS


def test_autocopy_generate_is_untouched():
    """The port must not change the existing auto-copy generator (fingerprints from test_wordgrammar)."""
    assert _digest(N.autocopy_generate(SEED, 2000, seed=11)) == "a14bc553fde69ca0"


# --------------------------------------------------------------------------- #
# experiment 018: memory across pages (p_prev_page / page_memory)
# --------------------------------------------------------------------------- #

_PAGED = [(f"p{p}", n, i == 0, False) for p in range(6) for i, n in enumerate([7, 9, 8, 6, 10])]


def test_memory_defaults_are_bit_identical():
    """``p_prev_page=0`` must not touch the random stream: output identical to the unextended generator."""
    for seed in (0, 1, 7):
        base = SC.selfcitation_generate(_PAGED, seed=seed)
        same = SC.selfcitation_generate(_PAGED, seed=seed, p_prev_page=0.0, page_memory=8)
        assert same == base


def test_memory_changes_output_and_keeps_layout():
    base = SC.selfcitation_generate(_PAGED, seed=3)
    mem = SC.selfcitation_generate(_PAGED, seed=3, p_prev_page=0.5, page_memory=2)
    assert [len(l) for l in mem] == [len(l) for l in base]
    assert mem != base
    assert mem == SC.selfcitation_generate(_PAGED, seed=3, p_prev_page=0.5, page_memory=2)
    # the first page has no earlier page: identical to the base run
    assert mem[:5] == base[:5]


def test_memory_parameter_validation():

    with pytest.raises(ValueError):
        SC.selfcitation_generate(_PAGED, seed=0, p_prev_page=1.5)
    with pytest.raises(ValueError):
        SC.selfcitation_generate(_PAGED, seed=0, p_prev_page=0.2, page_memory=0.5)


def test_page_offset_is_geometric_and_clamped():
    from voytools.selfcitation import SelfCitation

    gen = SelfCitation(seed=0, p_prev_page=0.5, page_memory=4)
    gen.prev_pages = [[]] * 3
    ks = [gen.page_offset() for _ in range(2000)]
    assert min(ks) == 1 and max(ks) == 3
    gen.prev_pages = [[]] * 1000
    ks = [gen.page_offset() for _ in range(4000)]
    assert 3.5 < sum(ks) / len(ks) < 4.5
    gen1 = SelfCitation(seed=0, p_prev_page=0.5, page_memory=1)
    gen1.prev_pages = [[]] * 5
    assert {gen1.page_offset() for _ in range(50)} == {1}
