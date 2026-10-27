"""Shared fixtures: a synthetic mini-IVTFF file and locations of the real data files."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from voytools.ivtff import parse_ivtff_text
from voytools.table import build_tables

# --------------------------------------------------------------------------- #
# synthetic IVTFF 2.0 file
# --------------------------------------------------------------------------- #

SYNTH = """\
#=IVTFF Eva- 2.0 M 1
# synthetic test file
#
<f1r>      <! $Q=A $P=A $F=a $B=1 $I=T $L=A $H=1 $C=1 $X=V>
# first page
<f1r.1,@P0>       <%>fachys.ykal.ar<!strange r>.[cth:oto]res.y,kor.sholdy
<f1r.2,+P0>       sory.{ck}hy.@150;ar.o?x.???.d<->kor.a<~>b
<f1r.3,*P0>       <@L=B>qokeedy.dai[n:r]<!x>.ol,/
/ chedy.shy
<f1r.4,=Pt>       ydaraishy<$>
<f1r.5,*P0>       <%>solo.line<$>
<f1r.6,&Cc>       oro.r,
<f1r.7,~Lz>       zod
<f1r.8,/Ls>       <!label>sta[r:s]
<f1r.9,!L0>       none
<f1v>      <! $Q=A $P=B $L=@ $H=1 $I=H>
<f1v.1,@Ri>       <!10:00>okeo,dy
<f1v.2,+Ro>       <@L=A>ko.ro
<f1v.3,+Lp>       plant
"""

SYNTH_INTERLINEAR = """\
#=IVTFF Eva- 1.5
# interlinear, old style fillers
<f1r>      <! $Q=A $L=A>
<f1r.1,@P0;H>      fachys.ykal.ar!.cth%res
<f1r.1,@P0;C>      fachys.ykal.ar.cthorys.
<f1r.2,+P0;H>      sory .ckhar. o!r
<f1r.2,+P0;C>      sory.ckhar.o.r
<f1r.2,+P0;U>      %%%%
<f1r.2,+P0;Z>      sory.ckhar.or<$>
<f2r>      <! $Q=A $L=B>
<f2r.44a,+P0;H>    odd.number
"""

SYNTH_V101 = """\
#=IVTFF v101 2.0 M 1
<f1r>      <! $Q=A $L=A $H=1>
<f1r.1,@P0>       <%>fa!s.%a9,hae.ch9<$>
"""


@pytest.fixture(scope="session")
def synth_text() -> str:
    return SYNTH


@pytest.fixture(scope="session")
def synth_path(tmp_path_factory) -> Path:
    p = tmp_path_factory.mktemp("synth") / "synth.txt"
    p.write_text(SYNTH, encoding="utf-8")
    return p


@pytest.fixture(scope="session")
def synth(synth_text):
    return parse_ivtff_text(synth_text)


@pytest.fixture(scope="session")
def synth_tables(synth):
    return build_tables(synth)


# --------------------------------------------------------------------------- #
# real data (tests are skipped when the files are not found)
# --------------------------------------------------------------------------- #

_CANDIDATES = [
    os.environ.get("VOYTOOLS_DATA"),
    str(Path(__file__).resolve().parents[1] / "data"),  # copy committed to the repository
    str(Path(__file__).resolve().parents[2] / "transliterations"),  # D:\Work\Voynich\transliterations
    str(Path(__file__).resolve().parents[1] / "transliterations"),
    "/mnt/user-data/uploads/Voynich/transliterations",
]


def data_dir() -> Path | None:
    for c in _CANDIDATES:
        if c and (Path(c) / "ZL3b-n.txt").exists():
            return Path(c)
    return None


def _need(name: str) -> Path:
    d = data_dir()
    if d is None or not (d / name).exists():
        pytest.skip(f"{name} not found (set VOYTOOLS_DATA to the transliterations folder)")
    return d / name


@pytest.fixture(scope="session")
def zl_path() -> Path:
    return _need("ZL3b-n.txt")


@pytest.fixture(scope="session")
def gc_path() -> Path:
    return _need("GC2a-n.txt")


@pytest.fixture(scope="session")
def it_path() -> Path:
    return _need("IT2a-n.txt")


@pytest.fixture(scope="session")
def lsi_path() -> Path:
    return _need("LSI_ivtff_0d.txt")


@pytest.fixture(scope="session")
def zl_tables(zl_path):
    from voytools.table import load_tables

    return load_tables(zl_path)


# --------------------------------------------------------------------------- #
# helpers for building small token tables
# --------------------------------------------------------------------------- #


def make_tokens(lines, *, page_size=None, header="#=IVTFF Eva- 2.0 M 1", locus_type="+P0", **kw):
    """Token table from a list of lines, each line a list of tokens (joined with '.').

    ``page_size`` starts a new page after that many lines (default: one page).
    """
    out = [header]
    page = 0
    n = 0
    for i, toks in enumerate(lines):
        if i == 0 or (page_size and i % page_size == 0):
            page += 1
            n = 0
            out.append(f"<f{page}r>      <! $Q=A $L=A>")
        n += 1
        loc = "@P0" if n == 1 else locus_type
        out.append(f"<f{page}r.{n},{loc}>   " + ".".join(toks))
    tables = build_tables(parse_ivtff_text("\n".join(out) + "\n"), **kw)
    return tables.tokens
