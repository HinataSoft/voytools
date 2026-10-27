"""Token, locus and page tables (pandas) built from a parsed IVTFF file.

:func:`build_tables` turns a :class:`~voytools.ivtff.Transliteration` into three
DataFrames:

``tokens``
    one row per token in file order; see the column dictionary in the README.
``loci``
    one row per locus (text line, label, circular/radial text).
``pages``
    one row per page with its page variables and counts.

Conventions
-----------
* ``folio`` holds the IVTFF *page name* (``f1r``, ``f67r1``, ``fRos``).
* Page variables (``$Q`` quire, ``$P`` page in quire, ``$F`` folio in quire, ``$B``
  bifolio, ``$I`` illustration, ``$L`` Currier language, ``$H`` hand, ``$C`` Currier's
  hand, ``$X`` extraneous writing) become string columns; a missing variable (or the
  value ``@`` before the first text tag sets it) is NA.  Text tags ``<@H=2>`` are
  applied to the whole line they stand on and to the following loci of the page.
* Paragraphs come from ``<%>`` / ``<$>``.  ``par_index`` counts paragraphs on the
  page (0-based) and is NA outside paragraphs and for non-``P`` loci.  Files without
  marks (e.g. RF1b) can use ``infer_paragraphs=True``: a ``P`` locus outside an
  explicit paragraph starts an implicit one when its locator is ``@`` or ``*``.
* In interlinear files every transcriber has its own line numbering; select one
  transcriber (``transcriber='H'``) before running statistics.
* The tables carry ``df.attrs['alphabet']`` and ``df.attrs['glyph_preset']``; the
  glyph-level statistics use them for ``preset='auto'``.
"""

from __future__ import annotations

import warnings as _warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .glyphs import collapse_eva, default_preset
from .ivtff import (
    UNCERTAIN_SPACE_MARK,
    Header,
    Locus,
    ParseWarning,
    Transliteration,
    parse_ivtff,
    render_high_ascii,
    tokenize_locus,
)

__all__ = [
    "TOKEN_COLUMNS",
    "Tables",
    "build_tables",
    "load_tables",
    "filter_tokens",
    "export",
    "most_complete_transcriber",
]

TOKEN_COLUMNS = [
    "idx", "folio", "quire", "page_in_quire", "folio_in_quire", "bifolio", "illustration",
    "language", "hand", "currier_hand", "extraneous", "locus_id", "locus_num", "locator",
    "locus_type", "generic_type", "subtype", "transcriber", "par_index", "line_in_par",
    "is_par_first_line", "is_par_last_line", "line_index", "token_index", "n_tokens",
    "pos_from_end", "rel_pos", "is_first", "is_last", "token", "token_raw", "len_chars",
    "len_glyphs", "has_unreadable", "has_alternative", "has_ligature", "has_high_ascii",
    "after_uncertain_space", "before_uncertain_space", "after_drawing_break",
    "before_drawing_break",
]

_PAGE_VARS = [
    ("quire", "Q"),
    ("page_in_quire", "P"),
    ("folio_in_quire", "F"),
    ("bifolio", "B"),
    ("illustration", "I"),
    ("language", "L"),
    ("hand", "H"),
    ("currier_hand", "C"),
    ("extraneous", "X"),
]

_DEFAULT_CLEAN = ("unreadable", "alternative", "high_ascii")
_CLEAN_FLAGS = {
    "unreadable": "has_unreadable",
    "alternative": "has_alternative",
    "high_ascii": "has_high_ascii",
    "ligature": "has_ligature",
}


@dataclass
class Tables:
    """The three tables plus provenance."""

    tokens: pd.DataFrame
    loci: pd.DataFrame
    pages: pd.DataFrame
    header: Header
    options: dict[str, Any] = field(default_factory=dict)
    warnings: list[ParseWarning] = field(default_factory=list)
    transcribers: list[str] = field(default_factory=list)

    @property
    def alphabet(self) -> str | None:
        """Alphabet code of the source file."""
        return self.header.alphabet

    def main_transcriber(self) -> str | None:
        """Transcriber with the most loci (``None`` for non-interlinear files)."""
        return most_complete_transcriber(self.loci)

    def export(self, directory: str | Path, fmt: str = "csv") -> list[Path]:
        """Write ``tokens``, ``loci`` and ``pages`` into ``directory``; returns the paths."""
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        ext = "parquet" if fmt == "parquet" else "csv"
        return [export(getattr(self, n), d / f"{n}.{ext}", fmt) for n in ("tokens", "loci", "pages")]


# --------------------------------------------------------------------------- #
# building
# --------------------------------------------------------------------------- #


class _ParState:
    __slots__ = ("n_lines", "par", "lip", "open", "explicit")

    def __init__(self) -> None:
        self.n_lines = 0
        self.par = -1
        self.lip = -1
        self.open = False
        self.explicit = False


def most_complete_transcriber(loci: pd.DataFrame | Transliteration) -> str | None:
    """Return the transcriber ID with the most loci (``None`` if there are none)."""
    if isinstance(loci, Transliteration):
        ids = [lc.transcriber for lc in loci.loci if lc.transcriber]
    else:
        ids = [t for t in loci["transcriber"].tolist() if isinstance(t, str)]
    if not ids:
        return None
    counts: dict[str, int] = {}
    for t in ids:
        counts[t] = counts.get(t, 0) + 1
    return sorted(counts, key=lambda t: (-counts[t], t))[0]


def _var(lc: Locus, key: str) -> str | None:
    v = lc.variables.get(key)
    return None if v in (None, "@") else v


def build_tables(
    translit: Transliteration,
    *,
    transcriber: str | None = None,
    glyph_preset: str | None = "auto",
    infer_paragraphs: bool = False,
    **tokenize_options: Any,
) -> Tables:
    """Build the token, locus and page tables.

    Parameters
    ----------
    translit:
        Result of :func:`~voytools.ivtff.parse_ivtff`.
    transcriber:
        Keep only loci of this transcriber ID (interlinear files).  ``None`` keeps all;
        then every transcriber's lines appear (with separate ``line_index`` numbering).
    glyph_preset:
        Preset used for ``len_glyphs`` (``'basic'``, ``'extended'``, ``None``).
        ``'auto'`` = ``'basic'`` for Eva alphabets and ``None`` (no collapsing) otherwise.
    infer_paragraphs:
        Infer paragraph structure for ``P`` loci that lie outside ``<%>...<$>`` marks.
    **tokenize_options:
        ``space_policy`` (``'split'``, ``'join'``, ``'strict'``) and ``alternatives``
        (``'first'``, ``'mark'``); see :func:`voytools.ivtff.tokenize_text`.
    """
    space_policy = tokenize_options.get("space_policy", "split")
    alternatives = tokenize_options.get("alternatives", "first")
    unknown = set(tokenize_options) - {"space_policy", "alternatives"}
    if unknown:
        raise TypeError(f"unknown tokenizer option(s): {sorted(unknown)}")
    alphabet = translit.header.alphabet
    preset = default_preset(alphabet) if glyph_preset == "auto" else glyph_preset

    loci_sel = [lc for lc in translit.loci if transcriber is None or lc.transcriber == transcriber]
    states: dict[tuple[int, str | None], _ParState] = {}
    glyph_cache: dict[str, int] = {}

    def n_glyphs(tok: str) -> int:
        v = glyph_cache.get(tok)
        if v is None:
            v = len(collapse_eva(tok, preset))
            glyph_cache[tok] = v
        return v

    # ---- pass 1: loci -------------------------------------------------------
    locus_rows: list[dict[str, Any]] = []
    locus_tokens: list[list] = []
    for lc in loci_sel:
        key = (lc.page_index, lc.transcriber)
        st = states.get(key)
        if st is None:
            st = states[key] = _ParState()
        line_index = st.n_lines
        st.n_lines += 1
        par_index: int | None = None
        line_in_par: int | None = None
        if lc.generic_type == "P":
            if lc.paragraph_start:
                st.par += 1
                st.lip = 0
                st.open = True
                st.explicit = True
            elif st.open and st.explicit:
                st.lip += 1
            elif infer_paragraphs:
                if st.open and lc.locator not in "@*":
                    st.lip += 1
                else:
                    st.par += 1
                    st.lip = 0
                    st.open = True
                    st.explicit = False
            if st.open:
                par_index, line_in_par = st.par, st.lip
            if lc.paragraph_end:
                st.open = False
        toks = tokenize_locus(lc, space_policy, alternatives)
        row: dict[str, Any] = {
            "folio": lc.folio,
            **{name: _var(lc, code) for name, code in _PAGE_VARS},
            "locus_id": lc.id,
            "locus_num": lc.num,
            "num_suffix": lc.num_suffix,
            "locator": lc.locator,
            "locus_type": lc.locus_type,
            "generic_type": lc.generic_type,
            "subtype": lc.subtype,
            "transcriber": lc.transcriber,
            "line_index": line_index,
            "par_index": par_index,
            "line_in_par": line_in_par,
            "is_par_first_line": line_in_par == 0 if line_in_par is not None else False,
            "is_par_last_line": False,
            "n_tokens": len(toks),
            "n_chars": sum(len(t.text.replace(UNCERTAIN_SPACE_MARK, "")) for t in toks),
            "paragraph_start": lc.paragraph_start,
            "paragraph_end": lc.paragraph_end,
            "comments": " | ".join(lc.comments),
            "inline_tags": ";".join(f"{k}={v}" for k, v in lc.inline_tags),
            "raw_text": lc.raw_text,
            "text": lc.text,
            "line_no": lc.line_no,
            "continued": lc.continued,
        }
        locus_rows.append(row)
        locus_tokens.append(toks)

    # last lines of paragraphs: explicit <$> mark, or the next P locus starts another paragraph
    next_par: dict[tuple[int, str | None], tuple[int, Any]] = {}
    for i in range(len(loci_sel) - 1, -1, -1):
        lc, row = loci_sel[i], locus_rows[i]
        if lc.generic_type != "P" or row["par_index"] is None:
            continue
        key = (lc.page_index, lc.transcriber)
        nxt = next_par.get(key)
        row["is_par_last_line"] = bool(
            lc.paragraph_end or nxt is None or nxt[1] != row["par_index"]
        )
        next_par[key] = (i, row["par_index"])

    # ---- pass 2: tokens -------------------------------------------------------
    cols: dict[str, list] = {c: [] for c in TOKEN_COLUMNS}
    idx = 0
    for row, toks in zip(locus_rows, locus_tokens):
        n = len(toks)
        for t in toks:
            ti = t.index
            cols["idx"].append(idx)
            idx += 1
            for c in ("folio", "quire", "page_in_quire", "folio_in_quire", "bifolio", "illustration",
                      "language", "hand", "currier_hand", "extraneous", "locus_id", "locus_num",
                      "locator", "locus_type", "generic_type", "subtype", "transcriber",
                      "par_index", "line_in_par", "is_par_first_line", "is_par_last_line",
                      "line_index"):
                cols[c].append(row[c])
            cols["token_index"].append(ti)
            cols["n_tokens"].append(n)
            cols["pos_from_end"].append(n - 1 - ti)
            cols["rel_pos"].append(ti / (n - 1) if n > 1 else 0.0)
            cols["is_first"].append(ti == 0)
            cols["is_last"].append(ti == n - 1)
            cols["token"].append(t.text)
            cols["token_raw"].append(t.raw)
            cols["len_chars"].append(len(t.text) - t.n_uncertain_inside)
            cols["len_glyphs"].append(n_glyphs(t.text))
            cols["has_unreadable"].append(t.has_unreadable)
            cols["has_alternative"].append(t.has_alternative)
            cols["has_ligature"].append(t.has_ligature)
            cols["has_high_ascii"].append(t.has_high_ascii)
            cols["after_uncertain_space"].append(t.after_uncertain_space)
            cols["before_uncertain_space"].append(t.before_uncertain_space)
            cols["after_drawing_break"].append(t.after_drawing_break)
            cols["before_drawing_break"].append(t.before_drawing_break)
    tokens = _frame(cols, TOKEN_COLUMNS)

    loci_cols = [
        "folio", *[n for n, _ in _PAGE_VARS], "locus_id", "locus_num", "num_suffix", "locator",
        "locus_type", "generic_type", "subtype", "transcriber", "line_index", "par_index",
        "line_in_par", "is_par_first_line", "is_par_last_line", "n_tokens", "n_chars",
        "paragraph_start", "paragraph_end", "comments", "inline_tags", "raw_text", "text",
        "line_no", "continued",
    ]
    loci_df = _frame({c: [r[c] for r in locus_rows] for c in loci_cols}, loci_cols)

    # ---- pages ------------------------------------------------------------------
    counts_loci: dict[int, int] = {}
    counts_p: dict[int, int] = {}
    counts_tok: dict[int, int] = {}
    for lc, row in zip(loci_sel, locus_rows):
        counts_loci[lc.page_index] = counts_loci.get(lc.page_index, 0) + 1
        counts_tok[lc.page_index] = counts_tok.get(lc.page_index, 0) + row["n_tokens"]
        if lc.generic_type == "P":
            counts_p[lc.page_index] = counts_p.get(lc.page_index, 0) + 1
    page_cols = ["folio", "page_index", "line_no", *[n for n, _ in _PAGE_VARS],
                 "n_loci", "n_loci_P", "n_tokens", "n_comment_lines"]
    page_data: dict[str, list] = {c: [] for c in page_cols}
    for p in translit.pages:
        if transcriber is not None and p.index not in counts_loci:
            continue
        page_data["folio"].append(p.name)
        page_data["page_index"].append(p.index)
        page_data["line_no"].append(p.line_no)
        for name, code in _PAGE_VARS:
            page_data[name].append(p.variables.get(code))
        page_data["n_loci"].append(counts_loci.get(p.index, 0))
        page_data["n_loci_P"].append(counts_p.get(p.index, 0))
        page_data["n_tokens"].append(counts_tok.get(p.index, 0))
        page_data["n_comment_lines"].append(len(p.comments))
    pages_df = _frame(page_data, page_cols)

    attrs = {
        "alphabet": alphabet,
        "glyph_preset": preset,
        "space_policy": space_policy,
        "alternatives": alternatives,
        "source": translit.path,
    }
    for df in (tokens, loci_df, pages_df):
        df.attrs.update(attrs)
    return Tables(
        tokens=tokens,
        loci=loci_df,
        pages=pages_df,
        header=translit.header,
        options={**attrs, "transcriber": transcriber, "infer_paragraphs": infer_paragraphs},
        warnings=list(translit.warnings),
        transcribers=translit.transcribers,
    )


def _frame(cols: dict[str, list], order: list[str]) -> pd.DataFrame:
    """DataFrame from column lists with explicit dtypes for the integer-with-NA columns."""
    data: dict[str, Any] = {}
    for c in order:
        values = cols[c]
        if c in ("par_index", "line_in_par"):
            data[c] = pd.array(values, dtype="Int64")
        else:
            data[c] = values
    df = pd.DataFrame(data, columns=order)
    for c in order:
        if c.startswith(("is_", "has_", "after_", "before_")) or c in ("paragraph_start", "paragraph_end", "continued"):
            df[c] = df[c].astype(bool)
    return df


def load_tables(path: str | Path, *, encoding: str = "utf-8", errors: str = "replace", **kwargs: Any) -> Tables:
    """Parse ``path`` and build the tables (``kwargs`` go to :func:`build_tables`)."""
    return build_tables(parse_ivtff(path, encoding=encoding, errors=errors), **kwargs)


# --------------------------------------------------------------------------- #
# filtering and export
# --------------------------------------------------------------------------- #


def _as_list(value: Any) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        if value.strip().lower() == "all":
            return None
        return [v.strip() for v in value.split(",") if v.strip()]
    if isinstance(value, (int, np.integer)):
        return [str(int(value))]
    return [str(v) for v in value]


def filter_tokens(
    df: pd.DataFrame,
    generic_type: str | Iterable[str] | None = None,
    language: str | Iterable[str] | None = None,
    illustration: str | Iterable[str] | None = None,
    hand: str | int | Iterable[str | int] | None = None,
    *,
    transcriber: str | Iterable[str] | None = None,
    quire: str | Iterable[str] | None = None,
    locus_type: str | Iterable[str] | None = None,
    subtype: str | Iterable[str] | None = None,
    folio: str | Iterable[str] | None = None,
    clean: bool | Iterable[str] = True,
) -> pd.DataFrame:
    """Select tokens by metadata and (optionally) drop poorly readable tokens.

    Every selector takes one value, a list, or a comma-separated string
    (``'P,L'``); ``'all'`` or ``None`` means no restriction.  Rows with a missing
    value in a selected column never match.

    Parameters
    ----------
    clean:
        ``True`` (default) drops tokens flagged ``has_unreadable``, ``has_alternative`` or
        ``has_high_ascii``.  ``False`` keeps everything.  An iterable chooses flags
        explicitly from ``'unreadable'``, ``'alternative'``, ``'high_ascii'``,
        ``'ligature'``.

    Returns
    -------
    DataFrame
        A copy with a fresh 0..n-1 index (``idx`` keeps the original token number)
        and ``attrs`` preserved.
    """
    mask = np.ones(len(df), dtype=bool)
    for col, val in (
        ("generic_type", generic_type), ("language", language), ("illustration", illustration),
        ("hand", hand), ("transcriber", transcriber), ("quire", quire),
        ("locus_type", locus_type), ("subtype", subtype), ("folio", folio),
    ):
        vals = _as_list(val)
        if vals is None:
            continue
        mask &= df[col].isin(vals).to_numpy(dtype=bool)
    if clean is True:
        flags = _DEFAULT_CLEAN
    elif clean is False or clean is None:
        flags = ()
    else:
        flags = tuple(clean)
        bad = set(flags) - set(_CLEAN_FLAGS)
        if bad:
            raise ValueError(f"unknown clean flag(s) {sorted(bad)}; use {sorted(_CLEAN_FLAGS)}")
    for f in flags:
        mask &= ~df[_CLEAN_FLAGS[f]].to_numpy(dtype=bool)
    out = df.loc[mask].reset_index(drop=True)
    out.attrs.update(df.attrs)
    return out


_RENDER_COLUMNS = ("token", "prefix", "stem", "suffix")


def export(df: pd.DataFrame, path: str | Path, fmt: str = "auto", *, render_high_ascii_codes: bool = True) -> Path:
    """Write a table to CSV or Parquet and return the path actually written.

    ``fmt='auto'`` chooses by suffix (``.parquet`` or anything else = CSV).  Parquet needs
    ``pyarrow`` (or ``fastparquet``); without it a warning is issued and a CSV file with
    the same stem is written instead.  In CSV the private-use symbols for high-ASCII
    characters are rendered back as ``@nnn;`` unless ``render_high_ascii_codes=False``.
    """
    p = Path(path)
    if fmt == "auto":
        fmt = "parquet" if p.suffix.lower() in (".parquet", ".pq") else "csv"
    if fmt not in ("csv", "parquet"):
        raise ValueError("fmt must be 'auto', 'csv' or 'parquet'")
    p.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "parquet":
        try:
            df.to_parquet(p, index=False)
            return p
        except ImportError:
            p = p.with_suffix(".csv")
            _warnings.warn(
                "pyarrow is not installed; writing CSV instead of Parquet: " + str(p), RuntimeWarning, stacklevel=2
            )
    out = df
    if render_high_ascii_codes:
        out = df.copy()
        for c in _RENDER_COLUMNS:
            if c in out.columns:
                out[c] = [render_high_ascii(v) if isinstance(v, str) else v for v in out[c].tolist()]
    out.to_csv(p, index=False, encoding="utf-8")
    return p
