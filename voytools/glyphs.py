"""Glyph handling: Eva composite collapsing and alphabet defaults.

Eva transliterates the Voynich script with a stroke-level alphabet in which a
few frequent letter groups are widely believed to be single glyphs
("composites"): the bench gallows ``cth ckh cph cfh`` and the bench glyphs
``ch sh``.  Entropy and morphology statistics are sensitive to whether such
groups are counted as one symbol or two (Rozanova & Temerev 2026 collapse them
before measuring glyph entropy; Zattera 2022 uses a similar glyph inventory).

Presets
-------
``'basic'``
    ``cth, ckh, cph, cfh, ch, sh`` (longest match first, scanning left to right).
``'extended'``
    ``basic`` plus ``iin, in, ee``; ``iin`` is tried before ``in``.
``None`` / ``'none'``
    Identity: every character is its own glyph.  This is the right choice for
    non-Eva alphabets such as v101 (Glen Claston), where ``c``/``h``/``s`` are
    unrelated letters.

Implementation note: matching is a single left-to-right scan with a
longest-alternative-first regular expression.  For the Eva presets this gives
the same result as sequential ``str.replace`` in the order composites, ``iin``,
``in``, ``ee``, because no composite shares a letter with ``iin``/``in``/``ee``
other than through ``i``/``n`` positions that cannot overlap.

A token that still contains an uncertain-space marker (``space_policy='strict'``)
is collapsed without it; an alternative reading ``[a:b]`` (``alternatives='mark'``)
is kept together as ONE pseudo-glyph.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Iterable, Sequence

__all__ = [
    "EVA_BASIC",
    "EVA_EXTENDED",
    "PRESETS",
    "is_eva_alphabet",
    "default_preset",
    "resolve_patterns",
    "collapse_eva",
    "glyph_sequence",
    "to_symbol_string",
    "to_symbol_strings",
]

EVA_BASIC: tuple[str, ...] = ("cth", "ckh", "cph", "cfh", "ch", "sh")
EVA_EXTENDED: tuple[str, ...] = EVA_BASIC + ("iin", "in", "ee")
PRESETS: dict[str, tuple[str, ...]] = {
    "basic": EVA_BASIC,
    "extended": EVA_EXTENDED,
    "none": (),
}

_MARK = ","


def is_eva_alphabet(alphabet: str | None) -> bool:
    """True for the Eva alphabet codes ``Eva-`` and ``EvaT`` (case-insensitive)."""
    return bool(alphabet) and alphabet.lower().startswith("eva")


def default_preset(alphabet: str | None) -> str | None:
    """Default collapsing preset for an IVTFF alphabet code: ``'basic'`` for Eva, else ``None``."""
    return "basic" if is_eva_alphabet(alphabet) else None


def resolve_patterns(preset: str | Iterable[str] | None) -> tuple[str, ...]:
    """Return the tuple of composite patterns for a preset name, a custom iterable, or ``None``."""
    if preset is None:
        return ()
    if isinstance(preset, str):
        try:
            return PRESETS[preset.lower()]
        except KeyError:
            raise ValueError(f"unknown glyph preset {preset!r}; use one of {sorted(PRESETS)} or None") from None
    return tuple(preset)


@lru_cache(maxsize=32)
def _scanner(patterns: tuple[str, ...]) -> re.Pattern[str]:
    ordered = sorted(set(patterns), key=lambda p: (-len(p), p))
    alts = [r"\[[^\]]*\]"] + [re.escape(p) for p in ordered] + ["."]
    return re.compile("|".join(alts), re.S)


def collapse_eva(token: str, preset: str | Iterable[str] | None = "basic") -> tuple[str, ...]:
    """Split ``token`` into glyphs, merging Eva composites according to ``preset``.

    Examples
    --------
    >>> collapse_eva("chedy")
    ('ch', 'e', 'd', 'y')
    >>> collapse_eva("qokeedy", "extended")
    ('q', 'o', 'k', 'ee', 'd', 'y')
    >>> collapse_eva("daiin", "extended")
    ('d', 'a', 'iin')
    >>> collapse_eva("cthol", None)
    ('c', 't', 'h', 'o', 'l')
    """
    if _MARK in token:
        token = token.replace(_MARK, "")
    return tuple(_scanner(resolve_patterns(preset)).findall(token))


def glyph_sequence(
    df,
    preset: str | Iterable[str] | None = "basic",
    separator: str | None = None,
    scope: str = "line",
) -> list[list[str]]:
    """Glyph sequences of a token table, one list per line (locus).

    Parameters
    ----------
    df:
        Token table with a ``token`` column (and ``locus_id``/``folio`` so lines can be
        recognised).  Rows must be in text order.
    preset:
        Collapsing preset (see module docs); ``'auto'`` uses the preset stored in
        ``df.attrs['glyph_preset']`` (``'basic'`` if absent).
    separator:
        Optional boundary symbol inserted between consecutive tokens (never at the
        start or end of a sequence).  ``None`` removes spaces, as in Rozanova & Temerev.
    scope:
        ``'line'``: one sequence per locus; ``'all'``: a single concatenated sequence
        (with ``separator`` also between lines).
    """
    from ._common import glyph_lists, line_ids

    if scope not in ("line", "all"):
        raise ValueError("scope must be 'line' or 'all'")
    glyphs = glyph_lists(df, preset)
    if scope == "all":
        out: list[str] = []
        for i, g in enumerate(glyphs):
            if separator is not None and i > 0:
                out.append(separator)
            out.extend(g)
        return [out]
    ids = line_ids(df)
    lines: list[list[str]] = []
    prev = None
    for lid, g in zip(ids.tolist(), glyphs):
        if lid != prev:
            lines.append([])
            prev = lid
        cur = lines[-1]
        if separator is not None and cur:
            cur.append(separator)
        cur.extend(g)
    return lines


def _assign_symbols(glyphs: Iterable[str], mapping: dict[str, str] | None) -> dict[str, str]:
    mapping = dict(mapping or {})
    used = set(mapping.values())
    new = sorted({g for g in glyphs if g not in mapping})
    singles = [g for g in new if len(g) == 1 and g not in used]
    for g in singles:
        mapping[g] = g
        used.add(g)
    nxt = 0xE200  # private-use area, clear of the high-ASCII block 0xE080-0xE0FF
    for g in new:
        if g in mapping:
            continue
        while chr(nxt) in used:
            nxt += 1
        mapping[g] = chr(nxt)
        used.add(chr(nxt))
    return mapping


def to_symbol_string(
    glyphs: Iterable[str], mapping: dict[str, str] | None = None
) -> tuple[str, dict[str, str]]:
    """Map a flat glyph sequence to a string with ONE character per distinct glyph.

    Single-character glyphs map to themselves; multi-character glyphs (``'ch'``,
    ``'iin'`` ...) get a private-use character.  Pass the returned mapping back in to
    encode further sequences consistently.

    Returns
    -------
    (string, mapping)
    """
    seq = list(glyphs)
    mapping = _assign_symbols(seq, mapping)
    return "".join(mapping[g] for g in seq), mapping


def to_symbol_strings(
    sequences: Iterable[Sequence[str]], mapping: dict[str, str] | None = None
) -> tuple[list[str], dict[str, str]]:
    """Like :func:`to_symbol_string` for many sequences with one shared mapping."""
    seqs = [list(s) for s in sequences]
    mapping = _assign_symbols((g for s in seqs for g in s), mapping)
    return ["".join(mapping[g] for g in s) for s in seqs], mapping
