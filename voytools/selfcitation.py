"""Self-citation (auto-copying) text generator of Timm & Schinner (2020), ported.

Provenance and licence
----------------------
The mechanism is Timm & Schinner's: T. Timm, A. Schinner, "A possible generating algorithm of the
Voynich manuscript", Cryptologia 44(1):1-19 (2020); Java reference implementation
github.com/TorstenTimm/SelfCitationTextgenerator (MIT).  The code below is a port of the Python
re-implementation by L. Rozanova & D. Temerev (2026, arXiv:2608.17096), file
``analysis/reproduce_selfcitation_control.py`` of github.com/lrozanova/voynich-units (MIT licence,
Copyright (c) 2026 Liudmila Rozanova; port of commit 956a7c4).  Taken over unchanged in logic:
the data tables (``LIGATURES``, the glyph similarity tables ``SUBST`` / ``FINAL_SUBST``, prefix,
line-final and in-group replacement tables, the curve/line shape grammar ``can_before`` /
``can_after``), the token class :class:`G` and the generator class :class:`SelfCitation`.  Left out:
the driver parts (entropy, BPE curve, token-order MI, Levenshtein similarity, decipherment attack,
calibration) that voytools already has in :mod:`voytools.stats` / :mod:`voytools.morphology`.

MIT licence of the original code: Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
to whom the Software is furnished to do so, subject to the following conditions: The above
copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

Mechanism (summary of the original docstring)
---------------------------------------------
The scribe repeatedly copies a token already written nearby (same line to the left, or one of
the lines above on the same page, preferring the same writing position) and modifies it by
(I) replacing glyphs with similarly shaped ones, (II) adding or (III) deleting glyphs,
(IV) combining two tokens, (V) splitting one, (VI) duplicating a (sub)group, or (VII) copying it
unchanged; paragraph-initial lines carry extra gallows (VIII).  Lines are generated to a token
count template (page, number of tokens, paragraph-initial flag), so the output has the line and
page layout of the template.

Output alphabet: the generator works on Eva strings parsed into Timm's ligature units; a
generated token is the concatenation of its units, i.e. a plain basic-Eva string (letters
``a``-``z`` only), so :func:`voytools.glyphs.collapse_eva` treats it exactly like a real token.
The tables are Eva-specific: there is no v101 / GC version.

Public API: :func:`selfcitation_generate`, :func:`selfcitation_template`, :func:`selfcitation_like`
(and the ported :class:`SelfCitation` for direct use).
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Any, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "LIGATURES",
    "SUBST",
    "FINAL_SUBST",
    "DEFAULT_SEED_LINE",
    "G",
    "SelfCitation",
    "can_before",
    "can_after",
    "tokenize",
    "selfcitation_generate",
    "selfcitation_template",
    "selfcitation_like",
]

# ---------------------------------------------------------------------------
# Data tables transcribed from the Timm & Schinner reference implementation
# (Glyph.java, CurveLineCanFollow.java, SlimGroupMorpher.java).  These are the
# DOCUMENTED parts of the mechanism.
# ---------------------------------------------------------------------------

# multi-glyph units ("ligatures") that the generator treats as one token
LIGATURES = sorted(
    """ol or al ar dy qo ch sh cs eee ee cth cthh ckh ckhh cph cfh ith ikh iph
    ifh eke ete in iin iiin ir iir iiir is iis iiis il iil iiil im iim iiim om
    am og ag""".split(),
    key=len,
    reverse=True,
)

GALLOWS = ("k", "t", "p", "f")
LINE_INITIAL = ("o", "y", "d", "s")
COMBINABLE = ("ol", "or", "al", "ar")          # Glyph.combinableLigature
ALL_COMBINABLE = ("ol", "or", "al", "ar", "om", "am")
LIGATURE_TAIL_DROP = {"ol": "lrs", "al": "lrs", "or": "lrs", "ar": "lrs",
                      "om": "lrs", "am": "lrs"}

# similarity table: token -> [(replacement tokens, cumulative percent)]
_S = {
    "k": [("t", 77), ("p", 94), ("f", 100)],
    "t": [("k", 84), ("p", 96), ("f", 100)],
    "p": [("k", 59), ("t", 97), ("f", 100)],
    "f": [("k", 56), ("t", 92), ("p", 100)],
    "in": [("n", 8), ("iin", 84), ("iiin", 87), ("ir", 97), ("iir", 98),
           ("iis", 99), ("il", 100)],
    "iin": [("n", 13), ("in", 70), ("iiin", 74), ("ir", 92), ("iir", 97),
            ("is", 99), ("il", 100)],
    "iiin": [("n", 6), ("in", 31), ("iin", 90), ("ir", 98), ("iir", 100)],
    "ir": [("r", 6), ("in", 33), ("iin", 95), ("iiin", 97), ("iir", 99),
           ("iiir", 100)],
    "iir": [("r", 6), ("in", 30), ("iin", 89), ("iiin", 91), ("ir", 99),
            ("iiir", 100)],
    "iiir": [("ir", 90), ("iir", 100)],
    "is": [("in", 50), ("iis", 100)],
    "iis": [("iin", 50), ("is", 100)],
    "il": [("in", 50), ("iil", 100)],
    "iil": [("iin", 50), ("il", 99), ("iiil", 100)],
    "iiil": [("il", 66), ("iil", 100)],
    "im": [("in", 30), ("iin", 100)],
    "iim": [("in", 30), ("iin", 100)],
    "iiim": [("in", 30), ("iin", 100)],
    "om": [("ol", 45), ("or", 94), ("og", 98), ("omg", 100)],
    "am": [("al", 45), ("ar", 94), ("ag", 98), ("amg", 100)],
    "og": [("or", 30), ("al", 63), ("ar", 100)],
    "ag": [("ol", 48), ("or", 72), ("ar", 100)],
    "ol": [("or", 30), ("al", 63), ("ar", 100)],
    "or": [("ol", 47), ("al", 73), ("ar", 100)],
    "al": [("ol", 48), ("or", 72), ("ar", 100)],
    "ar": [("ol", 49), ("or", 73), ("al", 100)],
    "e": [("e", 50), ("ee", 99), ("eee", 100)],
    "ee": [("ch", 40), (("ch", "e"), 50), ("e", 98), ("eee", 100)],
    "eee": [("ch", 10), (("ch", "e"), 20), ("ee", 65), ("e", 100)],
    "ch": [("ee", 10), (("ch", "e"), 20), ("sh", 90), ("ckh", 97), ("cth", 100)],
    "sh": [("ee", 10), ("ch", 90), ("ckh", 97), ("cth", 100)],
    "ckh": [("cth", 30), (("k", "ch"), 50), (("t", "ch"), 70), ("eke", 72),
            ("ete", 74), ("cph", 78), ("ch", 100)],
    "cth": [("ckh", 30), (("k", "ch"), 50), (("t", "ch"), 70), ("eke", 72),
            ("ete", 74), ("cph", 78), ("cfh", 80), ("ch", 100)],
    "cs": [("sh", 100)],
    "ckhh": [("ckh", 100)],
    "cthh": [("cth", 100)],
    "ikh": [("ckh", 100)],
    "ith": [("cth", 100)],
    "iph": [("cph", 100)],
    "ifh": [("cfh", 100)],
    "eke": [("ckh", 30), ("cth", 50), (("k", "ee"), 56), (("t", "ee"), 60),
            ("ete", 65), ("ee", 100)],
    "ete": [("ckh", 30), ("cth", 50), (("k", "ee"), 56), (("t", "ee"), 60),
            ("eke", 65), ("ee", 100)],
    "cph": [("ckh", 40), ("cth", 75), ("cfh", 80), ("ch", 100)],
    "cfh": [("ckh", 40), ("cth", 75), ("cph", 80), ("ch", 100)],
    "y": [("o", 100)],
    "o": [("y", 100)],
    "n": [("r", 50), ("in", 63), ("iin", 100)],
    "l": [("r", 100)],
    "r": [("r", 50), ("s", 100)],
    "g": [("m", 100)],
    "s": [("r", 75), ("d", 100)],
    "d": [("d", 90), ("s", 100)],
    "a": [("a", 98), ("o", 100)],
    "qo": [("o", 80), ("y", 100)],
}
_F = {  # word-final substitutions (Config.use_word_final_substitutions)
    "om": [("o", 8), ("y", 30), ("ol", 75), ("or", 100)],
    "am": [("o", 4), ("y", 30), ("al", 70), ("ar", 100)],
    "og": [("o", 8), ("y", 30), ("or", 55), ("al", 80), ("ar", 100)],
    "ag": [("o", 4), ("y", 30), ("ol", 65), ("or", 80), ("ar", 100)],
    "ol": [("o", 8), ("y", 30), ("or", 54), ("al", 81), ("ar", 100)],
    "or": [("o", 8), ("y", 30), ("ol", 63), ("al", 81), ("ar", 100)],
    "al": [("o", 4), ("y", 30), ("ol", 64), ("or", 81), ("ar", 100)],
    "ar": [("o", 4), ("y", 30), ("ol", 64), ("or", 82), ("al", 100)],
    "y": [("o", 20), ("ol", 44), ("or", 60), ("al", 75), ("ar", 100)],
    "o": [("y", 20), ("ol", 44), ("or", 60), ("al", 75), ("ar", 100)],
    "d": [("dy", 70), ("d", 100)],
    "dy": [("d", 10), ("dy", 100)],
}


def _norm(table):
    return {
        key: [((v,) if isinstance(v, str) else tuple(v), p) for v, p in rows]
        for key, rows in table.items()
    }


SUBST = _norm(_S)
FINAL_SUBST = _norm(_F)

# tokens to which a given prefix may be attached (Glyph.prefixGlyphs)
PREFIX_ALLOWED = {
    "l": ("k", "t", "p", "f", "d", "ch", "sh", "o", "a", "e", "i"),
    "o": ("k", "t", "p", "f", "d", "ch", "sh"),
    "y": ("k", "t", "p", "f", "d", "ch", "sh"),
    "ch": ("k", "t", "p", "f", "d", "ol", "or", "al", "ar"),
    "sh": ("k", "t", "p", "f", "d", "ol", "or", "al", "ar"),
    "q": (),
    "d": ("a",),
    "x": ("ol", "or", "al", "ar"),
}
LINE_FINAL_REPLACEMENT = {"ol": "om", "or": "om", "al": "am", "ar": "am",
                          "om": "og", "am": "ag", "im": "mg", "in": "n",
                          "iin": "im", "iiin": "im"}
INGROUP_REPLACEMENT = {"y": "o", "m": "r", "g": "r", "n": "r", "in": "ir",
                       "iin": "iir", "dy": "da"}
GROUP_FINAL_REPLACEMENT = {"o": "y", "a": "y", "k": "t", "t": "k",
                           "p": "f", "f": "p"}

# curve/line glyph-shape grammar (CurveLineCanFollow)
CF_START = ("qo", "a", "o", "y", "c", "s", "d", "k", "t", "p", "f", "x")
CF_FINAL = ("y", "n", "l", "r", "s", "d", "m", "g", "x")
CF_CC = ("e", "h", "d", "s", "y", "o", "ch", "sh", "ckh", "cth", "cph", "cfh",
         "al", "ol", "x", "l")
CF_CFINAL = ("d", "g", "o", "y", "dy", "s", "om", "am", "og", "ag", "al",
             "ar", "ol", "or", "x")
CF_CL = ("a",)
CF_LC = ("ikh", "ith", "iph", "ifh")
CF_LL = ("i",)
CF_LFINAL = ("n", "in", "iin", "iiin", "r", "ir", "iir", "iiir", "m", "im",
             "iiil", "iil", "il", "iis", "is")
CF_AFTER_GALLOW = ("a", "e", "o", "y", "h", "ch", "sh")
CF_BEFORE_GALLOW = ("a", "e", "o", "l", "y", "h")
CF_AOY = ("a", "o", "y")
CF_RMNG = ("r", "m", "n", "g")

# Timm's default seed: line f103v.P.9 (Currier B)
DEFAULT_SEED_LINE = "pchal shal shorchdy okeor okain shedy pchedy qotchedy qotar ol lkar"


def tokenize(word: str) -> list[str]:
    """Split an EVA string into Timm's ligature tokens (longest match)."""
    out, position = [], 0
    while position < len(word):
        for lig in LIGATURES:
            if word.startswith(lig, position):
                out.append(lig)
                position += len(lig)
                break
        else:
            out.append(word[position])
            position += 1
    return out


def _end_type(token: str) -> str:
    if token == "":
        return "EMPTY"
    if token.endswith(CF_CL) or token.endswith(CF_LL):
        return "LINE"
    if token.endswith(CF_CC) or token.endswith(CF_LC):
        return "CURVE"
    if token.endswith(CF_FINAL):
        return "FINAL"
    if token.endswith(GALLOWS):
        return "GALLOW"
    return "NONE"


def _start_type(token: str) -> str:
    if token == "":
        return "EMPTY"
    if token.startswith(CF_LL) or token.startswith(CF_LC) or token.startswith(CF_LFINAL):
        return "LINE"
    if token.startswith(CF_CC) or token.startswith(CF_CL) or token.startswith(CF_CFINAL):
        return "CURVE"
    if token.startswith(GALLOWS):
        return "GALLOW"
    return "NONE"


def can_before(add: str, group2: str) -> bool:
    """May ``add`` stand immediately in front of ``group2``?"""
    if add in ALL_COMBINABLE and group2.startswith(CF_AOY):
        return True
    if group2.startswith(ALL_COMBINABLE) and add.endswith(CF_RMNG):
        return True
    if add and group2.startswith(add):
        return False
    if _end_type(group2) == "EMPTY":
        return add in CF_LFINAL or add in CF_CFINAL
    kind = _end_type(add)
    if kind == "EMPTY":
        return group2.startswith(CF_START)
    if kind == "LINE":
        return group2.startswith(CF_LC + CF_LFINAL + CF_LL)
    if kind == "CURVE":
        return group2.startswith(CF_CC + CF_CL + CF_CFINAL + GALLOWS)
    if kind == "GALLOW":
        return group2.startswith(CF_AFTER_GALLOW)
    if kind == "FINAL":
        return group2 == ""
    return False


def can_after(group1: str, add: str) -> bool:
    """May ``add`` follow immediately after ``group1``?"""
    if add in ALL_COMBINABLE and group1.endswith(CF_RMNG):
        return True
    if group1.endswith(ALL_COMBINABLE) and add.startswith(CF_AOY):
        return True
    if add and group1.endswith(add):
        return False
    if _end_type(group1) == "EMPTY":
        return add in CF_START
    kind = _start_type(add)
    if kind == "EMPTY":
        return group1.endswith(CF_FINAL)
    if kind == "LINE":
        return group1.endswith(CF_LL + CF_CL)
    if kind == "CURVE":
        return group1.endswith(CF_CC + CF_LC + GALLOWS)
    if kind == "GALLOW":
        return group1 == "" or group1.endswith(CF_BEFORE_GALLOW)
    return False


class G:
    """A generated glyph group (token) with its ligature parse and origin."""

    __slots__ = ("s", "toks", "kind")

    def __init__(self, toks, kind="INITIAL"):
        if isinstance(toks, str):
            self.s = toks
            self.toks = tokenize(toks)
        else:
            self.s = "".join(toks)
            self.toks = tokenize(self.s)
        self.kind = kind

    def __eq__(self, other):
        return isinstance(other, G) and other.s == self.s

    def __hash__(self):
        return hash(self.s)

    def __repr__(self):
        return f"G({self.s})"

    def is_type_i(self):
        return "i" in self.s

    def is_type_ol(self):
        return any(t in COMBINABLE for t in self.toks)

    def is_type_dy(self):
        return "dy" in self.toks or self.s[-1] in "yd"

    def contains_gallow(self):
        return any(c in self.s for c in GALLOWS)


class SelfCitation:
    """Python port (simplified) of the Timm & Schinner self-citation generator.

    Free parameters (calibrated by the driver, defaults are the Java values):
      p_copy          probability of an unchanged copy of the source (rule VII;
                      the Java code has no explicit copy path, so 0 there);
      replace_weights weights of 1x / 2x / 3x chained similar-glyph
                      replacements (Java: 30/50/20);
      p_add_remove, p_combine_split   morph-method mixture (Java: 20 / 30, rest
                      is replace);
      p_same_position probability of taking the source from the same writing
                      position of a previous line (Java: 28);
      p_line_final    probability that the last token of a line receives the
                      line-final substitution ol->om etc. (Java applies it only
                      when a token does not fit into the 55-character line; we
                      generate to a token-count template, so this is a
                      simplification);
      p_prev_page     (experiment 018, voytools extension) probability that the
                      source line of the "local" mode is taken from an earlier
                      page instead of the current one; 0 = original behaviour
                      (sources only from the current page, bit-identical output);
      page_memory     mean of the geometric distribution (in pages back, >= 1) of
                      the earlier page when ``p_prev_page`` fires; an offset past
                      the first page is clamped to the first page.
    """

    def __init__(self, seed: int, p_copy=0.0, replace_weights=(30, 50, 20),
                 p_add_remove=20, p_combine_split=30, p_same_position=28,
                 p_reuse_last=10, p_line_final=0.35, p_suggest=40,
                 max_repeat=3, seed_line=DEFAULT_SEED_LINE,
                 p_prev_page=0.0, page_memory=1.0):
        if not 0.0 <= float(p_prev_page) <= 1.0:
            raise ValueError("p_prev_page must be in [0, 1]")
        if float(page_memory) < 1.0:
            raise ValueError("page_memory must be >= 1 (mean number of pages back)")
        self.rng = random.Random(seed)
        self.p_prev_page = float(p_prev_page)
        self.page_memory = float(page_memory)
        self.p_copy = p_copy
        self.replace_weights = tuple(replace_weights)
        self.p_add_remove = p_add_remove
        self.p_combine_split = p_combine_split
        self.p_same_position = p_same_position
        self.p_reuse_last = p_reuse_last
        self.p_line_final = p_line_final
        self.p_suggest = p_suggest
        self.max_repeat = max_repeat
        self.seed_line = [G(w) for w in seed_line.split()]
        # state
        self.lines: list[list[G]] = []
        self.page_lines: list[list[G]] = []
        self.prev_pages: list[list[list[G]]] = []   # finished pages (018 memory extension)
        self.par_initial_lines: list[list[G]] = []
        self.freq: Counter = Counter()
        self.n_tokens = 0
        self.n_i = self.n_dy = self.n_ol = 0
        self.repeat_kind = None
        self.repeat_count = 0

    # ---------------------------------------------------------------- utils
    def rand(self, n: int) -> int:
        return self.rng.randrange(n) if n > 0 else 0

    def chance(self, percent: float) -> bool:
        return self.rng.random() * 100 < percent

    def random_gallow(self, first_line: bool) -> str:
        r = self.rand(100)
        if r <= 34:
            return "k"
        if r <= 49:
            return "t"
        if r <= 89:
            return "p" if first_line else "k"
        return "f" if first_line else "t"

    def random_line_initial(self) -> str:
        r = self.rand(100)
        return "o" if r <= 45 else "y" if r <= 75 else "d" if r <= 90 else "s"

    def remember(self, g: G) -> None:
        self.freq[g.s] += 1
        self.n_tokens += 1
        if g.is_type_i():
            kind = "i"
            self.n_i += 1
        elif g.is_type_dy():
            kind = "dy"
            self.n_dy += 1
        elif g.is_type_ol():
            kind = "ol"
            self.n_ol += 1
        else:
            kind = "other"
        if kind == self.repeat_kind:
            self.repeat_count += 1
        else:
            self.repeat_kind, self.repeat_count = kind, 1

    # --------------------------------------------------------- source choice
    def choose_source(self, current: list[G], is_par_initial: bool,
                      is_line_initial: bool) -> list[G]:
        """Return [source, neighbour] following PageSourceGroupChooser."""
        mode = "local"
        if not self.lines:
            mode = "random"
        if (is_par_initial and len(self.par_initial_lines) > 1
                and (is_line_initial or (len(current) > 1 and self.chance(70)))):
            mode = "paragraph_initial"
        if self.suggestion_needed() and self.chance(self.p_suggest):
            mode = "suggestion"
        if mode == "random":
            pool = self.seed_line
            first = self.rng.choice(pool)
            second = self.rng.choice(pool)
            ret = [first, second]
        elif mode == "paragraph_initial":
            source_line = self.rng.choice(self.par_initial_lines)
            pos = self.rand(len(source_line))
            ret = self._pair(source_line, pos)
        elif mode == "suggestion":
            ret = [self.suggest(), self.rng.choice(self.seed_line)]
        else:
            # one of the previous lines on this page (or the last line of the
            # previous page at a page top)
            candidates = self.page_lines if self.page_lines else self.lines[-1:]
            if self.p_prev_page > 0 and self.prev_pages and self.rng.random() < self.p_prev_page:
                candidates = self.prev_pages[-self.page_offset()]
            source_line = candidates[self.rand(len(candidates))]
            probability = self.p_same_position / 2 if is_line_initial else self.p_same_position
            if self.rand(100) <= probability:
                pos = self._same_position(current, source_line)
            else:
                pos = self.rand(len(source_line))
            ret = self._pair(source_line, pos)
        # ChooserHelper.removeInitialGallow
        out = []
        for g in ret:
            if g.toks and g.toks[0] in GALLOWS and len(g.toks) > 1:
                out.append(G(g.toks[1:], g.kind))
            else:
                out.append(g)
        return out

    def page_offset(self) -> int:
        """Pages back (>= 1) of the memory source: geometric with mean ``page_memory``, clamped to the first page."""
        q = 1.0 / self.page_memory
        if q >= 1.0:
            k = 1
        else:
            u = self.rng.random()
            k = 1 + int(math.log(max(u, 1e-300)) / math.log(1.0 - q))
        return min(max(k, 1), len(self.prev_pages))

    @staticmethod
    def _pair(line: list[G], pos: int) -> list[G]:
        if pos < len(line) - 1:
            return [line[pos], line[pos + 1]]
        if pos > 0:
            return [line[pos], line[pos - 1]]
        return [line[pos], line[pos]]

    @staticmethod
    def _same_position(current: list[G], source_line: list[G]) -> int:
        writing = sum(len(g.s) + 1 for g in current)
        position = 0
        for index, g in enumerate(source_line):
            position += len(g.s) + 1
            if writing <= position:
                return index
        return min(len(source_line) - 1, len(current))

    def suggestion_needed(self) -> bool:
        if self.n_tokens < 50:
            return False
        return (self.n_i / self.n_tokens < 0.20) or (self.n_dy / self.n_tokens < 0.25)

    def suggest(self) -> G:
        """Most frequent i-type or dy-type token so far (StatisticHelper)."""
        want = "i" if self.n_i / self.n_tokens < 0.20 else "dy"
        for word, _ in self.freq.most_common():
            g = G(word)
            if (want == "i" and g.is_type_i()) or (want == "dy" and g.is_type_dy()):
                return g
        return G("daiin") if want == "i" else G("chedy")

    # ------------------------------------------------------------- morphing
    def morph(self, sources: list[G], previous: G | None, is_par_initial: bool,
              is_line_initial: bool) -> list[G]:
        source = sources[0]
        length = len(source.toks)
        ret: list[G] = []
        if self.p_copy > 0 and self.rng.random() < self.p_copy:
            ret = [G(source.toks, "COPY")]
        else:
            r = 0 if (is_par_initial and is_line_initial) else self.rand(100)
            if r <= self.p_add_remove and source.kind != "COMBINE":
                method = "add_remove"
            elif r <= self.p_add_remove + self.p_combine_split:
                method = "combine_split"
            else:
                method = "replace"
            if method == "add_remove":
                morphed = (self.add_random_glyph(source, previous, is_par_initial, is_line_initial)
                           if length < 6 else source)
                if morphed != source:
                    ret.append(morphed)
                else:
                    morphed = self.delete_prefix(source)
                    if morphed != source:
                        ret.append(morphed)
            elif method == "combine_split":
                r = self.rand(100)
                if source.kind != "COMBINE":
                    if length < 6:
                        combine = length <= 2 or r < 96
                    else:
                        combine = r < 4 and length <= 8
                else:
                    combine = False
                if combine and len(sources) > 1:
                    morphed = self.combine(sources)
                    if morphed != source:
                        ret.append(morphed)
                else:
                    ret.extend(self.split(source))
            else:
                w1, w2, w3 = self.replace_weights
                r = self.rng.random() * (w1 + w2 + w3)
                repeats = 1 if r < w1 else 2 if r < w1 + w2 else 3
                temp = source
                for _ in range(repeats):
                    nxt = self.replace_random(temp, is_par_initial)
                    if nxt == source:
                        nxt = temp
                    temp = nxt
                if temp == source:  # one more try, as in the Java 1x branch
                    temp = self.replace_random(source, is_par_initial)
                if temp != source:
                    ret.append(temp)
        if not ret:
            return ret
        self.handle_gallows(is_par_initial, is_line_initial, ret)
        if self.p_reuse_last > 0:
            self.reuse_last(is_par_initial, ret)
        return ret

    def add_random_glyph(self, g: G, previous, is_par_initial, is_line_initial) -> G:
        if self.rand(100) < (80 if is_par_initial else 8):
            return self.add_gallow(g, is_par_initial, is_line_initial)
        return self.add_prefix(g, previous)

    def choose_prefix(self, j: int, previous, tried: list[str]) -> str:
        r = self.rand(100)
        if j == 0 and previous is not None and previous.toks and previous.toks[-1] == "dy" and r < 90:
            return "q"
        if r < 5 and "l" not in tried:
            return "l"
        if r < 34 and "o" not in tried:
            return "o"
        if r < 40 and "y" not in tried:
            return "y"
        if r < 61 and "ch" not in tried:
            return "ch"
        if r < 72 and "sh" not in tried:
            return "sh"
        if r < 88 and "q" not in tried:
            return "q"
        if r < 99:
            return "d"
        if j == 0:
            return "x"
        return "o"

    def add_prefix(self, g: G, previous) -> G:
        tried: list[str] = []
        for j in range(7):
            prefix = self.choose_prefix(j, previous, tried)
            tried.append(prefix)
            start = g.toks[0]
            if prefix == "q":
                if start in ("o", "y"):
                    return G(["qo"] + g.toks[1:], "ADD" if g.kind != "COMBINE" else "COMBINE")
            elif prefix == "x":
                if start in PREFIX_ALLOWED["x"]:
                    return G(["x"] + g.toks, "ADD" if g.kind != "COMBINE" else "COMBINE")
            else:
                if start in PREFIX_ALLOWED[prefix]:
                    toks = list(g.toks)
                    # o + daiin -> okaiin, y + chol -> ykol
                    if len(toks) > 2 and start in ("d", "ch", "sh"):
                        gallow = self.random_gallow(False)
                        if can_before(gallow, toks[1]) and can_after(prefix, gallow):
                            if self.rand(100) < 70:
                                toks[0] = gallow
                    return G([prefix] + toks, "ADD" if g.kind != "COMBINE" else "COMBINE")
        return g

    def add_gallow(self, g: G, is_par_initial: bool, is_line_initial: bool) -> G:
        length = len(g.toks)
        if is_par_initial and is_line_initial:
            gallow = self.random_gallow(True)
            out = self.place_gallow(gallow, g, 0)
            if out.toks[0] != gallow:
                toks = list(out.toks)
                if can_before("o", out.s):
                    toks.insert(0, "o")
                elif can_before("a", out.s):
                    toks.insert(0, "a")
                toks.insert(0, gallow)
                out = G(toks, out.kind)
            return out
        if length > 1:
            if g.contains_gallow() and self.rand(100) < 90:
                return g
            for _ in range(5):
                if is_par_initial:
                    pos = self.rand(length)
                elif length == 2:
                    pos = 1
                else:
                    pos = self.rand(length - 2) + 1
                gallow = self.random_gallow(is_par_initial)
                out = self.place_gallow(gallow, g, pos)
                if out != g:
                    return out
        return g

    def place_gallow(self, gallow: str, g: G, pos: int) -> G:
        toks = list(g.toks)
        if toks[pos] in GALLOWS:
            toks[pos] = gallow
            return G(toks, g.kind)
        last = "" if pos == 0 else toks[pos - 1]
        nxt = "" if pos == len(toks) else toks[pos]
        last_ok = can_after(last, gallow)
        next_ok = can_before(gallow, nxt)
        if last_ok and next_ok:
            toks.insert(pos, gallow)
            return G(toks, g.kind)
        if last_ok and len(toks) > pos + 1:
            alt = toks[pos + 1]
            toks[pos] = gallow
            if can_before(gallow, alt):
                return G(toks, g.kind)
            return g
        if next_ok and pos - 2 > 0:
            toks[pos - 1] = gallow
            if can_after(toks[pos - 2], gallow):
                return G(toks, g.kind)
        return g

    def delete_prefix(self, g: G) -> G:
        toks = list(g.toks)
        start = toks.pop(0)
        if len(g.toks) > 2 and can_before("", start):
            first, second = toks[0], toks[1]
            if first in GALLOWS:
                r = self.rand(100)
                if second.startswith("a") and r < 50:
                    toks[0] = "d"
                if (second.startswith("e") or second.startswith("o")) and r < 50:
                    toks[0] = "ch"
            return G(toks, "DELETE")
        return g

    def combine(self, sources: list[G]) -> G:
        g1, g2 = sources[0], sources[1]
        toks1 = self._subgroup_for_combine(g1, True)
        if toks1 and len("".join(toks1)) > 1:
            last1 = toks1[-1]
            toks2 = self._subgroup_for_combine(g2, False)
            drop = LIGATURE_TAIL_DROP.get(last1)
            if drop and toks2 and toks2[0] in drop:
                toks2 = toks2[1:]
            if toks2 and len("".join(toks2)) > 1:
                merged = toks1 + toks2
                if len("".join(merged)) < 9:
                    if ((last1 in COMBINABLE or can_before(last1, toks2[0]))
                            and can_after(toks2[-1], "")):
                        return G(merged, "COMBINE")
        return self._combine_fallback(g1, g2)

    def _subgroup_for_combine(self, g: G, is_first: bool) -> list[str]:
        n = len(g.toks)
        if (is_first and n <= 2) or (not is_first and n <= 3):
            return list(g.toks)
        pos = self._split_position(g)
        if pos >= 0:
            parts = self._split_at(g, pos)
            if parts:
                if is_first or len(parts) == 1:
                    return list(parts[0].toks)
                return list(parts[1].toks)
        return []

    def _combine_fallback(self, g1: G, g2: G) -> G:
        if len(g1.toks) == 2 or (len(g1.toks) == 1 and g1.s in COMBINABLE):
            if self.rand(100) < 60:
                out = self.self_combine(g1)
                if out != g1:
                    return out
        toks2 = self._choose_subgroups(g2)
        if toks2:
            toks1 = list(g1.toks)
            last1 = toks1[-1]
            if last1 in COMBINABLE:
                drop = LIGATURE_TAIL_DROP.get(last1)
                if drop and toks2 and toks2[0] in drop:
                    toks2 = toks2[1:]
            if toks2 and len(toks1) + len(toks2) < 9:
                merged = toks1 + toks2
                if can_before(last1, toks2[0]) and can_after(toks2[-1], ""):
                    return G(merged, "COMBINE")
        if len(g1.toks) <= 2:
            return self.self_combine(g1)
        return g1

    def _choose_subgroups(self, g: G) -> list[str]:
        r = self.rand(100)
        n = len(g.toks)
        if r <= 39:
            if n > 1:
                k = self.rand(n - 1) + 1
                return list(g.toks[:k])
            return []
        if r <= 59:
            return list(g.toks[:-1])
        for t in g.toks:
            if t in COMBINABLE:
                return [t]
        return list(g.toks[:-1])

    def self_combine(self, g: G) -> G:
        """olol, chochy: duplicate a short group (rule VI)."""
        if not g.toks:
            return g
        first, last = g.toks[0], g.toks[-1]
        if not can_before(last, first):
            return g
        toks = list(g.toks)
        toks1 = list(g.toks)
        final = GROUP_FINAL_REPLACEMENT.get(last, last)
        ingroup = INGROUP_REPLACEMENT.get(last, last)
        if final != last and self.rand(100) < 80:
            toks[-1] = final
        if ingroup != last:
            toks1[-1] = ingroup
        if can_before("k", first) and can_before(last, "k") and self.rand(100) < 30:
            toks.insert(1, self.random_gallow(False))
        return G(toks1 + toks, "COMBINE")

    def _split_position(self, g: G) -> int:
        for i in range(1, len(g.toks)):
            last, nxt = g.toks[i - 1], g.toks[i]
            if not can_before(last, nxt) or last in COMBINABLE or nxt in GALLOWS:
                if i > 1 or len(last) > 1:
                    return i
        return -1

    def _split_at(self, g: G, pos: int) -> list[G]:
        out = []
        part1 = G(g.toks[:pos], "SPLIT")
        if (len(part1.toks) > 1 or part1.s in COMBINABLE) and can_after(part1.toks[-1], ""):
            out.append(part1)
        part2 = G(g.toks[pos:], "SPLIT")
        if len(part2.toks) > 1:
            out.append(part2)
        return out

    def split(self, g: G) -> list[G]:
        pos = self._split_position(g)
        out = self._split_at(g, pos) if pos >= 0 else []
        return out

    def replace_random(self, g: G, is_par_initial: bool) -> G:
        positions = [i for i in range(len(g.toks)) for _ in range(2)]
        j = -1
        # NB: the Java loop condition ``j < posToTryList.size()`` is re-evaluated
        # while the list shrinks, so an n-token group gets n tries, and the
        # word-final table (j > 0) never applies to a one-token group.
        while True:
            j += 1
            if j >= len(positions):
                break
            pos = positions.pop(self.rand(len(positions))) if len(positions) > 1 else positions.pop()
            sub = g.toks[pos]
            is_final = j > 0 and pos == len(g.toks) - 1
            table = FINAL_SUBST.get(sub) if is_final else None
            if table is None:
                table = SUBST.get(sub)
            if table is None:
                continue
            r = self.rand(100)
            substitute = next(v for v, p in table if p > r)
            toks = list(g.toks)
            nxt = toks[pos + 1] if len(toks) > pos + 1 else ""
            remove_next = any(c in GALLOWS for c in "".join(substitute)) and nxt in GALLOWS
            if self._replacable(substitute, g, pos, remove_next):
                del toks[pos]
                if remove_next:
                    del toks[pos]
                toks[pos:pos] = list(substitute)
                out = G(toks, "COMBINE" if g.kind == "COMBINE" else "REPLACE")
                if out != g:
                    return out
        return g

    def _replacable(self, substitute, g: G, pos: int, remove_next: bool) -> bool:
        last = "" if pos == 0 else g.toks[pos - 1]
        next_pos = pos + 2 if remove_next else pos + 1
        nxt = "" if next_pos >= len(g.toks) else g.toks[next_pos]
        last_ok = can_after(last, substitute[0])
        next_ok = can_before(substitute[-1], nxt)
        # "multiple instances of the same token within a group are rare"
        # (Java tests String.contains on the whole group)
        if (last_ok and next_ok and substitute[0] in g.s
                and substitute[0] not in COMBINABLE):
            compare = 100 if (substitute[0] == last or substitute[-1] == nxt) else 80
            if self.rand(100) < compare:
                return False
        return last_ok and next_ok

    def handle_gallows(self, is_par_initial: bool, is_line_initial: bool, ret: list[G]) -> None:
        if is_line_initial:
            if is_par_initial and self.rand(100) < 94:
                ret[0] = self.add_gallow(ret[0], True, True)
            else:
                ret[0] = self.add_line_initial(ret[0])
        elif ret[0].toks[0] in LINE_INITIAL:
            if self.rand(len(PREFIX_ALLOWED)) < 3 and len(ret[0].toks) > 3:
                ret[0] = self.remove_line_initial(ret[0])
        if ret[0].contains_gallow():
            if not is_par_initial:
                if self.rand(100) < 94:
                    ret[0] = G(ret[0].s.replace("p", "k").replace("f", "t")
                               if self.rand(2) else
                               ret[0].s.replace("p", "t").replace("f", "k"), ret[0].kind)
            else:
                gallow = self.random_gallow(True)
                s = ret[0].s
                for old in GALLOWS:
                    s = s.replace(old, gallow)
                ret[0] = G(s, ret[0].kind)

    def add_line_initial(self, g: G) -> G:
        if g.toks[0] in LINE_INITIAL or self.rand(100) >= 30:
            return g
        glyph = self.random_line_initial()
        toks = list(g.toks)
        if len(toks) > 2 and glyph in ("o", "y") and toks[0] in ("d", "ch"):
            if self.rand(100) < 90:
                gallow = self.random_gallow(False)
                if can_before(gallow, toks[1]):
                    toks[0] = gallow
        toks.insert(0, glyph)
        out = G(toks, g.kind)
        if can_before(glyph, toks[1] if len(toks) > 1 else ""):
            return out
        return g

    def remove_line_initial(self, g: G) -> G:
        if g.toks[0] in LINE_INITIAL and len(g.toks) > 2:
            toks = list(g.toks[1:])
            first, second = toks[0], toks[1]
            if first in GALLOWS:
                if second.startswith("a"):
                    toks[0] = "d"
                if second.startswith("e"):
                    toks[0] = "ch"
            return G(toks, g.kind)
        return g

    def reuse_last(self, is_par_initial: bool, ret: list[G]) -> None:
        """chol -> chol or: derive a second token from the one just made."""
        if len(ret) > 1 or self.rand(100) >= self.p_reuse_last:
            return
        last = ret[0]
        lig = next((t for t in last.toks if t in COMBINABLE), None)
        if lig is not None:
            add = G(lig, "SPLIT")
            if self.rand(100) < 50:
                add = self.replace_random(add, is_par_initial)
            ret.append(add)
        else:
            morphed = self.delete_prefix(last)
            if morphed != last:
                if self.rand(100) < 50:
                    morphed = self.replace_random(morphed, is_par_initial)
                ret.append(morphed)

    # ------------------------------------------------------------ generation
    def generate(self, template: list[tuple[str, int, bool, bool]]) -> list[list[str]]:
        """Generate one token line per template row (page, n_tokens,
        is_paragraph_initial, is_paragraph_final)."""
        for g in self.seed_line:
            self.remember(g)
        output = []
        current_page = None
        for page, n_tokens, is_par_initial, is_par_final in template:
            if page != current_page:
                current_page = page
                if self.page_lines:
                    self.prev_pages.append(self.page_lines)
                self.page_lines = []
            line = self.generate_line(n_tokens, is_par_initial, is_par_final)
            self.lines.append(line)
            self.page_lines.append(line)
            if is_par_initial:
                self.par_initial_lines.append(line)
            output.append([g.s for g in line])
        return output

    def generate_line(self, n_tokens: int, is_par_initial: bool, is_par_final: bool) -> list[G]:
        line: list[G] = []
        last_sources = None
        last_generated = None
        attempts = 0
        while len(line) < n_tokens:
            attempts += 1
            is_line_initial = not line
            sources = self.choose_source(line, is_par_initial, is_line_initial)
            use = True
            if len(sources[0].toks) > 5:
                use = len(sources[0].s) < 4 + self.rand(4)
            if use and last_sources is not None and last_sources[0] == sources[0]:
                use = False
            if use and sources[0].kind == "COMBINE":
                if sources[1].kind == "COMBINE" or self.rand(100) < 30:
                    use = False
            if attempts > 100:
                sources = [self.suggest() if self.n_tokens else self.seed_line[0],
                           self.rng.choice(self.seed_line)]
                use = True
            force = attempts > 110
            if attempts > 400:
                # emergency: copy the seed verbatim rather than loop forever
                line.append(G(self.rng.choice(self.seed_line).toks, "COPY"))
                self.remember(line[-1])
                attempts = 0
                continue
            if not use:
                continue
            last_sources = sources
            morphed = self.morph(sources, last_generated, is_par_initial, is_line_initial)
            if not morphed:
                continue
            first = morphed[0]
            ok = can_before("", first.s)
            if line and first == line[-1]:
                ok = self.rand(100) < 50
            elif self.repeat_count > self.max_repeat:
                kind = ("i" if first.is_type_i() else "dy" if first.is_type_dy()
                        else "ol" if first.is_type_ol() else "other")
                if kind == self.repeat_kind:
                    ok = False
            if not (ok or force):
                continue
            for g in morphed:
                if len(line) >= n_tokens:
                    break
                line.append(g)
                self.remember(g)
                last_generated = g
            attempts = 0
        # line-final substitution (Java: tryToTrim when the line is full)
        if line and self.rng.random() < self.p_line_final:
            last = line[-1]
            new_last = LINE_FINAL_REPLACEMENT.get(last.toks[-1])
            if new_last is not None:
                line[-1] = G(last.toks[:-1] + [new_last], last.kind)
        return line


# ---------------------------------------------------------------------------
# voytools API
# ---------------------------------------------------------------------------


def _template_rows(template_lines: Sequence[Sequence[Any]]) -> list[tuple[str, int, bool, bool]]:
    rows = []
    for row in template_lines:
        row = tuple(row)
        if len(row) == 2:
            page, (n, first), last = "p", row, False
        elif len(row) == 3:
            (page, n, first), last = row, False
        elif len(row) == 4:
            page, n, first, last = row
        else:
            raise ValueError("template rows are (n_tokens, is_par_first_line), "
                             "(page, n_tokens, is_par_first_line) or (..., is_par_last_line)")
        if int(n) < 0:
            raise ValueError("n_tokens must be >= 0")
        rows.append((str(page), int(n), bool(first), bool(last)))
    return rows


def selfcitation_generate(
    template_lines: Sequence[Sequence[Any]],
    *,
    seed: int,
    p_copy: float = 0.0,
    replace_weights: tuple[float, float, float] = (30, 50, 20),
    p_add_remove: float = 20,
    p_combine_split: float = 30,
    p_same_position: float = 28,
    p_reuse_last: float = 10,
    p_line_final: float = 0.35,
    p_suggest: float = 40,
    max_repeat: int = 3,
    seed_line: str = DEFAULT_SEED_LINE,
    p_prev_page: float = 0.0,
    page_memory: float = 1.0,
) -> list[list[str]]:
    """Generate text with the Timm & Schinner self-citation mechanism (port, see module docstring).

    Parameters
    ----------
    template_lines:
        One row per output line, in text order: ``(page, n_tokens, is_par_first_line)``, optionally
        with a fourth element ``is_par_last_line`` (unused by the mechanism), or
        ``(n_tokens, is_par_first_line)`` (all lines on one page).  The page matters: sources are
        taken from the previous lines of the same page.  :func:`selfcitation_template` builds the
        rows from a real token table.
    seed:
        Seed of the generator's ``random.Random`` (output is deterministic for a given seed).
    p_copy:
        Probability of an unchanged copy of the source (rule VII; the Java code has none, 0).
    replace_weights:
        Weights of 1x / 2x / 3x chained similar-glyph replacements (Java: 30/50/20).
    p_add_remove, p_combine_split:
        Percent of the morph-method mixture (add/remove glyph, combine/split tokens; the rest is
        replacement).  Java: 20 / 30.
    p_same_position:
        Percent chance of taking the source from the same writing position of a previous line
        (halved at line starts).  Java: 28.
    p_reuse_last:
        Percent chance of deriving a second token from the one just made.
    p_line_final:
        Probability (0-1) of the line-final substitution (``ol -> om`` ...) on the last token.
    p_suggest:
        Percent chance of taking the most frequent i-/dy-type token when those are under-represented.
    max_repeat:
        Maximum run of tokens of the same type (i / dy / ol / other).
    seed_line:
        Initial line the generator starts copying from (Timm's default: f103v.P.9).
    p_prev_page, page_memory:
        Memory across pages (experiment 018; not in Timm & Schinner): with probability
        ``p_prev_page`` (0-1) the "local" source line is taken from an earlier page, the number of
        pages back being geometric with mean ``page_memory`` (>= 1; clamped to the first page).
        The defaults (0, 1) give exactly the original output for every seed.

    Returns
    -------
    list of lines (lists of basic-Eva token strings) with exactly the template's token counts.
    """
    rows = _template_rows(template_lines)
    gen = SelfCitation(seed=seed, p_copy=p_copy, replace_weights=tuple(replace_weights),
                       p_add_remove=p_add_remove, p_combine_split=p_combine_split,
                       p_same_position=p_same_position, p_reuse_last=p_reuse_last,
                       p_line_final=p_line_final, p_suggest=p_suggest, max_repeat=max_repeat,
                       seed_line=seed_line, p_prev_page=p_prev_page, page_memory=page_memory)
    return gen.generate(rows)


def selfcitation_template(df: pd.DataFrame) -> list[tuple[str, int, bool, bool]]:
    """Template rows ``(folio, n_tokens, is_par_first_line, is_par_last_line)`` of a token table.

    One row per line (locus) in text order; the token counts are those of ``df`` (after any
    filtering), so a generated sample has exactly the layout of ``df``.
    """
    from ._common import line_ids

    if not len(df):
        return []
    ids = line_ids(df)
    starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
    counts = np.diff(np.r_[starts, len(df)])
    folio = df["folio"].astype(str).to_numpy()
    first = df["is_par_first_line"].fillna(False).astype(bool).to_numpy() if "is_par_first_line" in df else np.zeros(len(df), bool)
    last = df["is_par_last_line"].fillna(False).astype(bool).to_numpy() if "is_par_last_line" in df else np.zeros(len(df), bool)
    return [(folio[s], int(c), bool(first[s]), bool(last[s])) for s, c in zip(starts, counts)]


def selfcitation_like(template_df: pd.DataFrame, **params: Any) -> pd.DataFrame:
    """Self-citation sample with exactly the layout and metadata of a real token table.

    Generates one line per line of ``template_df`` with the same number of tokens, page and
    paragraph-initial flag (:func:`selfcitation_template`), and returns
    :func:`voytools.nulls.as_table` of the result (like :func:`voytools.nulls.autocopy_like`).
    ``params`` go to :func:`selfcitation_generate` (``seed`` is required).  The parameters are
    stored in ``attrs['selfcitation_params']``; ``attrs['source']`` is ``'selfcitation'``.
    """
    from .nulls import as_table

    if "seed" not in params:
        raise TypeError("selfcitation_like() needs a 'seed'")
    lines = selfcitation_generate(selfcitation_template(template_df), **params)
    out = as_table(lines, template_df)
    out.attrs["source"] = "selfcitation"
    out.attrs["selfcitation_params"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in params.items()}
    return out
