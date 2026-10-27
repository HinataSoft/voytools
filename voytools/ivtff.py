"""Parser for the IVTFF (Intermediate Voynich MS Transliteration File Format) 2.0.

The parser follows the format definition by R. Zandbergen (IVTFF 2.0.1/2.0.2):

* line 1 is the file header ``#=IVTFF <alphabet> <version> <origin> ...``;
* ``#`` lines are comments;
* a *page header* ``<f1r>  <! $Q=A $P=A ...>`` opens a page and sets page variables;
* a *locus line* ``<f1r.1,@P0>  text`` carries the transliterated text of one
  locus (``@`` = locator, ``P0`` = locus type, optional ``;H`` = transcriber ID
  in interlinear files);
* lines that start with ``/`` continue the previous line when it ends in ``/``.

The text of a locus is cleaned in two stages.  :func:`parse_ivtff` removes free
comments ``<!...>``, text tags ``<@X=y>`` and the paragraph marks ``<%>`` and
``<$>`` (their information is kept in dedicated :class:`Locus` fields).
:func:`tokenize_locus` then splits the cleaned text into :class:`Token` records
according to a *space policy* (how the uncertain space ``,`` is treated) and an
*alternatives policy* (how ``[a:b]`` readings are treated).

Design notes
------------
* Unknown or odd constructs never raise; they are recorded in
  :attr:`Transliteration.warnings`.  Only a malformed locus identifier (or a
  line starting with ``<`` that has no closing ``>``) raises :class:`IVTFFError`.
* High-ASCII characters ``@nnn;`` (nnn = 128..255) become ONE private-use
  Unicode symbol ``chr(0xE000 + nnn)``; use :func:`render_high_ascii` to turn
  them back into ``@nnn;``.
* Files with format version < 2.0 and an Eva-like alphabet (e.g. the Landini-
  Stolfi interlinear, ``#=IVTFF Eva- 1.5``) use ``!`` and ``%`` as interlinear
  filler characters and spaces for alignment; these are removed by default
  (``strip_fillers='auto'``).
* Locus numbers normally are integers; the interlinear file also has numbers
  like ``44a``.  The numeric part is stored in :attr:`Locus.num`, the letter
  in :attr:`Locus.num_suffix`, and :attr:`Locus.id` reproduces the full text.

Format reference: ``docs/IVTFF_format_2.0.txt`` (Tables 3, 6-11, Section 6).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Iterable, Iterator

__all__ = [
    "HIGH_ASCII_BASE",
    "UNCERTAIN_SPACE_MARK",
    "IVTFFError",
    "ParseWarning",
    "Header",
    "Page",
    "Locus",
    "Token",
    "Transliteration",
    "parse_ivtff",
    "parse_ivtff_text",
    "parse_ivtff_lines",
    "tokenize_locus",
    "tokenize_text",
    "encode_high_ascii",
    "decode_high_ascii",
    "render_high_ascii",
    "is_high_ascii_char",
]

HIGH_ASCII_BASE = 0xE000
"""High-ASCII code ``nnn`` is represented by the single character ``chr(HIGH_ASCII_BASE + nnn)``."""

UNCERTAIN_SPACE_MARK = ","
"""Literal marker kept inside tokens under ``space_policy='strict'``."""

LOCATORS = "@+*=&~/!"
GENERIC_TYPES = "PLCR"
SPACE_POLICIES = ("split", "join", "strict")
ALTERNATIVE_POLICIES = ("first", "mark")

_PUA_RE = re.compile("[-]")
_HIGH_ASCII_RE = re.compile(r"@(\d{3});")


class IVTFFError(ValueError):
    """Raised for a malformed locus identifier or an unreadable file header."""


@dataclass(frozen=True)
class ParseWarning:
    """A non-fatal problem found while parsing."""

    line_no: int
    category: str
    message: str
    locus: str | None = None

    def __str__(self) -> str:
        where = f"line {self.line_no}" + (f" [{self.locus}]" if self.locus else "")
        return f"{where}: {self.category}: {self.message}"


# --------------------------------------------------------------------------- #
# high-ASCII helpers
# --------------------------------------------------------------------------- #


def encode_high_ascii(code: int) -> str:
    """Return the single private-use character that stands for ``@code;``."""
    if not 128 <= code <= 255:
        raise ValueError(f"high-ASCII code must be in 128..255, got {code}")
    return chr(HIGH_ASCII_BASE + code)


def is_high_ascii_char(ch: str) -> bool:
    """True if ``ch`` is one of the private-use symbols made by :func:`encode_high_ascii`."""
    return len(ch) == 1 and HIGH_ASCII_BASE + 128 <= ord(ch) <= HIGH_ASCII_BASE + 255


def decode_high_ascii(text: str) -> str:
    """Replace every ``@nnn;`` (128 <= nnn <= 255) in ``text`` with its private-use symbol."""

    def repl(m: re.Match[str]) -> str:
        n = int(m.group(1))
        return chr(HIGH_ASCII_BASE + n) if 128 <= n <= 255 else m.group(0)

    return _HIGH_ASCII_RE.sub(repl, text)


def render_high_ascii(text: str) -> str:
    """Inverse of :func:`decode_high_ascii`: private-use symbols become ``@nnn;``."""
    return _PUA_RE.sub(lambda m: f"@{ord(m.group()) - HIGH_ASCII_BASE:03d};", text)


# --------------------------------------------------------------------------- #
# data model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Header:
    """The ``#=IVTFF`` file header."""

    alphabet: str | None
    version: str | None
    origin: str | None
    extra: str = ""
    raw: str = ""

    @property
    def version_tuple(self) -> tuple[int, ...]:
        """Format version as a tuple of ints (``(0,)`` when absent)."""
        if not self.version:
            return (0,)
        try:
            return tuple(int(p) for p in self.version.split("."))
        except ValueError:
            return (0,)

    @property
    def is_eva_like(self) -> bool:
        """True for the Eva alphabets (``Eva-``, ``EvaT``)."""
        return bool(self.alphabet) and self.alphabet.lower().startswith("eva")


@dataclass
class Page:
    """A page (``<f1r>``) with its page variables (``$Q=A`` ...)."""

    name: str
    variables: dict[str, str]
    line_no: int
    index: int
    comments: list[str] = field(default_factory=list)


@dataclass
class Locus:
    """One transliteration item (a text line, label, circular or radial text)."""

    folio: str
    num: int
    num_suffix: str
    locator: str
    locus_type: str
    generic_type: str
    subtype: str
    transcriber: str | None
    raw_text: str
    text: str
    comments: list[str]
    paragraph_start: bool
    paragraph_end: bool
    inline_tags: list[tuple[str, str]]
    line_no: int
    page_index: int
    variables: dict[str, str | None]
    continued: bool = False

    @property
    def id(self) -> str:
        """Locus identifier without locator, type and transcriber, e.g. ``f1r.1``."""
        return f"{self.folio}.{self.num}{self.num_suffix}"


@dataclass(frozen=True)
class Token:
    """A token (word) cut out of a locus, with quality and spacing flags."""

    text: str
    raw: str
    index: int
    has_unreadable: bool = False
    has_alternative: bool = False
    has_ligature: bool = False
    has_high_ascii: bool = False
    after_uncertain_space: bool = False
    before_uncertain_space: bool = False
    after_drawing_break: bool = False
    before_drawing_break: bool = False
    n_uncertain_inside: int = 0


@dataclass
class Transliteration:
    """A parsed IVTFF file."""

    header: Header
    pages: list[Page]
    loci: list[Locus]
    warnings: list[ParseWarning] = field(default_factory=list)
    warning_counts: dict[str, int] = field(default_factory=dict)
    comment_lines: list[tuple[int, str]] = field(default_factory=list)
    path: str | None = None

    @property
    def alphabet(self) -> str | None:
        """Alphabet code from the header (``Eva-``, ``v101`` ...)."""
        return self.header.alphabet

    @property
    def transcribers(self) -> list[str]:
        """Sorted transcriber IDs present (empty for non-interlinear files)."""
        return sorted({lc.transcriber for lc in self.loci if lc.transcriber})

    @property
    def is_interlinear(self) -> bool:
        """True when more than one transcriber ID occurs."""
        return len(self.transcribers) > 1

    def page(self, name: str) -> Page:
        """Return the page called ``name`` (KeyError if absent)."""
        for p in self.pages:
            if p.name == name:
                return p
        raise KeyError(name)

    def iter_loci(
        self,
        generic_type: str | Iterable[str] | None = None,
        transcriber: str | None = None,
    ) -> Iterator[Locus]:
        """Iterate over loci, optionally restricted to generic type(s) and one transcriber."""
        types = None
        if generic_type is not None:
            types = {generic_type} if isinstance(generic_type, str) else set(generic_type)
        for lc in self.loci:
            if types is not None and lc.generic_type not in types:
                continue
            if transcriber is not None and lc.transcriber != transcriber:
                continue
            yield lc


# --------------------------------------------------------------------------- #
# warnings collector
# --------------------------------------------------------------------------- #


class _Warner:
    """Collects warnings, storing at most ``limit`` messages per category."""

    def __init__(self, limit: int = 20) -> None:
        self.limit = limit
        self.items: list[ParseWarning] = []
        self.counts: dict[str, int] = {}

    def add(self, line_no: int, category: str, message: str, locus: str | None = None) -> None:
        n = self.counts.get(category, 0) + 1
        self.counts[category] = n
        if n <= self.limit:
            self.items.append(ParseWarning(line_no, category, message, locus))

    def finish(self) -> list[ParseWarning]:
        out = list(self.items)
        for cat, n in self.counts.items():
            if n > self.limit:
                out.append(
                    ParseWarning(0, cat, f"{n - self.limit} further '{cat}' warnings suppressed")
                )
        return out


# --------------------------------------------------------------------------- #
# line-level parsing
# --------------------------------------------------------------------------- #

_HEADER_RE = re.compile(
    r"^#=IVTFF\s+(?P<alpha>\S+)"
    r"(?:\s+(?P<ver>\d+\.\d+(?:\.\d+)?))?"
    r"(?:\s+(?P<orig>[MDA])(?=\s|$))?"
    r"\s*(?P<extra>.*)$"
)
_LOCUS_ID_RE = re.compile(
    r"^(?P<page>[^\s.,;<>]+)\.(?P<num>\d+)(?P<suffix>[A-Za-z]?),"
    r"(?P<locator>\S)(?P<gen>[A-Z])(?P<sub>[a-z0-9])"
    r"(?:;(?P<tr>\S))?$"
)
_ALT_LOCUS_RE = re.compile(r"^[A-Z]{2}\d{3}(?:;\S)?$")
_PAGE_VAR_RE = re.compile(r"\$(\w)=(\S+?)(?=\s|$)")
_ANGLE_RE = re.compile(r"<[^>]*>")
_DRAWING_ONLY_RE = re.compile(r"<[-~]>")
_BALANCED_ALT_RE = re.compile(r"\[[^\]]*\]")
_BALANCED_LIG_RE = re.compile(r"\{[^}]*\}")


def _parse_header(line: str) -> Header:
    m = _HEADER_RE.match(line)
    if not m:
        raise IVTFFError(f"unreadable file header: {line!r}")
    return Header(
        alphabet=m.group("alpha"),
        version=m.group("ver"),
        origin=m.group("orig"),
        extra=(m.group("extra") or "").strip(),
        raw=line.rstrip(),
    )


def _split_markup(
    raw: str, warn, locus_label: str
) -> tuple[str, list[str], bool, bool, list[tuple[str, str]]]:
    """Strip comments, tags and paragraph marks from the text of a locus.

    Returns ``(text, comments, paragraph_start, paragraph_end, tags)``.
    Drawing breaks ``<->`` / ``<~>`` are kept in ``text``.
    """
    pieces: list[str] = []
    comments: list[str] = []
    tags: list[tuple[str, str]] = []
    para_start = para_end = False
    pos = 0
    for m in _ANGLE_RE.finditer(raw):
        pieces.append(raw[pos : m.start()])
        pos = m.end()
        token = m.group()
        inner = token[1:-1]
        if inner.startswith("!"):
            comments.append(inner[1:])
        elif inner == "%":
            para_start = True
        elif inner == "$":
            para_end = True
        elif inner in ("-", "~"):
            pieces.append(token)
        elif inner.startswith("@") and "=" in inner:
            name, _, value = inner[1:].partition("=")
            tags.append((name, value))
        else:
            warn("unknown-markup", f"unknown in-line construct {token!r} removed", locus_label)
    pieces.append(raw[pos:])
    text = "".join(pieces)
    leftover = _DRAWING_ONLY_RE.sub("", text)
    if "<" in leftover or ">" in leftover:
        warn("stray-angle", "stray '<' or '>' removed (unclosed in-line comment?)", locus_label)
        kept: list[str] = []
        i = 0
        while i < len(text):
            if text[i : i + 3] in ("<->", "<~>"):
                kept.append(text[i : i + 3])
                i += 3
            elif text[i] in "<>":
                i += 1
            else:
                kept.append(text[i])
                i += 1
        text = "".join(kept)
    return text, comments, para_start, para_end, tags


def _validate_text(text: str, warn, locus_label: str) -> None:
    """Emit warnings for unbalanced brackets and malformed high-ASCII codes."""
    rest = _BALANCED_ALT_RE.sub("", text)
    rest = _BALANCED_LIG_RE.sub("", rest)
    if any(ch in rest for ch in "[]{}"):
        warn("unbalanced-bracket", "unbalanced '[' ']' '{' or '}'", locus_label)
    for m in _BALANCED_ALT_RE.finditer(text):
        if ":" not in m.group():
            warn("alternative-without-colon", f"uncertain reading {m.group()!r} has no ':'", locus_label)
    stripped = _HIGH_ASCII_RE.sub("", text)
    if "@" in stripped:
        warn("bad-high-ascii", "'@' not followed by three digits and ';'", locus_label)
    for m in _HIGH_ASCII_RE.finditer(text):
        if not 128 <= int(m.group(1)) <= 255:
            warn("bad-high-ascii", f"high-ASCII code {m.group(1)} outside 128..255", locus_label)


def _use_fillers(strip_fillers: bool | str, header: Header) -> bool:
    if strip_fillers == "auto":
        return header.is_eva_like and header.version_tuple < (2, 0)
    return bool(strip_fillers)


# --------------------------------------------------------------------------- #
# main parser
# --------------------------------------------------------------------------- #


def parse_ivtff_lines(
    lines: Iterable[str],
    *,
    strip_fillers: bool | str = "auto",
    path: str | None = None,
) -> Transliteration:
    """Parse IVTFF content given as an iterable of lines (see :func:`parse_ivtff`)."""
    warner = _Warner()

    def warn(category: str, message: str, locus: str | None = None, line_no: int | None = None) -> None:
        warner.add(cur_line if line_no is None else line_no, category, message, locus)

    cur_line = 0
    header: Header | None = None
    pages: list[Page] = []
    loci: list[Locus] = []
    comment_lines: list[tuple[int, str]] = []
    current_page: Page | None = None
    tag_state: dict[str, str | None] = {}
    fillers = False
    # pending continuation: [line_no, ident_info(dict), [text parts]]
    pending: list | None = None

    def finish_locus(p: list) -> None:
        nonlocal tag_state
        line_no, info, parts = p
        raw_text = "".join(parts)
        label = f"{info['page']}.{info['num']}{info['suffix']}"

        def w(category: str, message: str, locus: str | None = None) -> None:
            warner.add(line_no, category, message, locus)

        text, comments, ps, pe, tags = _split_markup(raw_text.strip(), w, label)
        text = re.sub(r"\s+", "", text)
        if fillers:
            text = text.replace("!", "").replace("%", "")
        _validate_text(text, w, label)
        for name, value in tags:
            tag_state[name] = None if value == "@" else value
        loci.append(
            Locus(
                folio=info["page"],
                num=int(info["num"]),
                num_suffix=info["suffix"],
                locator=info["locator"],
                locus_type=info["gen"] + info["sub"],
                generic_type=info["gen"],
                subtype=info["sub"],
                transcriber=info["tr"],
                raw_text=raw_text.strip(),
                text=text,
                comments=comments,
                paragraph_start=ps,
                paragraph_end=pe,
                inline_tags=tags,
                line_no=line_no,
                page_index=pages[-1].index if pages else -1,
                variables=dict(tag_state),
                continued=len(parts) > 1,
            )
        )

    for line_no, line in enumerate(lines, 1):
        cur_line = line_no
        line = line.rstrip("\r\n")
        if line_no == 1:
            line = line.lstrip("﻿")
            if line.startswith("#=IVTFF"):
                header = _parse_header(line)
                fillers = _use_fillers(strip_fillers, header)
                continue
            warn("missing-header", "first line is not an '#=IVTFF' file header")
            header = Header(None, None, None)
            fillers = _use_fillers(strip_fillers, header)
        if not line.strip():
            continue
        first = line[0]
        if first == "#":
            if pending is not None:
                warn("comment-in-continuation", "comment line between continuation lines")
            comment_lines.append((line_no, line[1:].strip() if line[1:2] == " " else line[1:]))
            if current_page is not None:
                current_page.comments.append(line[1:].strip())
            continue
        if first == "/":
            if pending is None:
                warn("orphan-continuation", "line starts with '/' but previous line did not end with '/'")
                continue
            body = line[1:]
            if body.rstrip().endswith("/"):
                pending[2].append(body.rstrip()[:-1])
            else:
                pending[2].append(body)
                finish_locus(pending)
                pending = None
            continue
        if pending is not None:
            warn("missing-continuation", "line ended with '/' but next line is not a continuation")
            finish_locus(pending)
            pending = None
        if first != "<":
            warn("unexpected-line", f"line starts with {first!r}; ignored")
            continue
        close = line.find(">")
        if close < 0:
            raise IVTFFError(f"line {line_no}: malformed identifier (no closing '>'): {line[:60]!r}")
        ident = line[1:close]
        rest = line[close + 1 :]
        if re.search(r"\s", ident):
            raise IVTFFError(f"line {line_no}: whitespace inside identifier <{ident}>")
        if "." not in ident:
            if _ALT_LOCUS_RE.match(ident):
                warn("alternative-locus", f"alternative locus identifier <{ident}> is not supported; skipped")
                continue
            # page header
            if not ident:
                raise IVTFFError(f"line {line_no}: empty page name '<>'")
            variables: dict[str, str] = {}
            m = re.search(r"<!(.*?)>", rest)
            if m:
                variables = {k: v for k, v in _PAGE_VAR_RE.findall(m.group(1))}
            current_page = Page(ident, variables, line_no, len(pages))
            pages.append(current_page)
            tag_state = {k: (None if v == "@" else v) for k, v in variables.items()}
            continue
        m = _LOCUS_ID_RE.match(ident)
        if not m:
            raise IVTFFError(f"line {line_no}: malformed locus identifier <{ident}>")
        info = m.groupdict()
        label = ident
        if info["locator"] not in LOCATORS:
            warn("unknown-locator", f"locator {info['locator']!r} is not defined", label)
        if info["gen"] not in GENERIC_TYPES:
            warn("unknown-locus-type", f"generic locus type {info['gen']!r} is not defined", label)
        if current_page is None:
            current_page = Page(info["page"], {}, line_no, len(pages))
            pages.append(current_page)
            tag_state = {}
            warn("locus-before-page", "locus line before any page header; page created", label)
        elif current_page.name != info["page"]:
            warn("page-mismatch", f"locus belongs to {info['page']} but current page is {current_page.name}", label)
        body = rest
        if body.rstrip().endswith("/"):
            pending = [line_no, info, [body.rstrip()[:-1]]]
        else:
            finish_locus([line_no, info, [body]])
    if pending is not None:
        warn("missing-continuation", "file ends inside a continued line")
        finish_locus(pending)
    if header is None:
        raise IVTFFError("empty file")
    return Transliteration(
        header=header,
        pages=pages,
        loci=loci,
        warnings=warner.finish(),
        warning_counts=dict(warner.counts),
        comment_lines=comment_lines,
        path=path,
    )


def parse_ivtff(
    path: str | Path,
    encoding: str = "utf-8",
    errors: str = "replace",
    *,
    strip_fillers: bool | str = "auto",
) -> Transliteration:
    """Parse an IVTFF file.

    Parameters
    ----------
    path:
        Transliteration file.  It is treated as data and never executed.
    encoding, errors:
        Passed to :func:`open`.  The interlinear LSI file is not valid UTF-8 in
        some comments, hence the lenient default ``errors='replace'``.
    strip_fillers:
        ``'auto'`` (default) removes the interlinear filler characters ``!`` and
        ``%`` (and whitespace) for format versions below 2.0 with an Eva
        alphabet; ``True``/``False`` force the behaviour.

    Returns
    -------
    Transliteration
        With ``.header``, ``.pages``, ``.loci`` and ``.warnings``.

    Raises
    ------
    IVTFFError
        For a malformed locus identifier or an unreadable header.
    """
    p = Path(path)
    with p.open("r", encoding=encoding, errors=errors, newline=None) as fh:
        return parse_ivtff_lines(fh, strip_fillers=strip_fillers, path=str(p))


def parse_ivtff_text(text: str, *, strip_fillers: bool | str = "auto") -> Transliteration:
    """Parse IVTFF content given as one string (convenient for tests)."""
    return parse_ivtff_lines(text.splitlines(), strip_fillers=strip_fillers)


# --------------------------------------------------------------------------- #
# tokenizer
# --------------------------------------------------------------------------- #

_TOKEN_SCAN_RE = re.compile(r"<[-~]>|\[[^\]]*\]|\{[^}]*\}|@\d{3};|[.,]|\s+|.", re.S)
_SEP_RANK = {"U": 1, "S": 2, "D": 3}  # uncertain < definite < drawing break


def _expand(fragment: str) -> tuple[str, str, bool, bool]:
    """Expand ligature braces and high-ASCII codes in ``fragment``.

    Returns ``(normalized, raw, has_ligature, has_high_ascii)``.
    """
    has_lig = "{" in fragment and "}" in fragment
    norm = fragment.replace("{", "").replace("}", "") if has_lig else fragment
    has_hi = bool(_HIGH_ASCII_RE.search(norm))
    if has_hi:
        norm = decode_high_ascii(norm)
    return norm, fragment, has_lig, has_hi


def tokenize_text(
    text: str,
    space_policy: str = "split",
    alternatives: str = "first",
    *,
    warn=None,
) -> list[Token]:
    """Split cleaned locus text into tokens.

    Parameters
    ----------
    text:
        Text without comments, tags and paragraph marks (as in :attr:`Locus.text`).
    space_policy:
        ``'split'``  : ``.`` and ``,`` both split; tokens next to a ``,`` get the
        ``after_uncertain_space`` / ``before_uncertain_space`` flags.
        ``'join'``   : ``,`` is removed, the two parts are joined into one token.
        ``'strict'`` : only ``.`` splits; ``,`` stays inside the token as the literal
        character :data:`UNCERTAIN_SPACE_MARK` (``','``).  Downstream length and glyph
        functions ignore that character.
        In every policy a drawing break ``<->`` / ``<~>`` splits and sets
        ``after_drawing_break`` / ``before_drawing_break``.
    alternatives:
        ``'first'``: ``[a:b:c]`` is replaced by the most likely (first) option.
        ``'mark'`` : the token keeps ``[a:b:c]`` literally (so it becomes a distinct,
        easily filtered type).  ``has_alternative`` is set in both modes.
    warn:
        Optional callable ``warn(category, message)``.

    Notes
    -----
    A run of several separators between two tokens (possible in old files after
    filler removal) is reduced to the strongest one (drawing break > ``.`` > ``,``).
    Separators at the start or end of the text are ignored.
    """
    if space_policy not in SPACE_POLICIES:
        raise ValueError(f"space_policy must be one of {SPACE_POLICIES}, got {space_policy!r}")
    if alternatives not in ALTERNATIVE_POLICIES:
        raise ValueError(f"alternatives must be one of {ALTERNATIVE_POLICIES}, got {alternatives!r}")

    tokens: list[Token] = []
    norm: list[str] = []
    raw: list[str] = []
    flag_alt = flag_lig = flag_hi = False
    n_unc = 0
    started = False
    sep_before: str | None = None  # separator that preceded the token being built
    pending: str | None = None  # strongest separator seen since the last token ended

    def close_token() -> None:
        nonlocal norm, raw, flag_alt, flag_lig, flag_hi, n_unc, started, sep_before
        body = "".join(norm)
        if body:
            tokens.append(
                Token(
                    text=body,
                    raw="".join(raw),
                    index=len(tokens),
                    has_unreadable="?" in body,
                    has_alternative=flag_alt,
                    has_ligature=flag_lig,
                    has_high_ascii=flag_hi,
                    after_uncertain_space=sep_before == "U",
                    after_drawing_break=sep_before == "D",
                    n_uncertain_inside=n_unc,
                )
            )
        norm, raw = [], []
        flag_alt = flag_lig = flag_hi = False
        n_unc = 0
        started = False
        sep_before = None

    def start_token() -> None:
        nonlocal started, sep_before, pending
        if started:
            return
        started = True
        if pending is not None and tokens:
            sep_before = pending
            tokens[-1] = replace(
                tokens[-1],
                before_uncertain_space=pending == "U",
                before_drawing_break=pending == "D",
            )
        pending = None

    def separator(kind: str) -> None:
        nonlocal pending
        close_token()
        if pending is None or _SEP_RANK[kind] > _SEP_RANK[pending]:
            pending = kind

    for m in _TOKEN_SCAN_RE.finditer(text):
        s = m.group()
        c = s[0]
        if s in ("<->", "<~>"):
            separator("D")
        elif s == ".":
            separator("S")
        elif s == ",":
            if space_policy == "split":
                separator("U")
            elif space_policy == "strict":
                start_token()
                norm.append(UNCERTAIN_SPACE_MARK)
                raw.append(UNCERTAIN_SPACE_MARK)
                n_unc += 1
            # 'join': the comma simply disappears
        elif c == "[" and len(s) >= 2 and s[-1] == "]":
            start_token()
            options = s[1:-1].split(":")
            if len(options) == 1 and warn:
                warn("alternative-without-colon", f"uncertain reading {s!r} has no ':'")
            expanded = [_expand(o) for o in options]
            flag_alt = True
            if alternatives == "first":
                norm.append(expanded[0][0])
                flag_lig = flag_lig or expanded[0][2]
                flag_hi = flag_hi or expanded[0][3]
            else:
                norm.append("[" + ":".join(e[0] for e in expanded) + "]")
                flag_lig = flag_lig or any(e[2] for e in expanded)
                flag_hi = flag_hi or any(e[3] for e in expanded)
            raw.append(s)
        elif c == "{" and len(s) >= 2 and s[-1] == "}":
            start_token()
            inner, _, _, hi = _expand(s[1:-1])
            norm.append(inner)
            raw.append(s)
            flag_lig = True
            flag_hi = flag_hi or hi
        elif c == "@" and len(s) == 5 and s[-1] == ";":
            start_token()
            code = int(s[1:4])
            if 128 <= code <= 255:
                norm.append(chr(HIGH_ASCII_BASE + code))
                flag_hi = True
            else:
                norm.append(s)
            raw.append(s)
        elif c.isspace():
            continue
        elif c in "[]{}<>":
            if warn:
                warn("stray-bracket", f"stray {c!r} ignored")
        else:
            start_token()
            norm.append(s)
            raw.append(s)
    close_token()
    return tokens


def tokenize_locus(
    locus: Locus,
    space_policy: str = "split",
    alternatives: str = "first",
    *,
    warn=None,
) -> list[Token]:
    """Tokenize ``locus.text``; see :func:`tokenize_text` for the policies."""
    return tokenize_text(locus.text, space_policy, alternatives, warn=warn)
