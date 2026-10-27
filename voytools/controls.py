"""Real-language control texts laid out like the Voynich paragraph text (experiment 009).

A control is a list of *divisions* (books, chapters, botanical orders ...), each a list of *items*
(passages, paragraphs, entries), each a list of word tokens.  The functions here

* read the source files (Project Gutenberg, Perseus ``.tess``, CSV) into divisions
  (:func:`read_text`, :func:`gutenberg_body`, :func:`words`, :func:`parse_tess`,
  :func:`parse_headed`, :func:`parse_linnaeus`, :func:`parse_culpeper`, :func:`blocks`);
* choose divisions and crop the text to the size of a template (:func:`select_divisions`);
* set the words into the line and page layout of a real token table (:func:`layout_lines`,
  :func:`control_table`), either as running prose (lines broken mechanically) or with every line
  starting a new item (entries of a catalogue).

The resulting table has the template's metadata (pages, paragraphs, line positions, language,
hand) except ``illustration``, which holds the control's division (majority per page), and the
glyph units are letters (``attrs['glyph_preset'] = None``).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

Division = tuple[str, list[list[str]]]

_WORD = re.compile(r"[^\W\d_]+")


# --------------------------------------------------------------------------- #
# reading and tokenisation
# --------------------------------------------------------------------------- #


def read_text(path: str | Path) -> str:
    """File contents as text: UTF-8 (with or without BOM), else Latin-1; line ends normalised to ``\\n``."""
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def gutenberg_body(text: str) -> str:
    """Text between the ``*** START OF ...`` and ``*** END OF ...`` lines of a Project Gutenberg file.

    Text without these markers is returned unchanged.
    """
    lines = text.split("\n")
    start = next((i + 1 for i, ln in enumerate(lines) if ln.startswith("*** START OF")), 0)
    end = next((i for i, ln in enumerate(lines) if ln.startswith("*** END OF")), len(lines))
    return "\n".join(lines[start:end])


def words(text: str) -> list[str]:
    """Lower-cased word tokens: maximal runs of letters; punctuation, digits and symbols are dropped.

    >>> words("Gallia est omnis divisa in partes tres, quarum unam (1.2) incolunt_Belgae.")
    ['gallia', 'est', 'omnis', 'divisa', 'in', 'partes', 'tres', 'quarum', 'unam', 'incolunt', 'belgae']
    >>> words("dell'acqua, Æthiopia")
    ['dell', 'acqua', 'æthiopia']
    """
    return _WORD.findall(text.lower())


def paragraphs(text: str) -> list[str]:
    """Blocks of text separated by blank lines (inner line breaks become spaces)."""
    out, cur = [], []
    for ln in text.split("\n"):
        if ln.strip():
            cur.append(ln.strip())
        elif cur:
            out.append(" ".join(cur))
            cur = []
    if cur:
        out.append(" ".join(cur))
    return out


def _group(pairs: Sequence[tuple[str, list[str]]]) -> list[Division]:
    """Consecutive (label, item) pairs -> divisions; empty items are dropped."""
    out: list[Division] = []
    for label, item in pairs:
        if not item:
            continue
        if out and out[-1][0] == label:
            out[-1][1].append(item)
        else:
            out.append((label, [item]))
    return out


def parse_tess(text: str, level: int = 1) -> list[Division]:
    """Perseus/Tesserae ``.tess`` text (``<cels. med. 1.pr> words ...`` per line) as divisions.

    The division label is the first ``level`` dot-separated components of the citation's last field
    (``level=1``: the book); every line is one item.

    >>> parse_tess("<plin. nat. 16.1> Pomiferae arbores.\\n<plin. nat. 17.1> Arbores, quae.")
    [('16', [['pomiferae', 'arbores']]), ('17', [['arbores', 'quae']])]
    """
    pairs = []
    for ln in text.split("\n"):
        m = re.match(r"\s*<([^>]*)>\s*(.*)$", ln)
        if not m:
            continue
        ref = m.group(1).split()[-1]
        pairs.append((".".join(ref.split(".")[:level]), words(m.group(2))))
    return _group(pairs)


def parse_headed(text: str, heading: str, label: str | None = None) -> list[Division]:
    """Split a text at lines matching the regular expression ``heading``; items are paragraphs.

    The division label is ``label`` formatted with the match (``'{0}'`` = the whole heading line,
    ``'{1}'`` = group 1 ...; default: the stripped heading line).  Text before the first heading is
    dropped; heading lines themselves are not part of the items.

    >>> parse_headed("x\\nLiber I\\nalpha beta.\\n\\ngamma\\nLiber II\\ndelta", r"Liber (\\w+)", "{1}")
    [('I', [['alpha', 'beta'], ['gamma']]), ('II', [['delta']])]
    """
    rx = re.compile(heading)
    pairs, cur_label, buf = [], None, []

    def flush():
        if cur_label is not None:
            pairs.extend((cur_label, words(p)) for p in paragraphs("\n".join(buf)))

    for ln in text.split("\n"):
        m = rx.fullmatch(ln.strip())
        if m:
            flush()
            groups = (m.group(0),) + m.groups()
            cur_label = label.format(*groups) if label else ln.strip()
            buf = []
        else:
            buf.append(ln)
    flush()
    return _group(pairs)


_STARS = re.compile(r"\s*\*(\s+\*){4,}\s*")


def parse_linnaeus(text: str) -> list[Division]:
    """Linnaeus, *Species Plantarum* (Project Gutenberg e-texts 20771/20772) as divisions = orders.

    The body runs from the first ``_Classis I._``-style line to the first row of asterisks after it
    (transcriber's appendices follow).  Division label = class numeral + order name
    (``'I MONOGYNIA'``); items are paragraphs: genus headings, species diagnoses, synonyms and
    ``Habitat`` lines, the bracketed species sidenotes (``[_Indica._]``) included as short items.
    The class title lines (``MONANDRIA``) are dropped.
    """
    lines = gutenberg_body(text).split("\n")
    cls_rx = re.compile(r"_Classis ([IVXL]+)\._")
    ord_rx = re.compile(r"_([A-Z]+GYNIA)\._")
    start = next(i for i, ln in enumerate(lines) if cls_rx.fullmatch(ln.strip()))
    cls, order, pairs, buf = None, None, [], []

    def flush():
        if order is not None:
            pairs.extend((f"{cls} {order}", words(p)) for p in paragraphs("\n".join(buf)))

    for ln in lines[start:]:
        s = ln.strip()
        if _STARS.fullmatch(ln):
            break
        m = cls_rx.fullmatch(s)
        if m:
            flush()
            cls, order, buf = m.group(1), None, []
            continue
        m = ord_rx.fullmatch(s)
        if m:
            flush()
            order, buf = m.group(1), []
            continue
        if order is None:  # class title (e.g. MONANDRIA) before the first order
            continue
        buf.append(ln)
    flush()
    return _group(pairs)


_HERB = re.compile(r"    [A-Z][A-Z ,'\-]+\.?")
_SENT = re.compile(r"(?<=[.;:?!])\s+")


def sentences(text: str) -> list[str]:
    """Pieces of ``text`` ending at ``. ; : ? !`` (a crude clause split for item layouts)."""
    return [x for x in _SENT.split(text) if x.strip()]


def _smooth_letters(letters: list[str], window: int) -> list[str]:
    """Running (lower) median of the letters over ``window`` entries, i.e. the letter of the surrounding run."""
    h = window // 2
    out = []
    for i in range(len(letters)):
        w = sorted(letters[max(0, i - h): i + h + 1])
        out.append(w[(len(w) - 1) // 2])
    return out


def parse_culpeper(
    text: str,
    first: str = "AMARA DULCIS.",
    stop: str = "DIRECTIONS FOR MAKING",
    item: str = "paragraph",
    window: int = 9,
) -> list[Division]:
    """Culpeper, *The Complete Herbal* (Project Gutenberg 49513): herb entries, divisions = alphabet letter.

    The herbal runs from the heading ``first`` to the line starting with ``stop``.  An entry starts
    at an indented upper-case heading (``    ASH TREE.``; a heading may continue on the next line).
    The herbal is alphabetical by the herb's key word, but headings such as ``WILD TANSY`` or
    ``ENGLISH TOBACCO`` start with a modifier, so the division is the running median (over
    ``window`` entries) of the headings' first letters (a leading ``THE`` ignored).  Items are the
    heading and the paragraphs of the entry (``item='paragraph'``: ``Descript.``, ``Place.``,
    ``Time.``, ``Government and virtues.`` ...) or their sentence pieces (``item='sentence'``,
    :func:`sentences`).
    """
    if item not in ("paragraph", "sentence"):
        raise ValueError("item must be 'paragraph' or 'sentence'")
    lines = gutenberg_body(text).split("\n")
    start = next(i for i, ln in enumerate(lines) if ln.strip() == first)
    entries: list[tuple[str, list[str]]] = []
    prev_head = False
    for ln in lines[start:]:
        if ln.startswith(stop):
            break
        if _HERB.fullmatch(ln):
            if prev_head:  # continuation of a two-line heading
                entries[-1][1].insert(-1, ln)
                continue
            head = ln.strip()
            head = head[4:] if head.startswith("THE ") else head
            entries.append((head[0], [ln, ""]))
            prev_head = True
            continue
        if ln.strip():
            prev_head = False
        if entries:
            entries[-1][1].append(ln)
    letters = _smooth_letters([e[0] for e in entries], window)
    pairs = []
    for letter, (_, buf) in zip(letters, entries):
        for par in paragraphs("\n".join(buf)):
            pieces = [par] if item == "paragraph" else sentences(par)
            pairs.extend((letter, words(x)) for x in pieces)
    return _group(pairs)


def blocks(tokens: Sequence[str], n_blocks: int, limit: int | None = None, item_len: int = 60) -> list[Division]:
    """Text without usable divisions: the first ``limit`` tokens cut into ``n_blocks`` equal blocks.

    Items are consecutive runs of ``item_len`` tokens.  Labels are ``'b1'``, ``'b2'`` ...
    """
    toks = list(tokens[: limit or len(tokens)])
    edges = np.linspace(0, len(toks), n_blocks + 1).round().astype(int)
    out = []
    for b in range(n_blocks):
        part = toks[edges[b]: edges[b + 1]]
        out.append((f"b{b + 1}", [part[i: i + item_len] for i in range(0, len(part), item_len)]))
    return out


# --------------------------------------------------------------------------- #
# layout
# --------------------------------------------------------------------------- #


def n_tokens(divisions: Sequence[Division]) -> int:
    """Total number of tokens of a list of divisions."""
    return sum(len(it) for _, items in divisions for it in items)


def section_quotas(line_lengths: Sequence[int], k: int) -> np.ndarray:
    """Section index (0..k-1) of every line: lines are assigned in order to k blocks of equal token mass."""
    lengths = np.asarray(line_lengths, dtype=float)
    before = np.concatenate([[0.0], np.cumsum(lengths)[:-1]])
    return np.minimum((before / lengths.sum() * k).astype(np.int64), k - 1)


def _fill(items: Sequence[Sequence[str]], line_lengths: Sequence[int], mode: str) -> list[list[str]] | None:
    """Lines of the given lengths filled from ``items``; ``None`` if the items run out.

    ``mode='prose'``: one running token stream, lines broken mechanically.  ``mode='items'``: every
    line starts with a new item; a shorter item is followed by the next items, a longer one is cut
    at the end of the line (its rest is dropped).
    """
    if mode == "prose":
        flat = [t for it in items for t in it]
        if len(flat) < sum(line_lengths):
            return None
        out, pos = [], 0
        for n in line_lengths:
            out.append(flat[pos: pos + n])
            pos += n
        return out
    if mode != "items":
        raise ValueError("mode must be 'prose' or 'items'")
    out, i = [], 0
    for n in line_lengths:
        line: list[str] = []
        while len(line) < n:
            if i >= len(items):
                return None
            line.extend(items[i][: n - len(line)])
            i += 1
        out.append(line)
    return out


def merge_small(divisions: Sequence[Division], fits) -> list[Division]:
    """Merge consecutive divisions until each satisfies ``fits(items)``; a short tail joins the last group.

    Labels of merged divisions are joined with ``'+'``.
    """
    out: list[Division] = []
    cur_labels: list[str] = []
    cur_items: list[list[str]] = []
    for label, items in divisions:
        cur_labels.append(label)
        cur_items.extend(items)
        if fits(cur_items):
            out.append(("+".join(cur_labels), cur_items))
            cur_labels, cur_items = [], []
    if cur_items and out:
        label, items = out.pop()
        out.append(("+".join([label, *cur_labels]), items + cur_items))
    return out


def select_divisions(
    divisions: Sequence[Division],
    line_lengths: Sequence[int],
    k: int,
    mode: str = "prose",
    keep: Sequence[str] | None = None,
) -> list[Division]:
    """Choose at most ``k`` divisions that can each fill an equal share of the template's lines.

    The template lines are split into ``k`` consecutive blocks of equal token mass
    (:func:`section_quotas`).  Divisions (restricted to the labels in ``keep`` if given) too short to
    fill the largest block in ``mode`` are merged with the following ones (:func:`merge_small`); out
    of the resulting groups ``k`` are taken evenly spaced in source order.  If fewer than ``k`` groups
    exist, ``k`` is lowered until it fits (each block is then larger).
    """
    lens = np.asarray(line_lengths)
    pool0 = [d for d in divisions if keep is None or d[0] in keep]
    for kk in range(k, 0, -1):
        sec = section_quotas(line_lengths, kk)
        biggest = max((lens[sec == j].tolist() for j in range(kk)), key=sum)
        pool = merge_small(pool0, lambda items, b=biggest: _fill(items, b, mode) is not None)
        if len(pool) >= kk:
            idx = np.linspace(0, len(pool) - 1, kk).round().astype(int)
            return [pool[i] for i in idx]
    return []


def layout_lines(
    divisions: Sequence[Division], line_lengths: Sequence[int], mode: str = "prose"
) -> tuple[list[list[str]], list[str]]:
    """Fill the template lines from the divisions in order; block j of lines comes from division j.

    Returns the lines and the division label of every line.  Each division contributes from its
    beginning (contiguous text); raises ``ValueError`` if one is too short for its block.
    """
    k = len(divisions)
    if k == 0:
        raise ValueError("no divisions")
    sec = section_quotas(line_lengths, k)
    lens = np.asarray(line_lengths)
    lines, labels = [], []
    for j, (label, items) in enumerate(divisions):
        part = _fill(items, lens[sec == j].tolist(), mode)
        if part is None:
            raise ValueError(f"division {label!r} is too short for {int(lens[sec == j].sum())} tokens ({mode})")
        lines.extend(part)
        labels.extend([label] * len(part))
    return lines, labels


def template_lines(template: pd.DataFrame) -> list[int]:
    """Number of tokens of every line (locus) of a token table, in table order."""
    return template.groupby("locus_id", sort=False).size().tolist()


def control_table(
    lines: Sequence[Sequence[str]], line_labels: Sequence[str], template: pd.DataFrame, name: str = "control"
) -> pd.DataFrame:
    """Token table of a control text in the layout of ``template`` (see :func:`voytools.nulls.as_table`).

    ``illustration`` becomes the majority division label of each page (sections are page-level as in
    the manuscript); the line's own label is kept in ``control_division``.  Glyph units are letters
    (``attrs['glyph_preset'] = None``, ``attrs['alphabet'] = 'latin'``), ``attrs['source'] = name``.
    """
    from .nulls import as_table

    out = as_table(lines, template, alphabet="latin", glyph_preset=None)
    out.attrs.update({"alphabet": "latin", "glyph_preset": None, "source": name})
    lab = np.repeat(np.asarray(line_labels, dtype=object), [len(x) for x in lines])
    out["control_division"] = lab
    major = pd.Series(lab, index=out.index).groupby(out["folio"].to_numpy()).agg(lambda s: s.value_counts().index[0])
    out["illustration"] = out["folio"].map(major).to_numpy()
    return out


def layout_sequential(
    divisions: Sequence[Division], line_lengths: Sequence[int], mode: str = "prose", k: int = 7
) -> tuple[list[list[str]], list[str]]:
    """Fill the template lines from the start of the source, all divisions in order (no balancing).

    Used when the source is too short to give every section its own contiguous block (e.g. an item
    layout that drops the rest of every long item).  The divisions are first merged into about ``k``
    groups of equal source size (:func:`merge_small`); a line's label is the group of the item it
    starts with.  Raises ``ValueError`` if the source runs out.
    """
    total = n_tokens(divisions)
    groups = merge_small(divisions, lambda items: sum(map(len, items)) >= total / k) or list(divisions)
    items = [it for _, its in groups for it in its]
    owner = [g for g, (_, its) in enumerate(groups) for _ in its]
    lines = _fill(items, line_lengths, mode)
    if lines is None:
        raise ValueError(f"source too short for {sum(line_lengths)} tokens ({mode})")
    # label of each line: the group of the item (prose: token) at the start of the line
    if mode == "prose":
        tok_owner = np.repeat(owner, [len(it) for it in items])
        starts = np.concatenate([[0], np.cumsum(line_lengths)[:-1]]).astype(np.int64)
        labels = [groups[int(tok_owner[i])][0] for i in starts]
    else:
        labels, i = [], 0
        for n in line_lengths:
            labels.append(groups[owner[i]][0])
            got = 0
            while got < n:
                got += min(len(items[i]), n - got)
                i += 1
    return lines, labels


def build_control(
    divisions: Sequence[Division],
    template: pd.DataFrame,
    k: int,
    mode: str = "prose",
    name: str = "control",
    keep: Sequence[str] | None = None,
    strategy: str = "auto",
    min_k: int = 4,
) -> pd.DataFrame:
    """Lay a control text out in the lines of ``template`` and return its token table.

    ``strategy='balanced'``: :func:`select_divisions` + :func:`layout_lines` (``k`` sections of equal
    size, each a contiguous passage from the start of a division); ``'sequential'``:
    :func:`layout_sequential` (the text from its start, sections where the source has them);
    ``'auto'``: balanced if it yields at least ``min_k`` sections, else sequential.
    ``attrs['control']`` records name, mode, strategy, the section labels and the source size.
    """
    lens = template_lines(template)
    chosen = select_divisions(divisions, lens, k, mode, keep) if strategy in ("balanced", "auto") else []
    if strategy == "balanced" or (strategy == "auto" and len(chosen) >= min_k):
        used = "balanced"
        lines, labels = layout_lines(chosen, lens, mode)
    elif strategy in ("sequential", "auto"):
        used = "sequential"
        pool = [d for d in divisions if keep is None or d[0] in keep]
        lines, labels = layout_sequential(pool, lens, mode, k)
    else:
        raise ValueError("strategy must be 'auto', 'balanced' or 'sequential'")
    out = control_table(lines, labels, template, name)
    out.attrs["control"] = {"name": name, "mode": mode, "strategy": used,
                            "sections": list(dict.fromkeys(labels)), "source_tokens": n_tokens(divisions)}
    return out
