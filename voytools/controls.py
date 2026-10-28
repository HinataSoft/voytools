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


# --------------------------------------------------------------------------- #
# experiment 013: human gibberish (Gaskell & Bowern 2022) and Naibbe ciphertext (Greshko 2025)
# --------------------------------------------------------------------------- #

_PUNCT = re.compile(r"[^\w\s]|[\d_]")
PAGE_MARK = "__________"


def gibberish_words(line: str) -> list[str]:
    """Tokens of one line as cleaned by Gaskell & Bowern (2022): lower case, punctuation and digits deleted.

    Punctuation is deleted *inside* tokens (``co-nallol`` -> ``conallol``), as in their ``stats.py``.

    >>> gibberish_words("Co-nallol ningelritaid, ab2usem * nron")
    ['conallol', 'ningelritaid', 'abusem', 'nron']
    """
    return [w for w in (_PUNCT.sub("", t.lower()) for t in line.split()) if w]


def parse_gibberish(text: str) -> list[list[list[str]]]:
    """One writer's transcription as pages of lines of tokens.

    Lines are the writer's own lines; a line ``__________`` separates pages (Gaskell & Bowern's page
    marker).  Lines that are empty after cleaning are dropped, and so are empty pages.

    >>> parse_gibberish("ab cd\\nef\\n__________\\n* gh\\n")
    [[['ab', 'cd'], ['ef']], [['gh']]]
    """
    pages: list[list[list[str]]] = [[]]
    for ln in text.split("\n"):
        if ln.strip() == PAGE_MARK:
            pages.append([])
            continue
        toks = gibberish_words(ln)
        if toks:
            pages[-1].append(toks)
    return [p for p in pages if p]


def read_gibberish(path: str | Path) -> dict[str, list[list[list[str]]]]:
    """All writers of the Gaskell & Bowern gibberish data: writer id -> pages -> lines -> tokens.

    ``path`` is the folder of ``Gibberish - <id>.txt`` files or the zip archive
    ``data/gibberish_transcriptions.zip`` of github.com/danielgaskell/voynich.  Writers are sorted by id.
    """
    import zipfile

    path = Path(path)
    texts: dict[str, str] = {}
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                base = Path(name).name
                if base.startswith("Gibberish - ") and base.endswith(".txt"):
                    raw = z.read(name)
                    try:
                        texts[base[12:-4]] = raw.decode("utf-8-sig")
                    except UnicodeDecodeError:
                        texts[base[12:-4]] = raw.decode("latin-1")
    else:
        for p in sorted(path.glob("Gibberish - *.txt")):
            texts[p.name[12:-4]] = read_text(p)
    return {k: parse_gibberish(texts[k].replace("\r\n", "\n").replace("\r", "\n")) for k in sorted(texts)}


def read_lines(path: str | Path) -> list[list[str]]:
    """Whitespace-tokenised non-empty lines of a plain text file (no case or character changes)."""
    return [ln.split() for ln in read_text(path).split("\n") if ln.split()]


def lines_table(
    pages: Sequence[tuple[str, str, str, Sequence[Sequence[str]]]],
    alphabet: str = "latin",
    glyph_preset=None,
    name: str = "lines",
) -> pd.DataFrame:
    """Token table of a text in its *own* line and page layout.

    ``pages``: (folio, section, hand, lines) per page; every page is one paragraph, ``language`` is
    empty.  The columns are those of :func:`voytools.nulls.as_table` (all statistics accept the table).
    """
    from .nulls import as_table

    all_lines = [list(ln) for _, _, _, lines in pages for ln in lines]
    out = as_table(all_lines, lines_per_page=max(1, len(all_lines)), alphabet=alphabet, glyph_preset=glyph_preset)
    folio, sec, hand, lip, nlines = [], [], [], [], []
    for f, s, h, lines in pages:
        for i, ln in enumerate(lines):
            k = len(ln)
            folio += [f] * k
            sec += [s] * k
            hand += [h] * k
            lip += [i] * k
            nlines += [len(lines)] * k
    lip_a, n_a = np.asarray(lip), np.asarray(nlines)
    out["folio"] = folio
    out["illustration"] = sec
    out["hand"] = hand
    out["language"] = "?"
    out["line_in_par"] = lip_a
    out["line_index"] = lip_a
    out["locus_num"] = lip_a + 1
    out["locus_id"] = [f"{f}.{i + 1}" for f, i in zip(folio, lip)]
    out["is_par_first_line"] = lip_a == 0
    out["is_par_last_line"] = lip_a == n_a - 1
    out["locator"] = np.where(lip_a == 0, "@", "+")
    out["par_index"] = 0
    out.attrs.update({"alphabet": alphabet, "glyph_preset": glyph_preset, "source": name})
    return out


def crop_template(template: pd.DataFrame, n_tokens: int) -> pd.DataFrame:
    """Evenly spaced whole pages of a token table (in table order) with at most ``n_tokens`` tokens.

    Pages are taken with a constant stride (as many as fit), so every part of the manuscript is
    represented; then whole lines are dropped from the end until the size is at most ``n_tokens``.
    """
    pages = list(pd.unique(template["folio"]))
    sizes = template.groupby("folio", sort=False).size().reindex(pages).to_numpy()
    share = min(1.0, n_tokens / sizes.sum())
    keep, acc = [], 0.0
    for p, n in zip(pages, sizes):
        acc += share
        if acc >= 1 - 1e-9:
            keep.append(p)
            acc -= 1
    out = template[template["folio"].isin(keep)]
    if len(out) > n_tokens:
        lens = out.groupby("locus_id", sort=False).size()
        cum = lens.cumsum()
        ok = set(cum.index[cum <= n_tokens])
        out = out[out["locus_id"].isin(ok)]
    out = out.reset_index(drop=True)
    out.attrs.update(template.attrs)
    return out


def layout_stream(divisions: Sequence[Division], line_lengths: Sequence[int]) -> tuple[list[list[str]], list[str]]:
    """Fill the template lines with the divisions in order as one running text (prose, lines broken mechanically).

    The label of a line is the division of its first token.  The text is cut where the lines end;
    raises ``ValueError`` if it is too short.
    """
    items = [it for _, its in divisions for it in its]
    owner = np.repeat([lab for lab, its in divisions for it in its], [len(it) for _, its in divisions for it in its])
    lines = _fill(items, line_lengths, "prose")
    if lines is None:
        raise ValueError(f"text too short for {sum(line_lengths)} tokens")
    starts = np.concatenate([[0], np.cumsum(line_lengths)[:-1]]).astype(np.int64)
    return lines, [str(owner[i]) for i in starts]


def stream_control(
    divisions: Sequence[Division],
    template: pd.DataFrame,
    name: str = "control",
    writer_as_hand: bool = False,
    alphabet: str = "latin",
    glyph_preset=None,
) -> pd.DataFrame:
    """Control laid out as one running text in the template (:func:`layout_stream`), divisions kept as they are.

    ``illustration`` = majority division of the page (as :func:`control_table`); with
    ``writer_as_hand`` also ``hand`` = majority division (a division is a writer).  ``alphabet`` and
    ``glyph_preset`` set the glyph units (``'Eva-'`` / ``'extended'`` for an Eva ciphertext).
    """
    from .nulls import as_table

    lines, labels = layout_stream(divisions, template_lines(template))
    out = as_table(lines, template, alphabet=alphabet, glyph_preset=glyph_preset)
    out.attrs.update({"alphabet": alphabet, "glyph_preset": glyph_preset, "source": name})
    lab = np.repeat(np.asarray(labels, dtype=object), [len(x) for x in lines])
    out["control_division"] = lab
    major = pd.Series(lab, index=out.index).groupby(out["folio"].to_numpy()).agg(lambda s: s.value_counts().index[0])
    out["illustration"] = out["folio"].map(major).to_numpy()
    if writer_as_hand:
        out["hand"] = out["illustration"].to_numpy()
    out.attrs["control"] = {"name": name, "mode": "prose", "strategy": "stream",
                            "sections": list(dict.fromkeys(labels)), "source_tokens": n_tokens(divisions)}
    return out


# --------------------------------------------------------------------------- #
# experiment 020: controls with their *own* paragraphs in the template's lines and pages
# --------------------------------------------------------------------------- #


def tess_paragraphs(text: str, level: int = 1, target: int | None = None, min_tail: float = 0.5) -> list[Division]:
    """Perseus/Tesserae ``.tess`` text as divisions whose items are the source's own paragraphs.

    The division label is the first ``level`` components of the citation (``1``: the book); with
    ``target=None`` every ``.tess`` line (a chapter) is one item.  With ``target`` the sentences of a chapter
    (:func:`sentences`) are grouped in order until a group has at least ``target`` tokens (a chapter's last
    group shorter than ``min_tail * target`` joins the previous one) – the weak "sentence boundaries every
    ~ target tokens" paragraphs of spec 020; groups never cross a chapter.

    >>> tess_paragraphs("<c. m. 1.1> Aa bb. Cc dd. Ee ff.\\n<c. m. 2.1> Gg.", target=4)
    [('1', [['aa', 'bb', 'cc', 'dd'], ['ee', 'ff']]), ('2', [['gg']])]
    """
    pairs = []
    for ln in text.split("\n"):
        m = re.match(r"\s*<([^>]*)>\s*(.*)$", ln)
        if not m:
            continue
        ref = m.group(1).split()[-1]
        label = ".".join(ref.split(".")[:level])
        if target is None:
            pairs.append((label, words(m.group(2))))
            continue
        groups: list[list[str]] = []
        cur: list[str] = []
        for sent in sentences(m.group(2)):
            cur.extend(words(sent))
            if len(cur) >= target:
                groups.append(cur)
                cur = []
        if cur:
            if groups and len(cur) < min_tail * target:
                groups[-1].extend(cur)
            else:
                groups.append(cur)
        pairs.extend((label, g) for g in groups)
    return _group(pairs)


def merge_short_items(divisions: Sequence[Division], min_item: int) -> list[Division]:
    """Items shorter than ``min_item`` tokens are joined to the following item of their division (the last to the previous).

    >>> merge_short_items([("d", [["a"], ["b", "c"], ["d"]])], 2)
    [('d', [['a', 'b', 'c', 'd']])]
    """
    out: list[Division] = []
    for label, items in divisions:
        merged: list[list[str]] = []
        carry: list[str] = []
        for it in items:
            carry = carry + list(it)
            if len(carry) >= min_item:
                merged.append(carry)
                carry = []
        if carry:
            if merged:
                merged[-1].extend(carry)
            else:
                merged.append(carry)
        out.append((label, merged))
    return out


def layout_paragraphs(
    items: Sequence[Sequence[str]], line_lengths: Sequence[int]
) -> tuple[list[list[str]], list[bool], list[int]]:
    """Fill the template lines with the items in order; an item's end closes the line (shorter) and the next item
    starts on the next line (the paragraph layout of a manuscript: lines and pages of the template, paragraphs
    of the source).

    Returns the lines, a flag per line (``True`` = first line of an item) and the item index of every line.
    Empty items are skipped.  Raises ``ValueError`` if the items run out before the lines do.

    >>> layout_paragraphs([list("abcde"), list("fg")], [3, 3, 3])
    ([['a', 'b', 'c'], ['d', 'e'], ['f', 'g']], [True, False, True], [0, 0, 1])
    """
    lines, firsts, owner = [], [], []
    k, pos = 0, 0
    while k < len(items) and not items[k]:
        k += 1
    for n in line_lengths:
        if k >= len(items):
            raise ValueError(f"items run out after {len(lines)} of {len(line_lengths)} lines")
        it = items[k]
        take = it[pos: pos + n]
        lines.append(list(take))
        firsts.append(pos == 0)
        owner.append(k)
        pos += len(take)
        if pos >= len(it):
            k, pos = k + 1, 0
            while k < len(items) and not items[k]:
                k += 1
    return lines, firsts, owner


def paragraph_control(
    divisions: Sequence[Division], template: pd.DataFrame, name: str = "control", min_item: int = 1
) -> pd.DataFrame:
    """Token table of a control laid out with :func:`layout_paragraphs`: the template's lines and pages (lines
    cut short where a paragraph ends), paragraphs = the source's items.

    Items shorter than ``min_item`` tokens are merged first (:func:`merge_short_items`).  Paragraph columns
    (``par_index``, ``line_in_par``, ``is_par_first_line``, ``is_par_last_line``) are recomputed per page as in
    the manuscript tables (a page always starts a paragraph index 0); ``par_continues`` marks the tokens of a
    page whose last paragraph goes on on the next page (the "no end mark" case of spec 020), ``control_item``
    is the source item of the token.  ``illustration`` = majority division of the page, ``control_division``
    the line's division, as in :func:`control_table`.  ``attrs['control']`` records the layout.
    """
    from .nulls import as_table

    divs = merge_short_items(divisions, min_item) if min_item > 1 else list(divisions)
    items = [it for _, its in divs for it in its]
    labels = [lab for lab, its in divs for _ in its]
    lens = template_lines(template)
    lines, firsts, owner = layout_paragraphs(items, lens)
    # crop the template to the realised line lengths (first len(line) tokens of every locus)
    tok_in_line = template.groupby("locus_id", sort=False).cumcount().to_numpy()
    want = np.repeat([len(x) for x in lines], lens)
    tpl = template[tok_in_line < want].reset_index(drop=True)
    out = as_table(lines, tpl, alphabet="latin", glyph_preset=None)
    out.attrs.update({"alphabet": "latin", "glyph_preset": None, "source": name})
    n_tok = np.array([len(x) for x in lines])
    rep = np.repeat(np.arange(len(lines)), n_tok)
    first_a = np.asarray(firsts)
    folio_line = template.groupby("locus_id", sort=False)["folio"].first().to_numpy()
    page_start = np.ones(len(lines), dtype=bool)
    page_start[1:] = folio_line[1:] != folio_line[:-1]
    par_start = first_a | page_start
    par_idx = np.zeros(len(lines), dtype=int)
    cur = -1
    for i in range(len(lines)):
        cur = 0 if page_start[i] else (cur + 1 if par_start[i] else cur)
        par_idx[i] = cur
    page_end = np.ones(len(lines), dtype=bool)
    page_end[:-1] = page_start[1:]
    par_end = page_end.copy()
    par_end[:-1] |= par_start[1:]
    lip = np.zeros(len(lines), dtype=int)
    for i in range(1, len(lines)):
        lip[i] = 0 if par_start[i] else lip[i - 1] + 1
    continues = np.zeros(len(lines), dtype=bool)
    continues[:-1] = page_start[1:] & ~first_a[1:]
    cont_page = pd.Series(continues).groupby(folio_line).transform("any").to_numpy()
    out["par_index"] = par_idx[rep]
    out["line_in_par"] = lip[rep]
    out["is_par_first_line"] = par_start[rep]
    out["is_par_last_line"] = par_end[rep]
    out["locator"] = np.where(par_start[rep], "@", "+")
    out["par_continues"] = cont_page[rep]
    out["control_item"] = np.asarray(owner)[rep]
    tidx = out.groupby("locus_id", sort=False).cumcount().to_numpy()
    out["token_index"] = tidx
    out["n_tokens"] = n_tok[rep]
    out["pos_from_end"] = n_tok[rep] - 1 - tidx
    out["rel_pos"] = np.where(n_tok[rep] > 1, tidx / np.maximum(n_tok[rep] - 1, 1), 0.0)
    out["is_first"] = tidx == 0
    out["is_last"] = tidx == n_tok[rep] - 1
    lab = np.asarray([labels[k] for k in owner], dtype=object)[rep]
    out["control_division"] = lab
    major = pd.Series(lab, index=out.index).groupby(out["folio"].to_numpy()).agg(lambda s: s.value_counts().index[0])
    out["illustration"] = out["folio"].map(major).to_numpy()
    out.attrs["control"] = {"name": name, "mode": "paragraphs", "strategy": "own paragraphs in template lines",
                            "sections": list(dict.fromkeys(labels)), "source_tokens": n_tokens(divisions),
                            "n_items_used": int(max(owner) + 1), "min_item": min_item,
                            "tokens_laid": int(n_tok.sum()), "template_tokens": int(sum(lens))}
    return out
