"""Null models: token shuffles and a self-citation (auto-copy) text generator.

Shuffles
--------
:func:`shuffle_within_lines` and :func:`shuffle_within_pages` permute the tokens inside
each line (locus) or page.  They preserve the vocabulary and its frequencies, the line
lengths and all positional metadata, and destroy word order.  Everything that is a
property of the *token itself* (``token``, ``token_raw``, ``len_chars``, ``len_glyphs``,
the ``has_*`` quality flags and, if present, ``prefix`` / ``stem`` / ``suffix``) moves with
the token; positional columns stay where they are.  Consequently every statistic that
depends only on the bag of words or on the inside of words (vocabulary summary, length
distribution, affix inventories, decomposition, paradigm metrics, within-token glyph
statistics) is *identical* on a shuffled copy.  Shuffles are null models for order effects:
adjacent-token MI, edge-glyph MI, neighbour similarity, positional profiles and the glyph
entropy of the space-less text (whose glyph pairs at word boundaries depend on word order).

Auto-copy generator
-------------------
:func:`autocopy_generate` is a deliberately simple baseline in the spirit of Timm &
Schinner (2020), who proposed that Voynichese was produced by *self-citation*: each new
word is a slightly modified copy of one of the words written shortly before.  Here each
new token

1. with probability ``p_new`` is drawn from the seed vocabulary (frequency-weighted), otherwise
   copies a token chosen uniformly (or with geometric recency weight) from the previous
   ``window`` tokens;
2. then receives 0, 1 or 2 random glyph edits (substitution by a visually similar glyph from
   ``edit_table`` when that is legal, insertion of a glyph, deletion), with probability
   ``p_copy_exact`` of no edit and ``p_two_edits`` of a second edit among edited copies;
3. edits respect a glyph-bigram model estimated from reference words (``constrain_bigrams=True``):
   a new or substituted glyph is drawn from the glyphs that may follow its left neighbour in the
   reference words *and* may precede its right neighbour (word boundaries count as neighbours),
   with their reference frequencies; a deletion is only allowed if the two glyphs that become
   adjacent form an attested pair; words never exceed the longest reference word.  Without such
   a constraint pure copy-and-mutate dynamics drift into random glyph strings (words grow, the
   vocabulary explodes), because nothing keeps the words Voynich-like.  The constraint keeps the
   word-internal glyph pairs attested while leaving all *sequence* effects to the copying.
   If no legal edit is found in ``max_tries`` attempts the token stays an exact copy.

Grammar-bound edits (experiment 003): ``edit_model='trigram'`` additionally requires every
glyph trigram of the edited token (with boundary marks) to be attested in the reference words,
``edit_model='slot'`` requires the edited token to be admitted by a slot grammar learned from them
(see :mod:`voytools.wordgrammar`).  The candidate glyph of a substitution or insertion is drawn
only among the candidates that give an admitted token; a position/edit type without any such
candidate (or an inadmissible deletion) counts as a failed try, and after ``max_tries`` failures
the token stays an exact copy.  ``new_token_source='grammar'`` replaces the draw of a "new" token
from the seed vocabulary by a word sampled from the grammar (random walk over attested trigrams,
slot-grammatical, real length distribution; :meth:`voytools.wordgrammar.WordGrammar.sample`).
The defaults (``'bigram'``, ``'seed'``) reproduce the original generator exactly.

Lines are filled with token counts taken from a given list of real line lengths.  It is a
baseline for "how far do copy-and-modify dynamics alone get", not a reimplementation of the
published algorithm, whose rules are far more detailed.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from ._common import line_ids, page_ids
from .glyphs import collapse_eva, default_preset
from .table import TOKEN_COLUMNS
from .wordgrammar import WordGrammar, learn_word_grammar

__all__ = [
    "DEFAULT_EDIT_TABLE",
    "shuffle_within_lines",
    "shuffle_within_pages",
    "autocopy_generate",
    "as_table",
    "autocopy_like",
]

_TOKEN_PROPERTY_COLUMNS = (
    "token", "token_raw", "len_chars", "len_glyphs", "has_unreadable", "has_alternative",
    "has_ligature", "has_high_ascii", "prefix", "stem", "suffix",
)

DEFAULT_EDIT_TABLE: dict[str, tuple[str, ...]] = {
    # Illustrative confusion sets of similar-looking Eva glyphs (not fitted to data).
    "a": ("o", "e"), "o": ("a", "e", "y"), "e": ("o", "a", "ee"), "y": ("o", "a"),
    "l": ("r", "d"), "r": ("l", "s", "n"), "s": ("r", "l"), "d": ("l", "k"),
    "n": ("r", "m"), "m": ("n", "r"), "i": ("e", "ii"), "k": ("t", "d"), "t": ("k", "p"),
    "p": ("f", "t"), "f": ("p", "k"), "ch": ("sh", "ckh"), "sh": ("ch", "cth"),
    "q": ("k", "d"), "cth": ("ckh", "ch"), "ckh": ("cth", "cph"), "cph": ("cfh", "ckh"),
    "cfh": ("cph", "cth"), "g": ("l", "r"), "x": ("k", "t"),
}
"""Default glyph confusion table used for substitutions (Eva alphabet)."""


# --------------------------------------------------------------------------- #
# shuffles
# --------------------------------------------------------------------------- #


def _permute_within(df: pd.DataFrame, groups: np.ndarray, seed: int) -> pd.DataFrame:
    n = len(df)
    out = df.copy()
    if n == 0:
        return out
    rng = np.random.default_rng(seed)
    base = np.argsort(groups, kind="stable")
    order = np.lexsort((rng.random(n), groups))
    src = np.empty(n, dtype=np.int64)
    src[base] = order
    for c in _TOKEN_PROPERTY_COLUMNS:
        if c in df.columns:
            out[c] = df[c].to_numpy()[src]
    out.attrs.update(df.attrs)
    return out


def shuffle_within_lines(df: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """Copy of ``df`` with the tokens permuted inside every line (locus).

    Token properties travel with the token; position columns are unchanged (see the module
    docstring).  Deterministic for a given ``seed``.  The input must be in text order.
    """
    return _permute_within(df, line_ids(df), seed)


def shuffle_within_pages(df: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """Copy of ``df`` with the tokens permuted inside every page (across its lines)."""
    return _permute_within(df, page_ids(df), seed)


# --------------------------------------------------------------------------- #
# auto-copy generator
# --------------------------------------------------------------------------- #


def autocopy_generate(
    seed_tokens: Sequence[str],
    n_tokens: int,
    window: int = 8,
    edit_table: Mapping[str, Sequence[str]] | None = None,
    p_copy_exact: float = 0.35,
    *,
    p_two_edits: float = 0.25,
    edit_type_probs: tuple[float, float, float] = (0.6, 0.2, 0.2),
    p_new: float = 0.15,
    recency: float = 0.0,
    constrain_bigrams: bool = True,
    allowed_from: Sequence[str] | None = None,
    max_len: int | None = None,
    max_tries: int = 10,
    line_lengths: Sequence[int] | None = None,
    layout: str = "sample",
    preset: Any = "basic",
    seed: int = 0,
    edit_model: str = "bigram",
    new_token_source: str = "seed",
    word_grammar: WordGrammar | None = None,
    info: dict | None = None,
) -> list[list[str]]:
    """Generate ``n_tokens`` tokens by copying and modifying recent tokens.

    Parameters
    ----------
    seed_tokens:
        Real tokens that start the history and define the vocabulary for "new" tokens and the
        glyph frequencies for insertions.  Must not be empty.
    n_tokens:
        Number of tokens to generate (the seed is not part of the output).
    window:
        A copy source is chosen among the previous ``window`` tokens (seed tokens count as
        history at the start; the window runs across line breaks).
    edit_table:
        ``{glyph: candidates}`` used for substitutions; ``None`` means
        :data:`DEFAULT_EDIT_TABLE` for Eva presets and "any glyph of the seed" otherwise.
    p_copy_exact:
        Probability that the copy is left unchanged.
    p_two_edits:
        Among modified copies, the probability of two edits instead of one.
    edit_type_probs:
        Probabilities of ``(substitution, insertion, deletion)`` for each edit.
    p_new:
        Probability of taking a token from the seed vocabulary instead of copying.
    recency:
        0 = uniform choice in the window; ``r > 0`` weights the source at distance ``d``
        (1 = previous token) by ``exp(-r * (d - 1))``.
    constrain_bigrams, allowed_from, max_len, max_tries:
        Legality constraint for edits (see above).  ``allowed_from`` lists the reference words
        from which the glyph-bigram model is estimated (default: ``seed_tokens``); ``max_len`` is
        the maximum word length in glyphs (default: longest reference word).  The seed is usually
        short, so pass the full real vocabulary to ``allowed_from`` for a meaningful model.
    line_lengths:
        Real line lengths (tokens per line).  With ``layout='sample'`` lengths are drawn at random
        from this list; with ``layout='exact'`` they are used in order (cycled if needed).
        ``None`` draws ``1 + Poisson(7)``.
    preset:
        Glyph preset for splitting tokens into glyphs (``'basic'``, ``'extended'`` or ``None``).
    seed:
        Random seed (the output is deterministic for a given seed).
    edit_model:
        ``'bigram'`` (default, original behaviour), ``'trigram'`` or ``'slot'``: which word
        grammar an edited token must satisfy (see the module docstring).  ``'trigram'`` and
        ``'slot'`` require ``constrain_bigrams=True``.
    new_token_source:
        ``'seed'`` (default): a "new" token is drawn from ``seed_tokens``; ``'grammar'``: it is
        sampled from the slot/trigram grammar (falls back to a seed token if the sampler fails).
    word_grammar:
        A fitted :class:`~voytools.wordgrammar.WordGrammar` (same ``preset``) to use instead of
        learning one from the reference words; saves time when many samples share a reference.
    info:
        Optional dict that is filled with realised counts: ``n_tokens``, ``n_new``,
        ``n_edit_attempted`` (tokens drawn for editing), ``n_unchanged`` (output token equal to
        its source token, whether by ``p_copy_exact`` or by failed edits), ``n_edit_failed``
        (edit attempted but the token stayed unchanged), ``n_grammar_fallback`` and the shares
        ``p_unchanged`` = ``n_unchanged / n_tokens``, ``p_new_realised``.

    Returns
    -------
    list of lines, each a list of token strings; the last line may be shorter.
    """
    if not len(seed_tokens):
        raise ValueError("seed_tokens must not be empty")
    if layout not in ("sample", "exact"):
        raise ValueError("layout must be 'sample' or 'exact'")
    probs = np.asarray(edit_type_probs, dtype=float)
    if probs.shape != (3,) or probs.sum() <= 0:
        raise ValueError("edit_type_probs must be three non-negative numbers")
    probs = probs / probs.sum()
    if edit_model not in ("bigram", "trigram", "slot"):
        raise ValueError("edit_model must be 'bigram', 'trigram' or 'slot'")
    if new_token_source not in ("seed", "grammar"):
        raise ValueError("new_token_source must be 'seed' or 'grammar'")
    if edit_model != "bigram" and not constrain_bigrams:
        raise ValueError("edit_model 'trigram'/'slot' requires constrain_bigrams=True")
    rng = np.random.default_rng(seed)
    seeds = list(seed_tokens)
    glyph_cache: dict[str, tuple[str, ...]] = {}

    def glyphs(tok: str) -> tuple[str, ...]:
        g = glyph_cache.get(tok)
        if g is None:
            g = glyph_cache[tok] = collapse_eva(tok, preset)
        return g

    all_glyphs = [g for t in seeds for g in glyphs(t)]
    uniq_glyphs, g_counts = np.unique(np.asarray(all_glyphs, dtype=object), return_counts=True)
    g_p = g_counts / g_counts.sum()
    if edit_table is None:
        table: Mapping[str, Sequence[str]] = DEFAULT_EDIT_TABLE if preset is not None else {}
    else:
        table = edit_table

    ref_tokens = list(allowed_from) if allowed_from is not None else seeds
    pairs: dict[tuple[str, str], int] = {}
    successors: dict[str, dict[str, int]] = {}
    longest = 1
    for t in ref_tokens:
        gt = glyphs(t)
        longest = max(longest, len(gt))
        for a, b in zip(("^",) + gt, gt + ("$",)):
            pairs[(a, b)] = pairs.get((a, b), 0) + 1
            if b != "$":
                successors.setdefault(a, {})[b] = successors.get(a, {}).get(b, 0) + 1
    limit = max_len if max_len is not None else longest
    grammar = word_grammar
    if grammar is None and (edit_model != "bigram" or new_token_source == "grammar"):
        grammar = learn_word_grammar(ref_tokens, preset)
    admits = None
    if edit_model != "bigram":
        def admits(seq: Sequence[str]) -> bool:
            return grammar.accepts(seq, edit_model)

    def draw(
        left: str, right: str, exclude: str | None, prefer: Sequence[str] | None, build=None
    ) -> str | None:
        """A glyph that may stand between ``left`` and ``right`` (None if there is none).

        With a grammar-bound ``edit_model``, ``build(g)`` gives the edited token and only
        candidates whose token the grammar admits are kept.
        """
        if not constrain_bigrams:
            if prefer:
                return str(prefer[int(rng.integers(len(prefer)))])
            return str(rng.choice(uniq_glyphs, p=g_p))
        cands = [(g, c) for g, c in successors.get(left, {}).items() if g != exclude and (g, right) in pairs]
        if admits is not None:
            cands = [(g, c) for g, c in cands if admits(build(g))]
        if prefer:
            pref = [(g, c) for g, c in cands if g in prefer]
            cands = pref or cands
        if not cands:
            return None
        w = np.asarray([c for _, c in cands], dtype=float)
        return cands[int(rng.choice(len(cands), p=w / w.sum()))][0]

    def edit_once(seq: list[str]) -> list[str]:
        for _ in range(max_tries):
            kind = int(rng.choice(3, p=probs))
            if kind == 2 and len(seq) < 2:
                kind = 1
            if kind == 1:
                if len(seq) >= limit:
                    kind = 0
                else:
                    pos = int(rng.integers(len(seq) + 1))
                    left = seq[pos - 1] if pos > 0 else "^"
                    right = seq[pos] if pos < len(seq) else "$"
                    g = draw(left, right, None, None, lambda g, s=seq, i=pos: s[:i] + [g] + s[i:])
                    if g is not None:
                        return seq[:pos] + [g] + seq[pos:]
                    continue
            if kind == 2:
                pos = int(rng.integers(len(seq)))
                left = seq[pos - 1] if pos > 0 else "^"
                right = seq[pos + 1] if pos + 1 < len(seq) else "$"
                if not constrain_bigrams or (left, right) in pairs:
                    if admits is None or admits(seq[:pos] + seq[pos + 1 :]):
                        return seq[:pos] + seq[pos + 1 :]
                continue
            pos = int(rng.integers(len(seq)))
            left = seq[pos - 1] if pos > 0 else "^"
            right = seq[pos + 1] if pos + 1 < len(seq) else "$"
            g = draw(left, right, seq[pos], table.get(seq[pos]), lambda g, s=seq, i=pos: s[:i] + [g] + s[i + 1 :])
            if g is not None and g != seq[pos]:
                return seq[:pos] + [g] + seq[pos + 1 :]
        return seq

    history: list[str] = list(seeds[-max(window, 1):])
    out: list[str] = []
    w = np.exp(-recency * np.arange(window)) if recency > 0 else None
    n_new = n_attempted = n_unchanged = n_failed = n_fallback = 0
    for _ in range(n_tokens):
        if rng.random() < p_new:
            n_new += 1
            sampled = grammar.sample(rng) if new_token_source == "grammar" else None
            if sampled is not None:
                src = "".join(sampled)
                glyph_cache.setdefault(src, sampled)
            else:
                n_fallback += new_token_source == "grammar"
                src = seeds[int(rng.integers(len(seeds)))]
        else:
            recent = history[-window:][::-1]  # index 0 = previous token
            if w is None:
                src = recent[int(rng.integers(len(recent)))]
            else:
                ww = w[: len(recent)]
                src = recent[int(rng.choice(len(recent), p=ww / ww.sum()))]
        seq = list(glyphs(src))
        edited = False
        if rng.random() >= p_copy_exact:
            edited = True
            seq = edit_once(seq)
            if rng.random() < p_two_edits:
                seq = edit_once(seq)
        tok = "".join(seq)
        if not tok:
            tok = src
        n_attempted += edited
        if tok == src:
            n_unchanged += 1
            n_failed += edited
        out.append(tok)
        history.append(tok)
        if len(history) > 4 * max(window, 1):
            history = history[-window:]

    if info is not None:
        n = max(n_tokens, 1)
        info.update({
            "n_tokens": n_tokens, "n_new": n_new, "n_edit_attempted": n_attempted,
            "n_unchanged": n_unchanged, "n_edit_failed": n_failed, "n_grammar_fallback": n_fallback,
            "p_unchanged": n_unchanged / n, "p_new_realised": n_new / n,
        })

    # lay out into lines
    lines: list[list[str]] = []
    pos = 0
    k = 0
    lengths = [int(x) for x in line_lengths] if line_lengths is not None and len(line_lengths) else None
    while pos < len(out):
        if lengths is None:
            n = 1 + int(rng.poisson(7))
        elif layout == "exact":
            n = max(1, lengths[k % len(lengths)])
            k += 1
        else:
            n = max(1, lengths[int(rng.integers(len(lengths)))])
        lines.append(out[pos : pos + n])
        pos += n
    return lines


def as_table(
    lines: Sequence[Sequence[str]],
    template: pd.DataFrame | None = None,
    *,
    lines_per_page: int = 25,
    alphabet: str = "Eva-",
    glyph_preset: Any = "auto",
) -> pd.DataFrame:
    """Turn generated lines into a token table that all statistics functions accept.

    Without ``template`` the table has synthetic metadata: pages ``g001``, ``g002`` ...
    of ``lines_per_page`` lines, all loci of type ``P0`` forming one paragraph per page,
    no language/hand information.  With ``template`` (a real token table) the
    metadata and the line layout of the template are reused and only the token columns are
    replaced; ``lines`` must then contain exactly the template's line lengths in order
    (``autocopy_like`` guarantees this).

    The returned frame has ``attrs['alphabet']`` / ``attrs['glyph_preset']`` set.
    """
    preset = default_preset(alphabet) if glyph_preset == "auto" else glyph_preset
    flat = [t for line in lines for t in line]
    if template is not None:
        if len(flat) != len(template):
            raise ValueError(
                f"template has {len(template)} tokens but the generated lines have {len(flat)}; "
                "use autocopy_like() or layout='exact' with the template's line lengths"
            )
        out = template.copy()
        out["token"] = flat
        out["token_raw"] = flat
        out["len_chars"] = [len(t) for t in flat]
        out["len_glyphs"] = [len(collapse_eva(t, preset)) for t in flat]
        for c in ("has_unreadable", "has_alternative", "has_ligature", "has_high_ascii"):
            out[c] = False
        out = out.drop(columns=[c for c in ("prefix", "stem", "suffix") if c in out.columns])
        out.attrs.update(template.attrs)
        out.attrs["source"] = "autocopy"
        return out
    cols: dict[str, list] = {c: [] for c in TOKEN_COLUMNS}
    idx = 0
    for li, line in enumerate(lines):
        page = li // lines_per_page
        in_page = li % lines_per_page
        n = len(line)
        for ti, tok in enumerate(line):
            row = {
                "idx": idx, "folio": f"g{page + 1:03d}", "quire": None, "page_in_quire": None,
                "folio_in_quire": None, "bifolio": None, "illustration": None, "language": None,
                "hand": None, "currier_hand": None, "extraneous": None,
                "locus_id": f"g{page + 1:03d}.{in_page + 1}", "locus_num": in_page + 1,
                "locator": "@" if in_page == 0 else "+", "locus_type": "P0", "generic_type": "P",
                "subtype": "0", "transcriber": None, "par_index": 0, "line_in_par": in_page,
                "is_par_first_line": in_page == 0,
                "is_par_last_line": in_page == lines_per_page - 1 or li == len(lines) - 1,
                "line_index": in_page, "token_index": ti, "n_tokens": n, "pos_from_end": n - 1 - ti,
                "rel_pos": ti / (n - 1) if n > 1 else 0.0, "is_first": ti == 0, "is_last": ti == n - 1,
                "token": tok, "token_raw": tok, "len_chars": len(tok),
                "len_glyphs": len(collapse_eva(tok, preset)), "has_unreadable": False,
                "has_alternative": False, "has_ligature": False, "has_high_ascii": False,
                "after_uncertain_space": False, "before_uncertain_space": False,
                "after_drawing_break": False, "before_drawing_break": False,
            }
            for c in TOKEN_COLUMNS:
                cols[c].append(row[c])
            idx += 1
    df = pd.DataFrame(cols, columns=TOKEN_COLUMNS)
    df["par_index"] = pd.array(df["par_index"].tolist(), dtype="Int64")
    df["line_in_par"] = pd.array(df["line_in_par"].tolist(), dtype="Int64")
    df.attrs.update({"alphabet": alphabet, "glyph_preset": preset, "source": "autocopy"})
    return df


def autocopy_like(
    template: pd.DataFrame,
    seed_tokens: Sequence[str] | None = None,
    **params: Any,
) -> pd.DataFrame:
    """Auto-copy sample with exactly the layout and metadata of a real token table.

    Generates ``len(template)`` tokens with the template's line lengths (in order), seeded
    with ``seed_tokens`` (default: the first 100 tokens of the template), and returns
    :func:`as_table` of the result.  ``params`` go to :func:`autocopy_generate`
    (``window``, ``p_copy_exact``, ``seed``, ``edit_model`` ...).  The realised counts of
    ``info`` (see :func:`autocopy_generate`) are stored in ``attrs['autocopy_info']``.
    """
    if seed_tokens is None:
        seed_tokens = template["token"].tolist()[:100]
    lengths = np.bincount(line_ids(template)).tolist()
    params.setdefault("preset", template.attrs.get("glyph_preset", "basic"))
    info = params.pop("info", None)
    info = {} if info is None else info
    lines = autocopy_generate(
        seed_tokens, len(template), line_lengths=lengths, layout="exact", info=info, **params
    )
    out = as_table(lines, template)
    out.attrs["autocopy_info"] = dict(info)
    return out
