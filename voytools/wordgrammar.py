"""Word grammars learned from real tokens: glyph n-gram attestation and a slot grammar.

These grammars constrain the auto-copy generator (:func:`voytools.nulls.autocopy_generate`
with ``edit_model='trigram'`` / ``'slot'`` and ``new_token_source='grammar'``) and measure how
much of a text a grammar admits (:func:`grammaticality`).

Models
------
``'bigram'``
    Every glyph bigram of the token, with word-boundary marks, is attested in the reference
    words (the constraint of the original generator).
``'trigram'``
    Every glyph trigram of the token with one boundary mark on each side (``^ g1 g2``, ...,
    ``g(n-1) gn $``; a one-glyph token gives ``^ g1 $``) is attested in the reference words.
``'slot'``
    A slot grammar in the spirit of Zattera (2022), learned from the reference words:

    * the grammar is a sequence of slots, each holding one glyph; a glyph may occupy several
      slots (e.g. ``o`` early and in the middle of words);
    * a token is grammatical when its glyphs can be assigned to slots in strictly increasing
      order (each slot used at most once), except that runs of a *repeatable* glyph (one whose
      immediate repetition is common in the reference, e.g. ``ee``, ``ii`` in Eva) fill a single
      slot, up to the longest attested run; and when all its bigrams are attested.

    Learning: the initial slot order is the order of the glyphs' mean relative position in the
    token (``i / (n - 1)``; 0.5 for one-glyph tokens), one slot per glyph.  Extra slots are then
    added greedily: in each step the (glyph, insertion point) that makes the largest number of
    reference tokens grammatical is inserted, until ``target_coverage`` of the reference tokens
    (token-weighted) is grammatical or the best step gains less than ``min_gain`` of them.  On the
    learning text itself the coverage is therefore about ``target_coverage``, not exactly 1.

Boundary marks are internal sentinel strings that cannot occur as glyphs (v101 uses ``$`` as a
glyph).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from .glyphs import collapse_eva

__all__ = ["WordGrammar", "learn_word_grammar", "grammaticality", "MODELS"]

MODELS = ("bigram", "trigram", "slot")
BEGIN = "\x02^"
END = "$\x03"


def _collapse_runs(seq: Sequence[str], repeatable: set[str]) -> tuple[str, ...]:
    out: list[str] = []
    for x in seq:
        if out and out[-1] == x and x in repeatable:
            continue
        out.append(x)
    return tuple(out)


def _forward(template: list[str], w: Sequence[str]) -> list[int]:
    """Greedy left-to-right match; ``f[i]`` = template position after matching ``w[:i]``."""
    f = [0]
    p = 0
    for x in w:
        try:
            p = template.index(x, p) + 1
        except ValueError:
            break
        f.append(p)
    return f


def _backward(template: list[str], w: Sequence[str]) -> list[int]:
    """Greedy right-to-left match; ``b[i]`` = latest start so that ``w[i:]`` matches ``template[b[i]:]``."""
    m = len(template)
    b = [m] * (len(w) + 1)
    p = m
    for i in range(len(w) - 1, -1, -1):
        q = p - 1
        while q >= 0 and template[q] != w[i]:
            q -= 1
        if q < 0:
            for j in range(i, -1, -1):
                b[j] = -1
            return b
        p = q
        b[i] = p
    return b


class WordGrammar:
    """Bigram, trigram and slot grammar of a set of reference words (see the module docstring).

    Use :func:`learn_word_grammar` (or :meth:`fit`) to build one.  Attributes: ``preset``,
    ``bigrams`` / ``trigrams`` (sets of attested n-grams with boundary marks), ``slots`` (the
    learned slot sequence), ``repeatable`` (glyphs whose runs fill one slot), ``max_run``,
    ``coverage`` (token-weighted share of the reference words the slot grammar admits),
    ``lengths`` (glyph lengths of the reference tokens, for sampling).
    """

    def __init__(
        self,
        preset: Any,
        bigrams: set,
        trigrams: dict,
        slots: list[str],
        repeatable: set[str],
        max_run: dict[str, int],
        lengths: np.ndarray,
        coverage: float,
    ):
        self.preset = preset
        self.bigrams = bigrams
        self.trigram_next = trigrams  # {(a, b): {c: count}}
        self.trigrams = {(a, b, c) for (a, b), d in trigrams.items() for c in d}
        self.slots = slots
        self.repeatable = repeatable
        self.max_run = max_run
        self.lengths = lengths
        self.coverage = coverage
        # nxt[p][g] = first slot index >= p holding glyph g
        m = len(slots)
        nxt: list[dict[str, int]] = [dict() for _ in range(m + 1)]
        for p in range(m - 1, -1, -1):
            nxt[p] = dict(nxt[p + 1])
            nxt[p][slots[p]] = p
        self._nxt = nxt
        self._cache: dict[tuple[str, tuple[str, ...]], bool] = {}
        start: dict[str, int] = {}
        for (x, y), d in trigrams.items():
            if x == BEGIN:
                start[y] = start.get(y, 0) + sum(d.values())
        self._start = start

    # ------------------------------------------------------------------ learning

    @classmethod
    def fit(
        cls,
        tokens: Iterable[str],
        preset: Any = "basic",
        *,
        target_coverage: float = 0.99,
        min_gain: float = 2e-4,
        repeat_share: float = 0.05,
        max_slots: int = 400,
    ) -> "WordGrammar":
        """Learn the grammar from real tokens (strings, split into glyphs with ``preset``).

        ``repeat_share``: a glyph is repeatable when it is followed by itself in at least this
        share of its within-word successors (and at least 10 times).  ``target_coverage``,
        ``min_gain`` (as a share of the reference tokens) and ``max_slots`` stop the greedy slot
        insertion.
        """
        cnt = Counter(collapse_eva(t, preset) for t in tokens)
        cnt.pop((), None)
        if not cnt:
            raise ValueError("no reference tokens")
        types = list(cnt)
        counts = [cnt[t] for t in types]
        total = float(sum(counts))
        bigrams: set = set()
        tri: dict = defaultdict(dict)
        succ: dict[str, Counter] = defaultdict(Counter)
        pos_sum: dict[str, float] = defaultdict(float)
        pos_n: Counter = Counter()
        lengths: list[int] = []
        for g, c in zip(types, counts):
            lengths += [len(g)] * c
            padded = (BEGIN,) + g + (END,)
            for a, b in zip(padded, padded[1:]):
                bigrams.add((a, b))
            for a, b, x in zip(padded, padded[1:], padded[2:]):
                d = tri[(a, b)]
                d[x] = d.get(x, 0) + c
            for a, b in zip(g, g[1:]):
                succ[a][b] += c
            n = len(g)
            for i, x in enumerate(g):
                pos_sum[x] += c * (i / (n - 1) if n > 1 else 0.5)
                pos_n[x] += c
        repeatable = {a for a, s in succ.items() if s[a] >= 10 and s[a] >= repeat_share * sum(s.values())}
        max_run: dict[str, int] = {}
        for g, c in zip(types, counts):
            run = 1
            for a, b in zip(g, g[1:]):
                run = run + 1 if a == b else 1
                if a == b and a in repeatable:
                    max_run[a] = max(max_run.get(a, 1), run)
        seqs = [_collapse_runs(g, repeatable) for g in types]
        slots = sorted(pos_n, key=lambda x: (pos_sum[x] / pos_n[x], x))
        coverage = 0.0
        while True:
            acc = 0.0
            gains: dict[str, np.ndarray] = {}
            for w, c in zip(seqs, counts):
                f = _forward(slots, w)
                if len(f) == len(w) + 1:
                    acc += c
                    continue
                b = _backward(slots, w)
                seen: set = set()
                for i in range(min(len(f), len(w))):
                    lo, hi = f[i], b[i + 1]
                    if hi < lo:
                        continue
                    x = w[i]
                    arr = gains.get(x)
                    if arr is None:
                        arr = gains[x] = np.zeros(len(slots) + 1)
                    for p in range(lo, hi + 1):
                        if (x, p) not in seen:
                            seen.add((x, p))
                            arr[p] += c
            coverage = acc / total
            if coverage >= target_coverage or not gains or len(slots) >= max_slots:
                break
            best = max(((float(a.max()), x, int(a.argmax())) for x, a in sorted(gains.items())), key=lambda t: t[0])
            if best[0] < min_gain * total:
                break
            slots = slots[: best[2]] + [best[1]] + slots[best[2]:]
        return cls(preset, bigrams, dict(tri), slots, repeatable, max_run, np.asarray(lengths), coverage)

    # ------------------------------------------------------------------ checks

    def slot_state(self, p: int, prev: str | None, run: int, x: str) -> tuple[int, int] | None:
        """Advance the slot parser by glyph ``x`` from state (slot position, previous glyph, run).

        Returns the new ``(position, run)`` or ``None`` if ``x`` cannot follow.
        """
        if x == prev and x in self.repeatable:
            if run + 1 > self.max_run.get(x, 1):
                return None
            return p, run + 1
        q = self._nxt[p].get(x)
        if q is None:
            return None
        return q + 1, 1

    def accepts(self, seq: Sequence[str], model: str = "slot") -> bool:
        """True if the glyph sequence ``seq`` is admitted by ``model`` (``'bigram'|'trigram'|'slot'``)."""
        seq = tuple(seq)
        key = (model, seq)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        if model not in MODELS:
            raise ValueError(f"model must be one of {MODELS}")
        ok = len(seq) > 0
        padded = (BEGIN,) + seq + (END,)
        if ok and model in ("bigram", "slot"):
            ok = all(pair in self.bigrams for pair in zip(padded, padded[1:]))
        if ok and model == "trigram":
            ok = all(t in self.trigrams for t in zip(padded, padded[1:], padded[2:]))
        if ok and model == "slot":
            p, prev, run = 0, None, 0
            for x in seq:
                st = self.slot_state(p, prev, run, x)
                if st is None:
                    ok = False
                    break
                p, run = st
                prev = x
        if len(self._cache) < 500_000:
            self._cache[key] = ok
        return ok

    def accepts_token(self, token: str, model: str = "slot") -> bool:
        """:meth:`accepts` for a token string (split with the grammar's preset)."""
        return self.accepts(collapse_eva(token, self.preset), model)

    # ------------------------------------------------------------------ sampling

    def sample(self, rng: np.random.Generator, max_restarts: int = 100) -> tuple[str, ...] | None:
        """A new word: random walk over attested trigrams from ``^`` to ``$``, slot-grammatical.

        The length is drawn from the reference length distribution; each step draws the next
        glyph among attested trigram continuations (weighted by their counts) that keep the
        word slot-grammatical, and the last glyph must be attested before ``$``.  A dead end
        restarts the walk (with a new length); ``None`` after ``max_restarts`` failures.
        """
        for _ in range(max_restarts):
            n = int(self.lengths[int(rng.integers(len(self.lengths)))])
            seq: list[str] = []
            a, b = BEGIN, None
            p, prev, run = 0, None, 0
            ok = True
            for k in range(n):
                # first glyph: bigram continuations of the word start
                nxt = self._start if b is None else self.trigram_next.get((a, b), {})
                last = k == n - 1
                cands = []
                for x, c in nxt.items():
                    if x == END:
                        continue
                    if last:
                        tail = (b if b is not None else BEGIN, x)
                        if END not in self.trigram_next.get(tail, {}):
                            continue
                    st = self.slot_state(p, prev, run, x)
                    if st is None:
                        continue
                    cands.append((x, c, st))
                if not cands:
                    ok = False
                    break
                w = np.asarray([c for _, c, _ in cands], dtype=float)
                x, _, (p, run) = cands[int(rng.choice(len(cands), p=w / w.sum()))]
                prev = x
                seq.append(x)
                a, b = (b if b is not None else BEGIN), x
            if ok and seq and self.accepts(seq, "slot") and self.accepts(seq, "trigram"):
                return tuple(seq)
        return None


def learn_word_grammar(tokens: Iterable[str], preset: Any = "basic", **kw: Any) -> WordGrammar:
    """Learn a :class:`WordGrammar` from real tokens (see :meth:`WordGrammar.fit` for options)."""
    return WordGrammar.fit(tokens, preset, **kw)


def grammaticality(
    df_or_tokens: pd.DataFrame | Sequence[str],
    model: str = "slot",
    *,
    grammar: WordGrammar | None = None,
    reference: pd.DataFrame | Sequence[str] | None = None,
    preset: Any = "auto",
    weight: str = "token",
) -> float:
    """Share of tokens that a word grammar admits.

    Parameters
    ----------
    df_or_tokens:
        Token table (``token`` column) or a sequence of token strings.
    model:
        ``'bigram'``, ``'trigram'`` or ``'slot'`` (see the module docstring).
    grammar:
        A fitted :class:`WordGrammar`; default: learned on ``reference``, or on the input
        itself when ``reference`` is ``None`` (then the result is about ``target_coverage`` for
        ``'slot'`` and 1.0 for the n-gram models).
    preset:
        Glyph preset; ``'auto'`` takes ``attrs['glyph_preset']`` of the table (``'basic'`` for a
        list).  Ignored when ``grammar`` is given (its preset is used).
    weight:
        ``'token'`` (share of running tokens) or ``'type'`` (share of distinct types).
    """
    if weight not in ("token", "type"):
        raise ValueError("weight must be 'token' or 'type'")

    def tokens_of(x) -> list[str]:
        return x["token"].tolist() if isinstance(x, pd.DataFrame) else list(x)

    def preset_of(x) -> Any:
        if preset != "auto":
            return preset
        if isinstance(x, pd.DataFrame):
            return x.attrs.get("glyph_preset", "basic")
        return "basic"

    toks = tokens_of(df_or_tokens)
    if grammar is None:
        ref = df_or_tokens if reference is None else reference
        grammar = learn_word_grammar(tokens_of(ref), preset_of(ref))
    if not toks:
        return float("nan")
    cnt = Counter(toks)
    ok = {t: grammar.accepts_token(t, model) for t in cnt}
    if weight == "type":
        return sum(ok.values()) / len(cnt)
    return sum(c for t, c in cnt.items() if ok[t]) / len(toks)
