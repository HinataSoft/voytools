"""Statistics on token tables: frequency, entropy, mutual information, position, overlap, keyness.

All functions take the token DataFrame produced by :func:`voytools.table.build_tables`
(usually after :func:`voytools.table.filter_tokens`) and, where it makes sense, an
optional ``by=`` column name (or list of names) that splits the table into groups.
Functions with ``by`` always return a DataFrame; with ``by=None`` it has one row
labelled ``'all'``.  Missing values in a grouping column form the group ``'NA'``.

Method notes
------------
* Entropies are plug-in (maximum-likelihood) estimates in bits.  They are biased
  downwards for small samples; compare corpora of similar size.
* Adjacent-token and edge-glyph mutual information (MI) are *shuffle-corrected*: the
  observed plug-in MI is compared with the MI obtained when tokens are permuted
  inside each line, which keeps the marginal distributions and the line structure
  but removes order.  The excess (observed minus mean of the shuffles) estimates the
  genuine dependence; the raw plug-in MI of sparse data is mostly estimation bias
  (Rozanova & Temerev 2026 use the same device).
* The neighbour-similarity test follows Timm & Schinner (2020): in Voynichese a token
  is unusually often identical to, or one edit away from, one of the few preceding
  tokens.  The baseline permutes tokens within each page.
* Keyness is Dunning's log-likelihood ratio G2 (the standard corpus-linguistic
  "keyness" statistic) of a type in one group against all other groups.
* Positional profiles contrast token-initial/final glyphs of the first, last and
  middle tokens of a line (Currier 1976 noticed line-position effects; the
  Jensen-Shannon divergence summarises how different two distributions are).

References: Currier (1976) *Some important new statistical findings*; Reddy & Knight
(2011) *What we know about the Voynich manuscript*; Timm & Schinner (2020) *A possible
generating algorithm of the Voynich manuscript*, Cryptologia; Rozanova & Temerev (2026);
Zattera (2022) *A new transliteration alphabet brings new evidence of word structure and
multiple writing systems in the Voynich manuscript*, Cryptologia.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from ._common import (
    assemble,
    glyph_lists,
    iter_groups,
    line_ids,
    page_ids,
    resolve_preset,
    segments,
)

__all__ = [
    "token_frequency",
    "vocabulary_summary",
    "zipf_slope",
    "zipf_fit",
    "length_distribution",
    "ngram_counts",
    "entropy",
    "perplexity",
    "conditional_entropy",
    "glyph_entropy_report",
    "jsd",
    "levenshtein",
    "within_edit_distance",
    "adjacent_token_mi",
    "edge_glyph_mi",
    "PositionalProfile",
    "positional_profile",
    "neighbor_similarity",
    "vocabulary_overlap",
    "labels_vs_text",
    "keyness",
    "top_keyness",
]


# --------------------------------------------------------------------------- #
# frequency and vocabulary
# --------------------------------------------------------------------------- #


def _count_items(values: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """Distinct items and counts, ordered by decreasing count then alphabetically."""
    if len(values) == 0:
        return np.asarray([], dtype=object), np.zeros(0, dtype=np.int64)
    uniq, counts = np.unique(np.asarray(values, dtype=object), return_counts=True)
    order = np.argsort(-counts, kind="stable")
    return uniq[order], counts[order]


def token_frequency(df: pd.DataFrame, by: str | Sequence[str] | None = None, column: str = "token") -> pd.DataFrame:
    """Token (type) frequencies with ranks.

    Returns columns ``[by..., column, 'count', 'freq', 'rank']`` sorted by decreasing
    count (ties alphabetical); ``rank`` starts at 1 and is ordinal within each group.
    """
    frames = []
    for key, sub in iter_groups(df, by):
        items, counts = _count_items(sub[column].tolist())
        n = int(counts.sum())
        data: dict[str, Any] = {}
        if by is not None:
            cols = [by] if isinstance(by, str) else list(by)
            for c, k in zip(cols, key):
                data[c] = [k] * len(items)
        data[column] = items.tolist()
        data["count"] = counts.astype(np.int64)
        data["freq"] = counts / n if n else counts.astype(float)
        data["rank"] = np.arange(1, len(items) + 1)
        frames.append(pd.DataFrame(data))
    frames = [f for f in frames if len(f)]
    if not frames:
        return pd.DataFrame(columns=[column, "count", "freq", "rank"])
    return pd.concat(frames, ignore_index=True)


def zipf_fit(freq: Any, *, min_rank: int = 1, max_rank: int | None = None) -> dict[str, float]:
    """Log-log least-squares fit of the rank-frequency curve.

    ``freq`` may be a sequence of counts, a Series, or the output of
    :func:`token_frequency` (its ``count`` column is used).  Counts are sorted in
    decreasing order and regressed as ``log10(count) = intercept + slope * log10(rank)``
    over ranks ``min_rank..max_rank`` (default: all ranks).  Returns ``slope``,
    ``intercept``, ``r2`` and ``n_ranks``.
    """
    counts = _counts_array(freq)
    counts = np.sort(counts[counts > 0])[::-1]
    hi = len(counts) if max_rank is None else min(max_rank, len(counts))
    ranks = np.arange(1, len(counts) + 1)[min_rank - 1 : hi]
    y = np.log10(counts[min_rank - 1 : hi])
    if len(ranks) < 2:
        return {"slope": float("nan"), "intercept": float("nan"), "r2": float("nan"), "n_ranks": len(ranks)}
    x = np.log10(ranks)
    xm, ym = x.mean(), y.mean()
    sxx = float(((x - xm) ** 2).sum())
    slope = float(((x - xm) * (y - ym)).sum() / sxx)
    intercept = float(ym - slope * xm)
    ss_tot = float(((y - ym) ** 2).sum())
    ss_res = float(((y - (intercept + slope * x)) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"slope": slope, "intercept": intercept, "r2": r2, "n_ranks": int(len(ranks))}


def zipf_slope(freq: Any, *, min_rank: int = 1, max_rank: int | None = None) -> float:
    """Slope of the log-log rank-frequency regression (about -1 for Zipfian data)."""
    return zipf_fit(freq, min_rank=min_rank, max_rank=max_rank)["slope"]


def _counts_array(freq: Any) -> np.ndarray:
    if isinstance(freq, pd.DataFrame):
        freq = freq["count"]
    if isinstance(freq, Mapping):
        freq = list(freq.values())
    return np.asarray(list(freq) if not isinstance(freq, (np.ndarray, pd.Series)) else freq, dtype=float)


def vocabulary_summary(df: pd.DataFrame, by: str | Sequence[str] | None = None) -> pd.DataFrame:
    """Size and diversity of the vocabulary.

    Columns: ``n_tokens``, ``n_types``, ``n_hapax`` (types seen once), ``hapax_share``
    (hapax types / types), ``top50_coverage`` (share of tokens covered by the 50 most
    frequent types), ``ttr`` (types / tokens), ``mean_len`` (characters) and
    ``mean_len_glyphs``.
    """
    rows = []
    for key, sub in iter_groups(df, by):
        _, counts = _count_items(sub["token"].tolist())
        n = int(counts.sum())
        types = len(counts)
        hapax = int((counts == 1).sum())
        rows.append(
            (
                key,
                {
                    "n_tokens": n,
                    "n_types": types,
                    "n_hapax": hapax,
                    "hapax_share": hapax / types if types else float("nan"),
                    "top50_coverage": float(counts[:50].sum() / n) if n else float("nan"),
                    "ttr": types / n if n else float("nan"),
                    "mean_len": float(sub["len_chars"].mean()) if n else float("nan"),
                    "mean_len_glyphs": float(sub["len_glyphs"].mean()) if n else float("nan"),
                },
            )
        )
    return assemble(rows, by)


def length_distribution(
    df: pd.DataFrame, by: str | Sequence[str] | None = None, unit: str = "chars"
) -> pd.DataFrame:
    """Distribution of token lengths.

    Parameters
    ----------
    unit:
        ``'chars'`` (``len_chars``) or ``'glyphs'`` (``len_glyphs``, composites collapsed).

    Returns
    -------
    DataFrame with columns ``[by..., 'length', 'count', 'share']`` (long format).
    """
    if unit not in ("chars", "glyphs"):
        raise ValueError("unit must be 'chars' or 'glyphs'")
    col = "len_chars" if unit == "chars" else "len_glyphs"
    frames = []
    cols = [] if by is None else ([by] if isinstance(by, str) else list(by))
    for key, sub in iter_groups(df, by):
        vals = sub[col].to_numpy(dtype=np.int64)
        if len(vals) == 0:
            continue
        lengths, counts = np.unique(vals, return_counts=True)
        data: dict[str, Any] = {c: [k] * len(lengths) for c, k in zip(cols, key)}
        data["length"] = lengths
        data["count"] = counts
        data["share"] = counts / counts.sum()
        frames.append(pd.DataFrame(data))
    if not frames:
        return pd.DataFrame(columns=[*cols, "length", "count", "share"])
    return pd.concat(frames, ignore_index=True)


def ngram_counts(sequences: Any, n: int) -> pd.DataFrame:
    """Count n-grams in a token or glyph sequence.

    ``sequences`` is either a flat sequence of symbols (list of strings) or a list of
    sequences (e.g. one per line); n-grams never cross the boundary between two
    inner sequences.  A plain ``str`` is treated as a sequence of characters.

    Returns
    -------
    DataFrame with columns ``ngram`` (tuple), ``count`` and ``share``, sorted by
    decreasing count.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    counts: dict[tuple, int] = {}
    total = 0
    for seg in segments(sequences):
        seg = list(seg)
        for i in range(len(seg) - n + 1):
            g = tuple(seg[i : i + n])
            counts[g] = counts.get(g, 0) + 1
            total += 1
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return pd.DataFrame(
        {
            "ngram": [k for k, _ in items],
            "count": [v for _, v in items],
            "share": [v / total if total else 0.0 for _, v in items],
        }
    )


# --------------------------------------------------------------------------- #
# entropy
# --------------------------------------------------------------------------- #


def entropy(counts: Any) -> float:
    """Shannon entropy in bits of a distribution given as counts.

    ``counts`` may be a mapping (e.g. ``Counter``), a Series, an array or a list of
    non-negative numbers.  Zeros are ignored; an empty distribution has entropy 0.
    """
    if isinstance(counts, Mapping):
        arr = np.asarray(list(counts.values()), dtype=float)
    else:
        arr = np.asarray(counts, dtype=float).ravel()
    arr = arr[arr > 0]
    total = arr.sum()
    if total <= 0:
        return 0.0
    p = arr / total
    return float(-(p * np.log2(p)).sum())


def perplexity(entropy_bits: float) -> float:
    """Perplexity ``2**H`` for an entropy ``H`` in bits."""
    return float(2.0**entropy_bits)


def _encode(sequences: Any) -> tuple[np.ndarray, np.ndarray, int]:
    """Integer-encode symbols; returns ``(codes, segment_id, alphabet_size)``."""
    segs = segments(sequences)
    flat: list[Any] = []
    seg_id: list[int] = []
    for i, s in enumerate(segs):
        flat.extend(s)
        seg_id.extend([i] * len(s))
    if not flat:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), 0
    uniq = {}
    codes = np.empty(len(flat), dtype=np.int64)
    for i, sym in enumerate(flat):
        c = uniq.get(sym)
        if c is None:
            c = uniq[sym] = len(uniq)
        codes[i] = c
    return codes, np.asarray(seg_id, dtype=np.int64), len(uniq)


def _window_keys(codes: np.ndarray, seg: np.ndarray, size: int, order: int) -> np.ndarray:
    """Mixed-radix keys of all windows of length ``order`` lying inside one segment."""
    n = len(codes)
    if n < order:
        return np.zeros(0, dtype=np.int64)
    m = n - order + 1
    valid = seg[:m] == seg[order - 1 : order - 1 + m]
    key = np.zeros(m, dtype=np.int64)
    for j in range(order):
        key = key * size + codes[j : j + m]
    return key[valid]


def conditional_entropy(sequence: Any, order: int = 1) -> float:
    """Plug-in conditional entropy ``H(X_t | X_{t-order} .. X_{t-1})`` in bits.

    ``sequence`` is a flat sequence of hashable symbols or a list of sequences (lines);
    contexts never span two inner sequences.  The context distribution is taken from
    the same windows as the joint distribution (so ``H = H(joint) - H(context)``, as in
    Rozanova & Temerev's glyph-bigram entropy).

    Examples
    --------
    >>> round(conditional_entropy("aaab"), 4)
    0.9183
    >>> conditional_entropy("abababab")
    0.0
    """
    if order < 1:
        raise ValueError("order must be >= 1")
    codes, seg, size = _encode(sequence)
    if size == 0:
        return 0.0
    if (order + 1) * math.log2(max(size, 2)) < 62:
        joint = _window_keys(codes, seg, size, order + 1)
        if len(joint) == 0:
            return 0.0
        _, cj = np.unique(joint, return_counts=True)
        _, cc = np.unique(joint // size, return_counts=True)
        return max(0.0, entropy(cj) - entropy(cc))
    # fallback for huge alphabets: tuple counting
    joint_c: dict[tuple, int] = {}
    ctx_c: dict[tuple, int] = {}
    for s in segments(sequence):
        s = list(s)
        for i in range(len(s) - order):
            g = tuple(s[i : i + order + 1])
            joint_c[g] = joint_c.get(g, 0) + 1
            ctx_c[g[:-1]] = ctx_c.get(g[:-1], 0) + 1
    return max(0.0, entropy(joint_c) - entropy(ctx_c))


def glyph_entropy_report(
    df: pd.DataFrame,
    preset: Any = "auto",
    by: str | Sequence[str] | None = None,
    boundary: str | None = None,
    scope: str = "line",
    include_total: bool = True,
) -> pd.DataFrame:
    """Glyph-level entropy of the text with spaces removed.

    Columns: ``n_tokens``, ``n_glyphs`` (symbols in the stream, incl. boundary
    symbols), ``n_symbols`` (alphabet size), ``H1`` = H(X), ``H2`` = H(X_{t+1}|X_t),
    ``H3`` = H(X_{t+2}|X_t,X_{t+1}), ``perplexity`` = 2**H2 and ``mutual_info`` =
    H1 - H2.

    Parameters
    ----------
    preset:
        Glyph preset (``'basic'``, ``'extended'``, ``None``; ``'auto'`` follows the table).
    boundary:
        ``None``: spaces are removed (Rozanova & Temerev).  A string (e.g. ``'_'``) is
        inserted as an extra symbol between consecutive tokens.
    scope:
        ``'line'``: glyph bigrams never cross a line break; ``'all'``: the whole table is
        one stream (this is what Rozanova & Temerev did).
    include_total:
        With ``by`` given, add a final row ``all`` computed on the whole table.
    """
    from .glyphs import glyph_sequence

    preset = resolve_preset(df, preset)

    def one(sub: pd.DataFrame) -> dict[str, Any]:
        seqs = glyph_sequence(sub, preset, separator=boundary, scope=scope)
        flat = [g for s in seqs for g in s]
        _, counts = np.unique(np.asarray(flat, dtype=object), return_counts=True) if flat else (None, np.zeros(0))
        h1 = entropy(counts)
        h2 = conditional_entropy(seqs, 1) if flat else 0.0
        h3 = conditional_entropy(seqs, 2) if flat else 0.0
        return {
            "n_tokens": len(sub),
            "n_glyphs": len(flat),
            "n_symbols": int(len(counts)),
            "H1": h1,
            "H2": h2,
            "H3": h3,
            "perplexity": perplexity(h2),
            "mutual_info": h1 - h2,
        }

    rows = [(key, one(sub)) for key, sub in iter_groups(df, by)]
    res = assemble(rows, by)
    if by is not None and include_total:
        tot = one(df)
        cols = [by] if isinstance(by, str) else list(by)
        rec = {c: "all" for c in cols}
        rec.update(tot)
        res = pd.concat([res, pd.DataFrame([rec])], ignore_index=True)
    return res


# --------------------------------------------------------------------------- #
# divergence and edit distance
# --------------------------------------------------------------------------- #


def jsd(p: Any, q: Any) -> float:
    """Jensen-Shannon divergence in bits (0 = identical, 1 = disjoint support).

    ``p`` and ``q`` are non-negative vectors over the same support; they are normalised.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.shape != q.shape:
        raise ValueError("p and q must have the same length")
    sp, sq = p.sum(), q.sum()
    if sp <= 0 or sq <= 0:
        return float("nan")
    p, q = p / sp, q / sq
    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0
        return float((a[mask] * np.log2(a[mask] / b[mask])).sum())

    return max(0.0, 0.5 * kl(p, m) + 0.5 * kl(q, m))


def levenshtein(a: str, b: str) -> int:
    """Levenshtein edit distance (insertion, deletion, substitution; unit costs).

    >>> levenshtein("kitten", "sitting")
    3
    """
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def within_edit_distance(a: str, b: str, k: int) -> bool:
    """True if ``levenshtein(a, b) <= k`` (fast paths for ``k`` = 0 and 1)."""
    if a == b:
        return True
    la, lb = len(a), len(b)
    if abs(la - lb) > k:
        return False
    if k == 0:
        return False
    if k == 1:
        if la == lb:
            diff = 0
            for x, y in zip(a, b):
                if x != y:
                    diff += 1
                    if diff > 1:
                        return False
            return True
        if la < lb:
            a, b, la, lb = b, a, lb, la
        # a is one longer than b: b must equal a with one character removed
        i = 0
        while i < lb and a[i] == b[i]:
            i += 1
        return a[i + 1 :] == b[i:]
    return levenshtein(a, b) <= k


# --------------------------------------------------------------------------- #
# shuffle-corrected mutual information
# --------------------------------------------------------------------------- #


def _mi_from_pairs(x: np.ndarray, y: np.ndarray, kx: int, ky: int) -> float:
    n = len(x)
    if n == 0:
        return 0.0
    key = x.astype(np.int64) * ky + y
    uk, cnt = np.unique(key, return_counts=True)
    ux, uy = uk // ky, uk % ky
    nx = np.bincount(x, minlength=kx)
    ny = np.bincount(y, minlength=ky)
    mi = (cnt / n) * np.log2(cnt * n / (nx[ux].astype(float) * ny[uy]))
    return float(mi.sum())


def _shuffle_mi(
    left: np.ndarray,
    right: np.ndarray,
    line: np.ndarray,
    kl: int,
    kr: int,
    shuffles: int,
    seed: int,
) -> dict[str, float]:
    """Observed and within-line shuffled MI between ``left[i]`` and ``right[i+1]`` (same line)."""
    n = len(left)
    if n < 2:
        return {"n_pairs": 0, "mi_obs": float("nan"), "mi_shuffled_mean": float("nan"),
                "mi_shuffled_sd": float("nan"), "excess_bits": float("nan"), "z": float("nan")}
    same = line[:-1] == line[1:]
    obs = _mi_from_pairs(left[:-1][same], right[1:][same], kl, kr)
    rng = np.random.default_rng(seed)
    sims = np.empty(shuffles, dtype=float)
    for s in range(shuffles):
        order = np.lexsort((rng.random(n), line))
        sims[s] = _mi_from_pairs(left[order][:-1][same], right[order][1:][same], kl, kr)
    mean = float(sims.mean()) if shuffles else float("nan")
    sd = float(sims.std(ddof=1)) if shuffles > 1 else float("nan")
    return {
        "n_pairs": int(same.sum()),
        "mi_obs": obs,
        "mi_shuffled_mean": mean,
        "mi_shuffled_sd": sd,
        "excess_bits": obs - mean if shuffles else float("nan"),
        "z": (obs - mean) / sd if shuffles > 1 and sd > 0 else float("nan"),
    }


def adjacent_token_mi(
    df: pd.DataFrame,
    cap: int = 2000,
    shuffles: int = 100,
    seed: int = 0,
    by: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Shuffle-corrected mutual information between adjacent tokens within lines.

    Token types beyond the ``cap`` most frequent are mapped to ``'<other>'`` (this limits
    the sparsity bias).  Pairs ``(token_t, token_{t+1})`` are taken inside each locus.

    Columns: ``n_pairs``, ``n_types`` (capped alphabet size), ``H_cap`` (entropy of the capped
    token distribution), ``mi_obs``, ``mi_shuffled_mean``, ``mi_shuffled_sd``,
    ``excess_bits`` (observed minus shuffled mean), ``excess_share`` (``excess_bits / H_cap``)
    and ``z``.  The shuffles permute tokens within each line.
    """
    rows = []
    for key, sub in iter_groups(df, by):
        toks = sub["token"].tolist()
        items, counts = _count_items(toks)
        keep = {t: i for i, t in enumerate(items[:cap].tolist())}
        other = len(keep) if len(items) > cap else None
        codes = np.fromiter((keep.get(t, other) for t in toks), dtype=np.int64, count=len(toks))
        k = len(keep) + (1 if other is not None else 0)
        capped_counts = np.bincount(codes, minlength=k) if len(codes) else np.zeros(0)
        h_cap = entropy(capped_counts)
        res = _shuffle_mi(codes, codes, line_ids(sub), k, k, shuffles, seed) if len(codes) else {}
        res = dict(res)
        res.update(n_types=k, H_cap=h_cap)
        res["excess_share"] = res.get("excess_bits", float("nan")) / h_cap if h_cap > 0 else float("nan")
        order = ["n_pairs", "n_types", "H_cap", "mi_obs", "mi_shuffled_mean", "mi_shuffled_sd",
                 "excess_bits", "excess_share", "z"]
        rows.append((key, {c: res.get(c, float("nan")) for c in order}))
    return assemble(rows, by)


def edge_glyph_mi(
    df: pd.DataFrame,
    shuffles: int = 100,
    seed: int = 0,
    preset: Any = "auto",
    by: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Shuffle-corrected MI between the last glyph of a token and the first glyph of the next.

    Only pairs inside one line are used; the shuffles permute whole tokens within lines.
    Columns: ``n_pairs``, ``mi_obs``, ``mi_shuffled_mean``, ``mi_shuffled_sd``,
    ``excess_bits``, ``z``.  Because the glyphs at word edges are strongly constrained
    inside a word, a small positive excess is expected even if words are independent
    in order; compare with a null corpus.
    """
    preset = resolve_preset(df, preset)
    rows = []
    for key, sub in iter_groups(df, by):
        gl = glyph_lists(sub, preset)
        firsts = [g[0] if g else "" for g in gl]
        lasts = [g[-1] if g else "" for g in gl]
        fu = {g: i for i, g in enumerate(sorted(set(firsts)))}
        lu = {g: i for i, g in enumerate(sorted(set(lasts)))}
        f = np.fromiter((fu[g] for g in firsts), dtype=np.int64, count=len(firsts))
        la = np.fromiter((lu[g] for g in lasts), dtype=np.int64, count=len(lasts))
        res = _shuffle_mi(la, f, line_ids(sub), len(lu), len(fu), shuffles, seed) if len(f) else {}
        order = ["n_pairs", "mi_obs", "mi_shuffled_mean", "mi_shuffled_sd", "excess_bits", "z"]
        rows.append((key, {c: res.get(c, float("nan")) for c in order}))
    return assemble(rows, by)


# --------------------------------------------------------------------------- #
# positional profile
# --------------------------------------------------------------------------- #

_POS_CLASSES = ["line_first", "mid", "line_last", "parfirst_line_first", "parfirst_mid", "parfirst_line_last"]


@dataclass
class PositionalProfile:
    """Result of :func:`positional_profile`.

    Attributes
    ----------
    counts : DataFrame
        Glyph (rows) by position class (columns) counts.
    shares : DataFrame
        Column-normalised ``counts``.
    n : Series
        Number of tokens per class.
    jsd : DataFrame
        Pairwise Jensen-Shannon divergence (bits) between classes.
    pvalues : DataFrame or None
        Permutation p-values for ``jsd`` (only if ``permutations > 0``).
    """

    which: str
    unit: str
    counts: pd.DataFrame
    shares: pd.DataFrame
    n: pd.Series
    jsd: pd.DataFrame
    pvalues: pd.DataFrame | None = None

    def top(self, cls: str, k: int = 10) -> pd.Series:
        """The ``k`` most frequent glyphs of one class with their shares."""
        return self.shares[cls].sort_values(ascending=False, kind="stable").head(k)


def positional_profile(
    df: pd.DataFrame,
    which: str = "first",
    unit: str = "glyph",
    preset: Any = "auto",
    min_tokens: int = 2,
    permutations: int = 0,
    seed: int = 0,
) -> PositionalProfile:
    """Token-initial (or final) glyph distribution by position in the line.

    Position classes (tokens of lines with fewer than ``min_tokens`` tokens are ignored,
    so a one-word line is neither first nor last):

    ``line_first`` / ``mid`` / ``line_last``
        first, interior and last token of a line that is not the first line of a paragraph;
    ``parfirst_line_first`` / ``parfirst_mid`` / ``parfirst_line_last``
        the same three positions on the first line of a paragraph (``is_par_first_line``).

    Parameters
    ----------
    which:
        ``'first'`` glyph of each token or ``'last'`` glyph.
    unit:
        ``'glyph'`` (composites collapsed per ``preset``) or ``'char'``.
    permutations:
        If > 0, a permutation test is run for each pair of classes: class labels of the
        pooled glyphs are shuffled and JSD recomputed; ``pvalues`` holds
        ``(1 + #{JSD_perm >= JSD_obs}) / (1 + permutations)``.
    """
    if which not in ("first", "last"):
        raise ValueError("which must be 'first' or 'last'")
    if unit not in ("glyph", "char"):
        raise ValueError("unit must be 'glyph' or 'char'")
    preset = resolve_preset(df, preset)
    if unit == "glyph":
        gl = glyph_lists(df, preset)
        sym = [(g[0] if which == "first" else g[-1]) if g else "" for g in gl]
    else:
        sym = [(t[0] if which == "first" else t[-1]) if t else "" for t in df["token"].tolist()]
    n_tok = df["n_tokens"].to_numpy()
    is_first = df["is_first"].to_numpy(dtype=bool)
    is_last = df["is_last"].to_numpy(dtype=bool)
    parfirst = df["is_par_first_line"].to_numpy(dtype=bool)
    ok = n_tok >= max(2, min_tokens)
    cls = np.full(len(df), "", dtype=object)
    base = np.where(is_first, "line_first", np.where(is_last, "line_last", "mid"))
    cls[ok & ~parfirst] = base[ok & ~parfirst]
    cls[ok & parfirst] = np.char.add("parfirst_", base[ok & parfirst].astype(str)).astype(object)
    symbols = sorted(set(sym))
    sym_idx = {s: i for i, s in enumerate(symbols)}
    codes = np.fromiter((sym_idx[s] for s in sym), dtype=np.int64, count=len(sym))
    counts = np.zeros((len(symbols), len(_POS_CLASSES)), dtype=np.int64)
    per_class: dict[str, np.ndarray] = {}
    for j, c in enumerate(_POS_CLASSES):
        sel = codes[cls == c]
        per_class[c] = sel
        counts[:, j] = np.bincount(sel, minlength=len(symbols))
    cdf = pd.DataFrame(counts, index=symbols, columns=_POS_CLASSES)
    totals = cdf.sum(axis=0)
    shares = cdf / totals.where(totals > 0, np.nan)
    m = len(_POS_CLASSES)
    jm = np.full((m, m), np.nan)
    for i in range(m):
        for j in range(m):
            if totals.iloc[i] > 0 and totals.iloc[j] > 0:
                jm[i, j] = jsd(counts[:, i], counts[:, j])
    jdf = pd.DataFrame(jm, index=_POS_CLASSES, columns=_POS_CLASSES)
    pdf = None
    if permutations > 0:
        rng = np.random.default_rng(seed)
        pm = np.full((m, m), np.nan)
        for i in range(m):
            for j in range(i + 1, m):
                a, b = per_class[_POS_CLASSES[i]], per_class[_POS_CLASSES[j]]
                if len(a) == 0 or len(b) == 0:
                    continue
                pooled = np.concatenate([a, b])
                ge = 0
                for _ in range(permutations):
                    rng.shuffle(pooled)
                    ca = np.bincount(pooled[: len(a)], minlength=len(symbols))
                    cb = np.bincount(pooled[len(a) :], minlength=len(symbols))
                    if jsd(ca, cb) >= jm[i, j] - 1e-12:
                        ge += 1
                pm[i, j] = pm[j, i] = (1 + ge) / (1 + permutations)
        pdf = pd.DataFrame(pm, index=_POS_CLASSES, columns=_POS_CLASSES)
    return PositionalProfile(which, unit, cdf, shares, totals.rename("n"), jdf, pdf)


# --------------------------------------------------------------------------- #
# neighbour similarity (Timm & Schinner)
# --------------------------------------------------------------------------- #


def _neighbor_pair_keys(types: list[str], k: int) -> np.ndarray:
    """Sorted array of keys ``a*T + b`` for all ordered type pairs with edit distance <= k.

    Includes the identical pairs ``(a, a)``.  Candidate pairs are found with the
    symmetric-delete trick (words that share a string obtainable by at most ``k``
    deletions from each) and then verified exactly.
    """
    T = len(types)
    keys: list[int] = [i * T + i for i in range(T)]
    if k > 0:
        index: dict[str, list[int]] = {}
        for i, w in enumerate(types):
            variants = {w}
            frontier = {w}
            for _ in range(k):
                nxt = set()
                for v in frontier:
                    for p in range(len(v)):
                        nxt.add(v[:p] + v[p + 1 :])
                variants |= nxt
                frontier = nxt
            for v in variants:
                index.setdefault(v, []).append(i)
        seen: set[int] = set()
        for ids in index.values():
            if len(ids) < 2:
                continue
            for x in range(len(ids)):
                a = ids[x]
                for y in range(x + 1, len(ids)):
                    b = ids[y]
                    key = a * T + b if a < b else b * T + a
                    if key in seen:
                        continue
                    seen.add(key)
                    if within_edit_distance(types[a], types[b], k):
                        keys.append(a * T + b)
                        keys.append(b * T + a)
    return np.unique(np.asarray(keys, dtype=np.int64))


def _near_flags(codes: np.ndarray, page: np.ndarray, window: int, pair_keys: np.ndarray, T: int):
    """Per token: (has a near/identical token among the previous ``window`` tokens on its page)."""
    n = len(codes)
    near = np.zeros(n, dtype=bool)
    ident = np.zeros(n, dtype=bool)
    eligible = np.zeros(n, dtype=bool)
    for lag in range(1, window + 1):
        if lag >= n:
            break
        same = page[lag:] == page[:-lag]
        a = codes[:-lag]
        b = codes[lag:]
        eligible[lag:] |= same
        ident[lag:] |= same & (a == b)
        key = a * T + b
        pos = np.minimum(np.searchsorted(pair_keys, key), len(pair_keys) - 1)  # pair_keys is sorted and unique
        hit = pair_keys[pos] == key
        near[lag:] |= same & hit
    return near, ident, eligible


def neighbor_similarity(
    df: pd.DataFrame,
    window: int = 5,
    max_edit: int = 1,
    shuffles: int = 20,
    seed: int = 0,
    by: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Share of tokens that resemble one of the preceding tokens on the same page.

    A token counts as *near* if at least one of the previous ``window`` tokens of the same
    page (lines of the page are concatenated) is identical to it or within Levenshtein
    distance ``max_edit``.  The baseline permutes the tokens within each page.
    This is the "adjacent similarity" effect that Timm & Schinner (2020) model by
    self-citation.

    Columns: ``n_tokens`` (tokens with at least one predecessor on their page), ``window``,
    ``max_edit``, ``share_near``, ``share_near_shuffled`` (mean of the shuffles), ``share_near_sd``,
    ``ratio``, ``z``, plus the same for exact repeats only (``share_ident*``) and
    ``share_edit_only`` (near but not identical; observed and shuffled).
    """
    rows = []
    for key, sub in iter_groups(df, by):
        toks = sub["token"].tolist()
        n = len(toks)
        if n < 2:
            rows.append((key, {"n_tokens": 0}))
            continue
        types, codes = np.unique(np.asarray(toks, dtype=object), return_inverse=True)
        codes = codes.astype(np.int64)
        T = len(types)
        pair_keys = _neighbor_pair_keys(types.tolist(), max_edit)
        page = page_ids(sub)
        near, ident, elig = _near_flags(codes, page, window, pair_keys, T)
        n_el = int(elig.sum())
        obs_near = near.sum() / n_el
        obs_ident = ident.sum() / n_el
        rng = np.random.default_rng(seed)
        s_near = np.empty(shuffles)
        s_ident = np.empty(shuffles)
        for s in range(shuffles):
            order = np.lexsort((rng.random(n), page))
            nr, idn, el = _near_flags(codes[order], page, window, pair_keys, T)
            d = max(int(el.sum()), 1)
            s_near[s] = nr.sum() / d
            s_ident[s] = idn.sum() / d
        m_near = float(s_near.mean()) if shuffles else float("nan")
        sd_near = float(s_near.std(ddof=1)) if shuffles > 1 else float("nan")
        m_id = float(s_ident.mean()) if shuffles else float("nan")
        sd_id = float(s_ident.std(ddof=1)) if shuffles > 1 else float("nan")
        rows.append(
            (
                key,
                {
                    "n_tokens": n_el,
                    "window": window,
                    "max_edit": max_edit,
                    "share_near": float(obs_near),
                    "share_near_shuffled": m_near,
                    "share_near_sd": sd_near,
                    "ratio": float(obs_near) / m_near if m_near else float("nan"),
                    "z": (obs_near - m_near) / sd_near if sd_near and sd_near > 0 else float("nan"),
                    "share_ident": float(obs_ident),
                    "share_ident_shuffled": m_id,
                    "share_ident_sd": sd_id,
                    "ratio_ident": float(obs_ident) / m_id if m_id else float("nan"),
                    "share_edit_only": float(obs_near - obs_ident),
                    "share_edit_only_shuffled": m_near - m_id,
                },
            )
        )
    return assemble(rows, by)


# --------------------------------------------------------------------------- #
# overlap and keyness
# --------------------------------------------------------------------------- #


def vocabulary_overlap(df: pd.DataFrame, group_col: str, a: Any, b: Any, column: str = "token") -> dict[str, float]:
    """Overlap of the vocabularies of two groups (rows with ``df[group_col] == a`` / ``b``).

    Returns a dict with ``n_tokens_a/b``, ``n_types_a/b``, ``shared_types``, ``jaccard``
    (shared / union of types), ``a_token_in_b`` (share of a's tokens whose type occurs
    in b), ``a_type_in_b`` (share of a's types that occur in b) and the mirror values
    ``b_token_in_a``, ``b_type_in_a``.
    """
    col = df[group_col]
    ta = df.loc[(col == a).to_numpy(dtype=bool), column].tolist()
    tb = df.loc[(col == b).to_numpy(dtype=bool), column].tolist()
    sa, sb = set(ta), set(tb)
    shared = sa & sb
    union = sa | sb

    def share(tokens: list[str], other: set[str]) -> float:
        return sum(1 for t in tokens if t in other) / len(tokens) if tokens else float("nan")

    return {
        "n_tokens_a": len(ta),
        "n_tokens_b": len(tb),
        "n_types_a": len(sa),
        "n_types_b": len(sb),
        "shared_types": len(shared),
        "jaccard": len(shared) / len(union) if union else float("nan"),
        "a_token_in_b": share(ta, sb),
        "a_type_in_b": len(shared) / len(sa) if sa else float("nan"),
        "b_token_in_a": share(tb, sa),
        "b_type_in_a": len(shared) / len(sb) if sb else float("nan"),
    }


def labels_vs_text(
    df: pd.DataFrame,
    label_types: str | Iterable[str] = "L",
    text_types: str | Iterable[str] = "P",
) -> dict[str, float]:
    """:func:`vocabulary_overlap` of labels (generic type ``L``, group *a*) against paragraph text (``P``, group *b*)."""
    lab = [label_types] if isinstance(label_types, str) else list(label_types)
    txt = [text_types] if isinstance(text_types, str) else list(text_types)
    gt = df["generic_type"]
    kind = np.where(gt.isin(lab).to_numpy(dtype=bool), "labels", np.where(gt.isin(txt).to_numpy(dtype=bool), "text", ""))
    tmp = pd.DataFrame({"token": df["token"].tolist(), "kind": kind})
    return vocabulary_overlap(tmp, "kind", "labels", "text")


def _xlogx_ratio(a: np.ndarray, e: np.ndarray) -> np.ndarray:
    out = np.zeros_like(a, dtype=float)
    mask = (a > 0) & (e > 0)
    out[mask] = a[mask] * np.log(a[mask] / e[mask])
    return out


def keyness(
    df: pd.DataFrame,
    group_col: str,
    unit: str = "token",
    min_count: int = 5,
) -> pd.DataFrame:
    """Dunning log-likelihood keyness of every type in every group against the other groups.

    Parameters
    ----------
    group_col:
        Column that defines the groups (e.g. ``'language'``).  Rows with a missing value are ignored.
    unit:
        ``'token'`` (whole tokens) or ``'stem'`` (needs the ``stem`` column added by
        :func:`voytools.morphology.decompose_table`).
    min_count:
        Minimum total frequency of an item over all groups.

    Returns
    -------
    DataFrame with columns ``group, item, count_in, count_out, total_in, total_out,
    rel_in, rel_out, log2_ratio, G2, sign, signed_G2``, ordered by group and decreasing
    ``signed_G2``.  ``sign`` is +1 if the item is over-represented in the group, -1 if
    under-represented.  Under the null hypothesis G2 is approximately chi-square with 1
    degree of freedom (3.84 ~ p<0.05, 10.83 ~ p<0.001, ignoring multiple testing).
    """
    if unit not in ("token", "stem"):
        raise ValueError("unit must be 'token' or 'stem'")
    col = "token" if unit == "token" else "stem"
    if col not in df.columns:
        raise ValueError("unit='stem' needs a 'stem' column: run voytools.morphology.decompose_table first")
    valid = ~df[group_col].isna().to_numpy(dtype=bool)
    sub = df.loc[valid]
    groups = [str(g) for g in sub[group_col].tolist()]
    items = sub[col].tolist()
    g_names = sorted(set(groups))
    i_names, i_codes = np.unique(np.asarray(items, dtype=object), return_inverse=True)
    g_codes = np.searchsorted(np.asarray(g_names, dtype=object), np.asarray(groups, dtype=object))
    table = np.zeros((len(g_names), len(i_names)), dtype=np.int64)
    np.add.at(table, (g_codes, i_codes), 1)
    total_item = table.sum(axis=0)
    grand = int(table.sum())
    frames = []
    keep = total_item >= min_count
    for gi, g in enumerate(g_names):
        a = table[gi][keep].astype(float)
        b = (total_item[keep] - table[gi][keep]).astype(float)
        c = float(table[gi].sum())
        d = float(grand - c)
        if c == 0 or d == 0:
            continue
        e1 = c * (a + b) / (c + d)
        e2 = d * (a + b) / (c + d)
        g2 = 2.0 * (_xlogx_ratio(a, e1) + _xlogx_ratio(b, e2))
        rel_in, rel_out = a / c, b / d
        sign = np.where(rel_in >= rel_out, 1, -1)
        with np.errstate(divide="ignore", invalid="ignore"):
            lr = np.log2(rel_in / rel_out)
        frames.append(
            pd.DataFrame(
                {
                    "group": g,
                    "item": i_names[keep].tolist(),
                    "count_in": a.astype(np.int64),
                    "count_out": b.astype(np.int64),
                    "total_in": int(c),
                    "total_out": int(d),
                    "rel_in": rel_in,
                    "rel_out": rel_out,
                    "log2_ratio": lr,
                    "G2": g2,
                    "sign": sign,
                    "signed_G2": g2 * sign,
                }
            )
        )
    if not frames:
        return pd.DataFrame(
            columns=["group", "item", "count_in", "count_out", "total_in", "total_out", "rel_in",
                     "rel_out", "log2_ratio", "G2", "sign", "signed_G2"]
        )
    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["group", "signed_G2", "item"], ascending=[True, False, True], kind="stable")
    return out.reset_index(drop=True)


def top_keyness(result: pd.DataFrame, k: int = 15, direction: str = "positive") -> pd.DataFrame:
    """Top-``k`` items per group from :func:`keyness` (``direction``: ``'positive'``, ``'negative'`` or ``'both'``)."""
    if direction not in ("positive", "negative", "both"):
        raise ValueError("direction must be 'positive', 'negative' or 'both'")
    parts = []
    for g, sub in result.groupby("group", sort=True):
        pos = sub[sub["sign"] > 0].sort_values("G2", ascending=False, kind="stable").head(k)
        neg = sub[sub["sign"] < 0].sort_values("G2", ascending=False, kind="stable").head(k)
        if direction == "positive":
            parts.append(pos)
        elif direction == "negative":
            parts.append(neg)
        else:
            parts.extend([pos, neg])
    parts = [p for p in parts if len(p)]
    if not parts:
        return result.iloc[0:0].copy()
    return pd.concat(parts, ignore_index=True)
