"""Paragraph and page as units of content: stem sharing between lines (experiment 020).

The unit of analysis is the *line* (locus of paragraph text).  Two lines of the same page share a
unit (a stem or a whole token) when both contain it; the number of shared distinct units is
compared with its expectation under independent draws from a frequency table (the page's own unit
frequencies, so that the result is structure *above* the page vocabulary, or the global ones).

Workflow

* :func:`line_table` turns a token table (with a ``stem`` column from
  :func:`voytools.morphology.decompose_table`) into one row per line with the unit ids and the
  paragraph structure by position (``pos``, ``par``, ``lip`` ...).
* :func:`pair_content` computes, for every pair of lines of a page, the observed and expected
  sharing and the form statistics; :func:`pair_layout` adds distance and same-paragraph flags
  from the positions.  Layout nulls change only the positions (:func:`permute_lines`) or the
  paragraph boundaries (:func:`shuffle_boundaries`), so the content statistics are computed once.
* Metrics: :func:`sharing_by_distance`, :func:`boundary_step`, :func:`within_par_decay`,
  :func:`first_line_sharing`, :func:`page_boundary`, :func:`form_by_boundary`;
  :func:`all_metrics` returns them as one flat mapping (one ``runs.csv`` row).

Sharing ratio ``s = sum(observed) / sum(expected)`` over the pairs of a group; the expectation of a
pair of lines with ``n_i`` and ``n_j`` kept tokens is ``sum_s (1-(1-p_s)^n_i) (1-(1-p_s)^n_j)``
(draws with replacement from the frequency table ``p``).  All functions take an optional
``weights`` per page (bootstrap over pages).
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "GALLOWS",
    "strip_gallows",
    "unit_ids",
    "line_table",
    "paragraph_lengths",
    "structure_from_lengths",
    "permute_lines",
    "shuffle_boundaries",
    "pair_content",
    "pair_layout",
    "aggregate",
    "sharing_by_distance",
    "boundary_step",
    "within_par_decay",
    "first_line_sharing",
    "page_boundary",
    "form_by_boundary",
    "all_metrics",
]

# leading gallows glyphs (and their benched forms) that may be stripped from a stem (artefact check of spec point 3)
GALLOWS = {"eva": re.compile(r"^(?:cth|ckh|cph|cfh|[ktpf])+"), "v101": re.compile(r"^[gjfukhGJFUKH]+")}
LINE_COLS = ["folio", "rank", "pos", "par", "lip", "par_len", "n_pars", "n_lines", "language", "illustration", "hand"]


def strip_gallows(units: Sequence[str], alphabet: str = "eva") -> list[str]:
    """Units with leading gallows glyphs removed (``alphabet`` ``'eva'`` or ``'v101'``).

    >>> strip_gallows(["kchol", "cthey", "daiin", "pol"])
    ['chol', 'ey', 'daiin', 'ol']
    """
    rx = GALLOWS[alphabet]
    return [rx.sub("", str(u)) for u in units]


def unit_ids(units: Sequence[str]) -> tuple[np.ndarray, list[str]]:
    """Integer id of every unit (ids in order of first appearance) and the vocabulary."""
    codes, uniques = pd.factorize(pd.Series(list(units), dtype=object), sort=False)
    return codes.astype(np.int64), list(uniques)


def _alphabet(df: pd.DataFrame) -> str:
    a = str(df.attrs.get("alphabet", "")).lower()
    return "v101" if a.startswith("v101") else "eva"


def line_table(
    df: pd.DataFrame,
    unit: str = "stem",
    min_count: int = 2,
    drop_top: int = 0,
    gallows: bool = False,
    rank: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """One row per line (locus) of a paragraph-text token table, with unit ids and paragraph structure.

    ``unit``: column holding the unit (``'stem'`` or ``'token'``).  Units with fewer than
    ``min_count`` occurrences in ``df`` and the ``drop_top`` most frequent units are not kept (a
    line's ``n`` counts the kept tokens; lines without kept tokens stay, with ``n = 0``).
    ``gallows=True`` strips leading gallows glyphs from the units first (:func:`strip_gallows`).
    ``rank``: page -> order in the manuscript (default: order of appearance; used by
    :func:`page_boundary` to find consecutive pages).

    Columns: ``folio, rank, pos`` (line position on the page, 0-based), ``par, lip`` (paragraph
    index on the page and line in paragraph), ``par_len`` (lines of the paragraph), ``n_pars``,
    ``n_lines`` (of the page), ``language, illustration, hand``, ``n``, ``units`` (sorted distinct
    unit ids), ``counts`` (``{unit: {form: count}}``), ``locus_id``.  ``attrs``: ``unit_freq``
    (global counts of the kept units, indexed by id), ``form_freq`` (``{unit: {form: count}}``),
    ``vocab``.
    """
    d = df
    if unit not in d.columns:
        raise ValueError(f"column {unit!r} missing (run decompose_table first)")
    units = d[unit].astype(str).to_numpy()
    if gallows:
        units = np.asarray(strip_gallows(units, _alphabet(df)), dtype=object)
    forms = d["token"].astype(str).to_numpy()
    ids, vocab = unit_ids(units)
    counts = np.bincount(ids, minlength=len(vocab))
    keep = (counts >= min_count) & (np.array([v != "" for v in vocab]) if len(vocab) else np.array([], dtype=bool))
    if drop_top:
        top = np.argsort(-counts, kind="stable")[:drop_top]
        keep[top] = False
    kept_mask = keep[ids]
    form_ids, form_vocab = unit_ids(forms)

    folio = d["folio"].astype(str).to_numpy()
    locus = d["locus_id"].astype(str).to_numpy()
    # lines in table order
    first_idx = pd.Series(np.arange(len(d))).groupby(locus, sort=False).first()
    order = first_idx.sort_values().index.tolist()
    loc_pos = {l: k for k, l in enumerate(order)}
    line_of = np.array([loc_pos[l] for l in locus])

    meta = d.groupby(locus, sort=False).agg(
        folio=("folio", "first"), par=("par_index", "first"), lip=("line_in_par", "first"),
        language=("language", "first"), illustration=("illustration", "first"), hand=("hand", "first"),
    ).reindex(order)
    lines = meta.reset_index().rename(columns={"index": "locus_id", "locus_id": "locus_id"})
    lines["locus_id"] = order
    lines["par"] = pd.to_numeric(lines["par"], errors="coerce").fillna(0).astype(int)
    lines["lip"] = pd.to_numeric(lines["lip"], errors="coerce").fillna(0).astype(int)
    pages = list(dict.fromkeys(lines["folio"].tolist()))
    if rank is None:
        rank = {p: i for i, p in enumerate(pages)}
    lines["rank"] = lines["folio"].map(rank).astype(int)
    lines["pos"] = lines.groupby("folio", sort=False).cumcount()
    # paragraph index 0.. in order of position: a run of lines with the same par_index (the table's values may
    # skip, and a page may return to an earlier index after other paragraphs, e.g. ZL f84v; runs keep paragraphs contiguous)
    run = np.ones(len(lines), dtype=int)
    same_prev = (lines["par"].to_numpy()[1:] == lines["par"].to_numpy()[:-1]) & (lines["folio"].to_numpy()[1:] == lines["folio"].to_numpy()[:-1])
    run[1:] = ~same_prev
    lines["par"] = run
    lines["par"] = lines.groupby("folio", sort=False)["par"].cumsum().astype(int) - 1
    lines["lip"] = lines.groupby(["folio", "par"], sort=False).cumcount().astype(int)
    lines["par_len"] = lines.groupby(["folio", "par"], sort=False)["pos"].transform("size").astype(int)
    lines["n_pars"] = lines.groupby("folio", sort=False)["par"].transform("nunique").astype(int)
    lines["n_lines"] = lines.groupby("folio", sort=False)["pos"].transform("size").astype(int)

    n_lines = len(order)
    n = np.zeros(n_lines, dtype=int)
    units_col: list[np.ndarray] = [np.empty(0, dtype=np.int64) for _ in range(n_lines)]
    cnt_col: list[dict] = [dict() for _ in range(n_lines)]
    sel = np.flatnonzero(kept_mask)
    np.add.at(n, line_of[sel], 1)
    per_line: dict[int, dict] = {}
    for k in sel:
        li = int(line_of[k])
        dct = per_line.setdefault(li, {})
        inner = dct.setdefault(int(ids[k]), {})
        inner[int(form_ids[k])] = inner.get(int(form_ids[k]), 0) + 1
    for li, dct in per_line.items():
        units_col[li] = np.array(sorted(dct), dtype=np.int64)
        cnt_col[li] = dct
    lines["n"] = n
    lines["units"] = units_col
    lines["counts"] = cnt_col
    out = lines[LINE_COLS + ["n", "units", "counts", "locus_id"]].reset_index(drop=True)
    unit_freq = np.where(keep, counts, 0)
    form_freq: dict[int, dict[int, int]] = {}
    for u, f in zip(ids[sel], form_ids[sel]):
        form_freq.setdefault(int(u), {})
        form_freq[int(u)][int(f)] = form_freq[int(u)].get(int(f), 0) + 1
    out.attrs.update({"unit": unit, "min_count": min_count, "drop_top": drop_top, "gallows": gallows,
                      "stats": LineStats(unit_freq, form_freq, vocab, form_vocab),
                      "n_kept": int(kept_mask.sum()), "n_units": int(keep.sum())})
    return out


class LineStats:
    """Corpus-level tables of a :func:`line_table` (``attrs['stats']``): ``unit_freq`` (counts of the kept
    units by id), ``form_freq`` (``{unit: {form: count}}``), ``vocab``, ``form_vocab``.

    Kept out of the frame's plain ``attrs`` entries so that pandas does not deep-copy them on every operation.
    """

    def __init__(self, unit_freq, form_freq, vocab, form_vocab):
        self.unit_freq = np.asarray(unit_freq)
        self.form_freq = form_freq
        self.vocab = vocab
        self.form_vocab = form_vocab

    def __deepcopy__(self, memo):
        return self

    def __copy__(self):
        return self


# --------------------------------------------------------------------------- #
# paragraph structure and layout nulls
# --------------------------------------------------------------------------- #


def paragraph_lengths(lines: pd.DataFrame) -> dict[str, list[int]]:
    """Lines per paragraph of every page, in order of position."""
    out = {}
    for folio, g in lines.groupby("folio", sort=False):
        g = g.sort_values("pos")
        out[folio] = g.groupby("par", sort=True).size().tolist()
    return out


def structure_from_lengths(lengths: Sequence[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``par``, ``lip`` and ``par_len`` of every position of a page with paragraphs of the given lengths.

    >>> [a.tolist() for a in structure_from_lengths([2, 1])]
    [[0, 0, 1], [0, 1, 0], [2, 2, 1]]
    """
    par = np.repeat(np.arange(len(lengths)), lengths)
    lip = np.concatenate([np.arange(k) for k in lengths]) if len(lengths) else np.empty(0, dtype=int)
    par_len = np.repeat(np.asarray(lengths), lengths)
    return par.astype(int), lip.astype(int), par_len.astype(int)


def _apply_structure(lines: pd.DataFrame, lengths: Mapping[str, Sequence[int]]) -> pd.DataFrame:
    out = lines.copy()
    par = np.zeros(len(out), dtype=int)
    lip = np.zeros(len(out), dtype=int)
    plen = np.zeros(len(out), dtype=int)
    for folio, idx in out.groupby("folio", sort=False).indices.items():
        p, l, pl = structure_from_lengths(lengths[folio])
        pos = out["pos"].to_numpy()[idx]
        par[idx], lip[idx], plen[idx] = p[pos], l[pos], pl[pos]
    out["par"], out["lip"], out["par_len"] = par, lip, plen
    return out


def permute_lines(lines: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Null (b): line contents permuted among the positions of their page; the paragraph structure stays by position."""
    rng = np.random.default_rng(seed)
    out = lines.copy()
    pos = out["pos"].to_numpy().copy()
    for _, idx in out.groupby("folio", sort=False).indices.items():
        pos[idx] = pos[idx][rng.permutation(len(idx))]
    out["pos"] = pos
    return _apply_structure(out, paragraph_lengths(lines))


def shuffle_boundaries(lines: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Null (a): the order of the paragraph lengths of every page is permuted (same number and lengths of paragraphs)."""
    rng = np.random.default_rng(seed)
    lengths = paragraph_lengths(lines)
    new = {f: [lens[i] for i in rng.permutation(len(lens))] for f, lens in lengths.items()}
    return _apply_structure(lines, new)


# --------------------------------------------------------------------------- #
# pair statistics
# --------------------------------------------------------------------------- #


def _expect_matrix(ns: np.ndarray, p: np.ndarray) -> np.ndarray:
    """E[shared distinct units] for all pairs of lines with ``ns`` tokens, draws from frequencies ``p``."""
    a = 1.0 - np.power(1.0 - p[None, :], ns[:, None])
    return a @ a.T


def _form_expect(form_freq: Mapping[int, Mapping[int, int]]) -> dict[int, float]:
    """P(two independent occurrences of the unit have the same form) from the corpus-wide form distribution."""
    out = {}
    for u, ff in form_freq.items():
        tot = float(sum(ff.values()))
        out[u] = sum((c / tot) ** 2 for c in ff.values()) if tot else 0.0
    return out


def pair_content(lines: pd.DataFrame, expect: str = "both") -> pd.DataFrame:
    """Observed and expected sharing for every pair of lines of the same page (content only, no layout).

    Columns: ``folio, a, b`` (row indices into ``lines``), ``obs`` (shared distinct units), ``e_page``,
    ``e_glob`` (expectations from the page's and the corpus' unit frequencies), ``n_occ`` (pairs of
    occurrences of shared units), ``n_ident`` (of them with identical form), ``e_ident`` (expected
    identical from the units' corpus paradigms).  ``expect``: ``'page'``, ``'glob'`` or ``'both'``.
    """
    st = lines.attrs["stats"]
    freq = np.asarray(st.unit_freq, dtype=float)
    p_glob = freq / freq.sum() if freq.sum() else freq
    f_ident = _form_expect(st.form_freq)
    units = lines["units"].tolist()
    counts = lines["counts"].tolist()
    ns = lines["n"].to_numpy()
    rows = []
    for folio, idx in lines.groupby("folio", sort=False).indices.items():
        idx = np.asarray(idx)
        if len(idx) < 2:
            continue
        n_i = ns[idx]
        if expect in ("page", "both"):
            page_cnt = np.zeros(len(freq))
            for li in idx:
                for u, ff in counts[li].items():
                    page_cnt[u] += sum(ff.values())
            tot = page_cnt.sum()
            present = np.flatnonzero(page_cnt)
            e_page = _expect_matrix(n_i, page_cnt[present] / tot) if tot else np.zeros((len(idx), len(idx)))
        if expect in ("glob", "both"):
            e_glob = _expect_matrix(n_i, p_glob)
        for x in range(len(idx)):
            ua, ca = units[idx[x]], counts[idx[x]]
            for y in range(x + 1, len(idx)):
                shared = np.intersect1d(ua, units[idx[y]], assume_unique=True)
                n_occ = n_id = 0
                e_id = 0.0
                if len(shared):
                    cb = counts[idx[y]]
                    for u in shared.tolist():
                        fa, fb = ca[u], cb[u]
                        ta, tb = sum(fa.values()), sum(fb.values())
                        n_occ += ta * tb
                        n_id += sum(c * fb.get(f, 0) for f, c in fa.items())
                        e_id += ta * tb * f_ident[u]
                rows.append((folio, int(idx[x]), int(idx[y]), int(len(shared)),
                             float(e_page[x, y]) if expect in ("page", "both") else np.nan,
                             float(e_glob[x, y]) if expect in ("glob", "both") else np.nan,
                             n_occ, n_id, e_id))
    out = pd.DataFrame(rows, columns=["folio", "a", "b", "obs", "e_page", "e_glob", "n_occ", "n_ident", "e_ident"])
    return out


def pair_layout(pairs: pd.DataFrame, lines: pd.DataFrame) -> pd.DataFrame:
    """Add the layout of ``lines`` to a pair table: ``i, j`` (earlier/later row), ``d``, ``same`` and the
    paragraph attributes of both lines (``par_i, par_j, lip_i, lip_j, parlen_i, n_pars``)."""
    pos = lines["pos"].to_numpy()
    par = lines["par"].to_numpy()
    lip = lines["lip"].to_numpy()
    plen = lines["par_len"].to_numpy()
    npars = lines["n_pars"].to_numpy()
    a, b = pairs["a"].to_numpy(), pairs["b"].to_numpy()
    swap = pos[a] > pos[b]
    i = np.where(swap, b, a)
    j = np.where(swap, a, b)
    out = pairs.copy()
    out["i"], out["j"] = i, j
    out["d"] = pos[j] - pos[i]
    out["same"] = par[i] == par[j]
    out["par_i"], out["par_j"] = par[i], par[j]
    out["lip_i"], out["lip_j"] = lip[i], lip[j]
    out["parlen_i"], out["parlen_j"] = plen[i], plen[j]
    out["n_pars"] = npars[i]
    return out


def _w(pairs: pd.DataFrame, weights: Mapping[str, float] | None) -> np.ndarray:
    if weights is None:
        return np.ones(len(pairs))
    return pairs["folio"].map(weights).fillna(0.0).to_numpy(dtype=float)


def aggregate(pairs: pd.DataFrame, by: Sequence[str], weights: Mapping[str, float] | None = None,
              expect: str = "page") -> pd.DataFrame:
    """Sharing ratio per group: ``n`` (weighted pairs), ``obs``, ``exp``, ``s = obs / exp``."""
    w = _w(pairs, weights)
    df = pd.DataFrame({k: pairs[k].to_numpy() for k in by})
    df["n"] = w
    df["obs"] = pairs["obs"].to_numpy() * w
    df["exp"] = pairs[f"e_{expect}"].to_numpy() * w
    g = df.groupby(list(by), sort=True, as_index=False).sum()
    g["s"] = np.where(g["exp"] > 0, g["obs"] / g["exp"].where(g["exp"] > 0, 1.0), np.nan)
    return g


def sharing_by_distance(pairs: pd.DataFrame, max_d: int = 8, min_pars: int = 2,
                        weights: Mapping[str, float] | None = None, expect: str = "page",
                        exclude_first: bool = False) -> pd.DataFrame:
    """Metric 1: sharing ratio by distance ``d = 1..max_d`` and ``same`` paragraph, pages with ``>= min_pars``
    paragraphs.  ``exclude_first`` drops every pair with a paragraph-first line (artefact check: a generator
    that treats first lines specially gives a step at the boundary through those lines alone)."""
    sel = pairs[(pairs["d"] <= max_d) & (pairs["n_pars"] >= min_pars)]
    if exclude_first:
        sel = sel[(sel["lip_i"] >= 1) & (sel["lip_j"] >= 1)]
    g = aggregate(sel, ["d", "same"], weights, expect)
    full = pd.MultiIndex.from_product([range(1, max_d + 1), [False, True]], names=["d", "same"]).to_frame(index=False)
    return full.merge(g, on=["d", "same"], how="left")


def boundary_step(table: pd.DataFrame, max_d: int = 4) -> dict[int, float]:
    """``Delta(d) = s_same(d) - s_diff(d)`` from a :func:`sharing_by_distance` table."""
    piv = table.pivot(index="d", columns="same", values="s")
    out = {}
    for d in range(1, max_d + 1):
        if d in piv.index and True in piv.columns and False in piv.columns:
            out[d] = float(piv.loc[d, True] - piv.loc[d, False])
        else:
            out[d] = float("nan")
    return out


def _slope(x: np.ndarray, y: np.ndarray, w: np.ndarray | None = None) -> float:
    ok = np.isfinite(y)
    x, y = x[ok], y[ok]
    w = np.ones(len(x)) if w is None else np.asarray(w)[ok]
    if len(x) < 2 or w.sum() <= 0:
        return float("nan")
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    den = np.sum(w * (x - xm) ** 2)
    return float(np.sum(w * (x - xm) * (y - ym)) / den) if den > 0 else float("nan")


def within_par_decay(pairs: pd.DataFrame, max_d: int = 6, min_pars: int = 1,
                     weights: Mapping[str, float] | None = None, expect: str = "page") -> tuple[pd.DataFrame, float]:
    """Metric 2: ``s(d)`` for pairs of the same paragraph, ``d = 1..max_d``, and the slope of ``s`` over ``d`` (weighted by pairs)."""
    sel = pairs[(pairs["d"] <= max_d) & pairs["same"] & (pairs["n_pars"] >= min_pars)]
    g = aggregate(sel, ["d"], weights, expect)
    full = pd.DataFrame({"d": range(1, max_d + 1)}).merge(g, on="d", how="left")
    slope = _slope(full["d"].to_numpy(dtype=float), full["s"].to_numpy(dtype=float), full["n"].fillna(0).to_numpy())
    return full, slope


def first_line_sharing(pairs: pd.DataFrame, max_d: int = 6, weights: Mapping[str, float] | None = None,
                       expect: str = "page") -> pd.DataFrame:
    """Metric 3: sharing of the paragraph's first line with later lines of its own paragraph (``first_own``)
    and with lines of the *next* paragraph (``first_next``), and the same for a middle line (``lip >= 1``,
    not the last of its paragraph: ``mid_own``, ``mid_next``), by distance ``d = 1..max_d``.

    Returns rows ``kind, d, n, obs, exp, s``.  "Announcing" = ``first_own - first_next`` larger than
    ``mid_own - mid_next`` at the same distances; a mere edge effect lowers ``first_own`` and ``first_next`` alike.
    """
    sel = pairs[pairs["d"] <= max_d]
    first_i = sel["lip_i"] == 0
    mid_i = (sel["lip_i"] >= 1) & (sel["lip_i"] < sel["parlen_i"] - 1)
    nxt = sel["par_j"] == sel["par_i"] + 1
    kinds = {
        "first_own": sel[first_i & sel["same"]],
        "first_next": sel[first_i & nxt],
        "mid_own": sel[mid_i & sel["same"]],
        "mid_next": sel[mid_i & nxt],
    }
    parts = []
    for kind, sub in kinds.items():
        g = aggregate(sub, ["d"], weights, expect) if len(sub) else pd.DataFrame(columns=["d", "n", "obs", "exp", "s"])
        full = pd.DataFrame({"d": range(1, max_d + 1)}).merge(g, on="d", how="left")
        full.insert(0, "kind", kind)
        parts.append(full)
    return pd.concat(parts, ignore_index=True)


def page_boundary(lines: pd.DataFrame, pairs: pd.DataFrame, max_d: int = 3,
                  weights: Mapping[str, float] | None = None, open_end: Mapping[str, bool] | None = None) -> pd.DataFrame:
    """Metric 4: sharing between the last lines of page ``i`` and the first lines of the next page in the
    manuscript (``rank + 1``) at line distance ``d = 1..max_d`` (last/first = 1), with the *global*
    expectation, against within-page pairs at the same ``d`` (also global expectation).

    Rows: ``group`` in ``across_all``, ``across_same`` (next page of the same language and section),
    ``across_diff``, ``across_open`` (pages whose last paragraph has no end mark, ``open_end[folio]``
    True; only when given), ``within``, ``within_cross`` (within page, across a paragraph boundary);
    columns ``d, n, obs, exp, s``.  Weights apply to the earlier page.
    """
    freq = np.asarray(lines.attrs["stats"].unit_freq, dtype=float)
    p_glob = freq / freq.sum() if freq.sum() else freq
    by_rank = {r: idx for r, idx in lines.groupby("rank", sort=False).indices.items()}
    pos = lines["pos"].to_numpy()
    ns = lines["n"].to_numpy()
    units = lines["units"].tolist()
    meta = lines[["folio", "language", "illustration"]].astype(str).to_numpy()
    rows = []
    for r, idx in by_rank.items():
        if r + 1 not in by_rank:
            continue
        nxt = by_rank[r + 1]
        idx = idx[np.argsort(pos[idx])]
        nxt = nxt[np.argsort(pos[nxt])]
        last = idx[::-1][:max_d]
        first = nxt[:max_d]
        both = np.concatenate([last, first])
        e = _expect_matrix(ns[both], p_glob)
        same_attr = (meta[idx[0], 1] == meta[nxt[0], 1]) and (meta[idx[0], 2] == meta[nxt[0], 2])
        folio = meta[idx[0], 0]
        w = 1.0 if weights is None else float(weights.get(folio, 0.0))
        is_open = bool(open_end.get(folio, False)) if open_end is not None else False
        for a in range(len(last)):
            for b in range(len(first)):
                d = a + b + 1
                if d > max_d:
                    continue
                obs = len(np.intersect1d(units[last[a]], units[first[b]], assume_unique=True))
                rows.append(("across_all", d, w, obs * w, e[a, len(last) + b] * w))
                rows.append(("across_same" if same_attr else "across_diff", d, w, obs * w, e[a, len(last) + b] * w))
                if is_open:
                    rows.append(("across_open", d, w, obs * w, e[a, len(last) + b] * w))
    across = pd.DataFrame(rows, columns=["group", "d", "n", "obs", "exp"])
    within = pairs[pairs["d"] <= max_d]
    g1 = aggregate(within, ["d"], weights, "glob").assign(group="within")
    g2 = aggregate(within[~within["same"]], ["d"], weights, "glob").assign(group="within_cross")
    agg = across.groupby(["group", "d"], as_index=False).sum() if len(across) else across
    out = pd.concat([agg, g1[["group", "d", "n", "obs", "exp"]], g2[["group", "d", "n", "obs", "exp"]]], ignore_index=True)
    groups = ["across_all", "across_same", "across_diff", "within", "within_cross"] + (["across_open"] if open_end is not None else [])
    full = pd.MultiIndex.from_product([groups, range(1, max_d + 1)], names=["group", "d"]).to_frame(index=False)
    out = full.merge(out, on=["group", "d"], how="left")
    out["s"] = np.where(out["exp"] > 0, out["obs"] / out["exp"].where(out["exp"] > 0, 1.0), np.nan)
    return out


def form_by_boundary(pairs: pd.DataFrame, max_d: int = 3, min_pars: int = 2,
                     weights: Mapping[str, float] | None = None) -> pd.DataFrame:
    """Metric 5: among occurrence pairs of shared units, the share with identical form minus the paradigm
    expectation, for pairs in the same paragraph and across a boundary, ``d = 1..max_d`` and pooled (``d = 0``)."""
    sel = pairs[(pairs["d"] <= max_d) & (pairs["n_pars"] >= min_pars) & (pairs["n_occ"] > 0)]
    w = _w(sel, weights)
    df = pd.DataFrame({"d": sel["d"].to_numpy(), "same": sel["same"].to_numpy(), "n_occ": sel["n_occ"].to_numpy() * w,
                       "n_ident": sel["n_ident"].to_numpy() * w, "e_ident": sel["e_ident"].to_numpy() * w})
    g = df.groupby(["d", "same"], as_index=False).sum()
    pooled = df.groupby("same", as_index=False).sum().assign(d=0)
    g = pd.concat([g, pooled[g.columns]], ignore_index=True)
    full = pd.MultiIndex.from_product([range(0, max_d + 1), [False, True]], names=["d", "same"]).to_frame(index=False)
    out = full.merge(g, on=["d", "same"], how="left")
    out["ident"] = out["n_ident"] / out["n_occ"]
    out["excess"] = (out["n_ident"] - out["e_ident"]) / out["n_occ"]
    return out


def all_metrics(lines: pd.DataFrame, pairs: pd.DataFrame, weights: Mapping[str, float] | None = None,
                open_end: Mapping[str, bool] | None = None, prefix: str = "") -> dict:
    """All metrics of one sample as a flat mapping (keys ``<prefix><metric>__<key>``), for one ``runs.csv`` row.

    ``pairs`` is :func:`pair_layout` of :func:`pair_content` of ``lines``.  Expectation ``page`` (``p``)
    and ``glob`` (``g``) where both apply.
    """
    out: dict = {}

    def put(key: str, val) -> None:
        out[f"{prefix}{key}"] = round(float(val), 5) if val is not None and np.isfinite(val) else np.nan

    for tag, expect in (("p", "page"), ("g", "glob")):
        t = sharing_by_distance(pairs, 8, 2, weights, expect)
        for _, r in t.iterrows():
            put(f"share_{tag}__d{int(r['d'])}_{'same' if r['same'] else 'diff'}", r["s"])
            if tag == "p":
                put(f"share_n__d{int(r['d'])}_{'same' if r['same'] else 'diff'}", r["n"])
        for d, v in boundary_step(t, 4).items():
            put(f"step_{tag}__d{d}", v)
        t_nf = sharing_by_distance(pairs, 4, 2, weights, expect, exclude_first=True)
        for _, r in t_nf.iterrows():
            put(f"sharenf_{tag}__d{int(r['d'])}_{'same' if r['same'] else 'diff'}", r["s"])
            if tag == "p":
                put(f"sharenf_n__d{int(r['d'])}_{'same' if r['same'] else 'diff'}", r["n"])
        for d, v in boundary_step(t_nf, 4).items():
            put(f"stepnf_{tag}__d{d}", v)
        dec, slope = within_par_decay(pairs, 6, 1, weights, expect)
        for _, r in dec.iterrows():
            put(f"decay_{tag}__d{int(r['d'])}", r["s"])
            if tag == "p":
                put(f"decay_n__d{int(r['d'])}", r["n"])
        put(f"decay_{tag}__slope", slope)
        fl = first_line_sharing(pairs, 6, weights, expect)
        for _, r in fl.iterrows():
            put(f"first_{tag}__{r['kind']}_d{int(r['d'])}", r["s"])
            if tag == "p":
                put(f"first_n__{r['kind']}_d{int(r['d'])}", r["n"])
    pb = page_boundary(lines, pairs, 3, weights, open_end)
    for _, r in pb.iterrows():
        put(f"page__{r['group']}_d{int(r['d'])}", r["s"])
        put(f"page_n__{r['group']}_d{int(r['d'])}", r["n"])
    piv = pb.pivot(index="d", columns="group", values="s")
    for d in range(1, 4):
        for grp in ("across_all", "across_same", "across_diff"):
            if grp in piv.columns and "within" in piv.columns and piv.loc[d, "within"] > 0:
                put(f"page_ratio__{grp}_d{d}", piv.loc[d, grp] / piv.loc[d, "within"])
            else:
                put(f"page_ratio__{grp}_d{d}", float("nan"))
    fb = form_by_boundary(pairs, 3, 2, weights)
    for _, r in fb.iterrows():
        tag = "same" if r["same"] else "diff"
        put(f"form__excess_{tag}_d{int(r['d'])}", r["excess"])
        put(f"form__ident_{tag}_d{int(r['d'])}", r["ident"])
        put(f"form_n__{tag}_d{int(r['d'])}", r["n_occ"])
    return out


def pairs_of(lines: pd.DataFrame) -> pd.DataFrame:
    """:func:`pair_layout` of :func:`pair_content` – the pair table of a sample."""
    return pair_layout(pair_content(lines), lines)


__all__ += ["pairs_of", "LineStats"]
