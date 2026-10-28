"""Recurrence of rare stems at a distance (experiment 015).

The unit of analysis is a *stem occurrence*: a paragraph token whose stem (``decompose_table``)
has between ``min_count`` and ``max_count`` occurrences in the corpus.  Every pair of occurrences
of the same stem has a page distance (``|rank_i - rank_j|`` in file order; 0 = the same page) and
covariates of the two pages (section, language, hand, quire, bifolio) and of the two tokens
(identical form or not).  The functions here build the occurrence table and its pairs, bin the
distances, and give the expectations under the token-to-page permutation null (every token of a
stem is assigned to a random page with probability proportional to page size, the stem's
frequency being kept); most expectations are exact (no sampling), the sampled permutation
(:func:`permute_pages`) is used for intervals and as a check.

Distance bins: 0 / 1 / 2-3 / 4-7 / 8-15 / 16-31 / 32+ (:data:`DISTANCE_EDGES`); "far" = >= 16
pages (:data:`FAR`); form classes 0 / 1-3 / 4-15 / 16+ (:data:`FORM_EDGES`).
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd

DISTANCE_EDGES = (0, 1, 2, 4, 8, 16, 32)
DISTANCE_LABELS = ("0", "1", "2-3", "4-7", "8-15", "16-31", "32+")
FORM_EDGES = (0, 1, 4, 16)
FORM_LABELS = ("0", "1-3", "4-15", "16+")
FAR = 16
NONLOCAL = 4
PAGE_COLUMNS = ("illustration", "language", "hand", "quire", "bifolio")
SAME_COLUMNS = ("section", "language", "hand", "quire", "bifolio")
LABEL_SUBTYPES = ("z", "f", "c", "s", "n")


# --------------------------------------------------------------------------- #
# pages
# --------------------------------------------------------------------------- #


def page_order(tokens: pd.DataFrame) -> dict:
    """``folio -> rank`` over all pages of a token table in file order (all locus types)."""
    return {p: i for i, p in enumerate(pd.unique(tokens["folio"].to_numpy()))}


def page_table(df: pd.DataFrame, rank: dict | None = None) -> pd.DataFrame:
    """One row per page of ``df`` (folio index): ``rank``, ``size`` (tokens), ``weight`` and metadata.

    ``rank`` is the file order of the pages (:func:`page_order`; default: order of appearance in
    ``df``).  ``section`` is the ``illustration`` code; ``bifolio`` is made unique as ``quire/bifolio``.
    Missing metadata columns are filled with ``''``.
    """
    pages = pd.unique(df["folio"].to_numpy())
    rank = rank or {p: i for i, p in enumerate(pages)}
    g = df.groupby("folio", sort=False)
    out = pd.DataFrame({"size": g.size().astype(float)}).loc[pages]
    out["rank"] = [int(rank[p]) for p in pages]
    for c in PAGE_COLUMNS:
        out[c] = g[c].first().astype(str).loc[pages] if c in df.columns else ""
    out["section"] = out.pop("illustration")
    out["bifolio"] = out["quire"] + "/" + out["bifolio"]
    out["weight"] = out["size"] / out["size"].sum()
    return out


def _edges_to_bins(x: np.ndarray, edges: Sequence[int], labels: Sequence[str]) -> np.ndarray:
    idx = np.searchsorted(np.asarray(edges), x, side="right") - 1
    return np.asarray(labels, dtype=object)[idx]


def distance_bin(dist: Iterable[int]) -> np.ndarray:
    """Labels of :data:`DISTANCE_LABELS` for page distances."""
    return _edges_to_bins(np.asarray(list(dist), dtype=int), DISTANCE_EDGES, DISTANCE_LABELS)


def form_bin(dist: Iterable[int]) -> np.ndarray:
    """Labels of :data:`FORM_LABELS` (coarser classes for the form-by-distance metric)."""
    return _edges_to_bins(np.asarray(list(dist), dtype=int), FORM_EDGES, FORM_LABELS)


# --------------------------------------------------------------------------- #
# occurrences and pairs
# --------------------------------------------------------------------------- #


def rare_stems(dfd: pd.DataFrame, min_count: int = 2, max_count: int = 5,
               within: Iterable[str] | None = None) -> pd.Series:
    """Counts of the non-empty stems with ``min_count <= count <= max_count`` (optionally only ``within``)."""
    vc = dfd.loc[dfd["stem"] != "", "stem"].value_counts()
    vc = vc[(vc >= min_count) & (vc <= max_count)]
    if within is not None:
        vc = vc[vc.index.isin(set(within))]
    return vc


def stem_occurrences(dfd: pd.DataFrame, pages: pd.DataFrame, min_count: int = 2, max_count: int = 5,
                     within: Iterable[str] | None = None) -> pd.DataFrame:
    """Occurrence table of the rare stems: one row per token, sorted by stem.

    Columns: ``stem``, ``token``, ``folio``, ``count`` (of the stem) and, from ``pages``
    (:func:`page_table`), ``rank``, ``section``, ``language``, ``hand``, ``quire``, ``bifolio``.
    """
    counts = rare_stems(dfd, min_count, max_count, within)
    keep = dfd["stem"].isin(counts.index)
    occ = dfd.loc[keep, ["stem", "token", "folio"]].reset_index(drop=True)
    occ["count"] = occ["stem"].map(counts).astype(int)
    occ = occ.sort_values(["stem", "folio"], kind="stable").reset_index(drop=True)
    return relocate(occ, pages, occ["folio"].to_numpy())


def relocate(occ: pd.DataFrame, pages: pd.DataFrame, folio: np.ndarray) -> pd.DataFrame:
    """Copy of ``occ`` with every occurrence moved to the page ``folio[i]`` (rank and metadata of ``pages``)."""
    out = occ.copy()
    out["folio"] = folio
    p = pages.loc[folio]
    out["rank"] = p["rank"].to_numpy()
    for c in SAME_COLUMNS:
        out[c] = p[c].to_numpy()
    return out


def permute_pages(occ: pd.DataFrame, pages: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """Permutation null: every occurrence gets a random page, probability proportional to page size.

    Stems keep their tokens (forms) and frequencies; pages are drawn with replacement, so two
    occurrences of a stem can land on the same page (distance 0) as they can in the real text.
    """
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(pages), size=len(occ), p=pages["weight"].to_numpy())
    return relocate(occ, pages, pages.index.to_numpy()[idx])


def occurrence_pairs(occ: pd.DataFrame) -> pd.DataFrame:
    """All pairs of occurrences of the same stem.

    Columns: ``stem``, ``i``, ``j`` (row numbers in ``occ``), ``dist`` (pages), ``same_form``
    (identical token), ``same_section`` / ``same_language`` / ``same_hand`` / ``same_quire`` /
    ``same_bifolio``.  ``occ`` must be sorted by stem (:func:`stem_occurrences` is).
    """
    stems = occ["stem"].to_numpy()
    n = len(occ)
    if n == 0:
        return pd.DataFrame(columns=["stem", "i", "j", "dist", "same_form"] + [f"same_{c}" for c in SAME_COLUMNS])
    change = np.flatnonzero(np.r_[True, stems[1:] != stems[:-1]])
    starts, sizes = change, np.diff(np.r_[change, n])
    ii, jj = [], []
    for k in np.unique(sizes):
        if k < 2:
            continue
        a, b = np.triu_indices(int(k), 1)
        s = starts[sizes == k][:, None]
        ii.append((s + a[None, :]).ravel())
        jj.append((s + b[None, :]).ravel())
    if not ii:
        return pd.DataFrame(columns=["stem", "i", "j", "dist", "same_form"] + [f"same_{c}" for c in SAME_COLUMNS])
    i, j = np.concatenate(ii), np.concatenate(jj)
    rank = occ["rank"].to_numpy()
    out = pd.DataFrame({"stem": stems[i], "i": i, "j": j, "dist": np.abs(rank[i] - rank[j])})
    tok = occ["token"].to_numpy()
    out["same_form"] = tok[i] == tok[j]
    for c in SAME_COLUMNS:
        v = occ[c].to_numpy()
        out[f"same_{c}"] = v[i] == v[j]
    return out


# --------------------------------------------------------------------------- #
# exact expectations under the page permutation
# --------------------------------------------------------------------------- #


def _pair_weights(pages: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Distance and weight ``w_p w_q`` of every ordered page pair (p, q), p = q included."""
    r, w = pages["rank"].to_numpy(), pages["weight"].to_numpy()
    d = np.abs(r[:, None] - r[None, :]).ravel()
    return d, (w[:, None] * w[None, :]).ravel()


def expected_profile(pages: pd.DataFrame) -> dict:
    """Expected share of pairs per distance bin when both occurrences are random size-weighted pages.

    Also ``far`` (share >= :data:`FAR`) and ``far_nonlocal`` (share >= FAR among pairs >= :data:`NONLOCAL`).
    """
    d, w = _pair_weights(pages)
    bins = distance_bin(d)
    out = {b: float(w[bins == b].sum()) for b in DISTANCE_LABELS}
    out["far"] = float(w[d >= FAR].sum())
    nl = float(w[d >= NONLOCAL].sum())
    out["far_nonlocal"] = out["far"] / nl if nl > 0 else np.nan
    return out


def observed_profile(pairs: pd.DataFrame) -> dict:
    """Share of pairs per distance bin, ``far`` and ``far_nonlocal`` (see :func:`expected_profile`)."""
    d = pairs["dist"].to_numpy()
    n = len(d)
    if n == 0:
        return {**{b: np.nan for b in DISTANCE_LABELS}, "far": np.nan, "far_nonlocal": np.nan}
    bins = distance_bin(d)
    out = {b: float((bins == b).mean()) for b in DISTANCE_LABELS}
    out["far"] = float((d >= FAR).mean())
    nl = int((d >= NONLOCAL).sum())
    out["far_nonlocal"] = float((d >= FAR).sum() / nl) if nl else np.nan
    return out


def same_by_distance(pages: pd.DataFrame) -> pd.DataFrame:
    """Exact ``P(same section/language/hand/quire/bifolio | distance d)`` for random size-weighted page pairs.

    Indexed by distance; the distance-matched expectation of :func:`far_pairs`.
    """
    d, w = _pair_weights(pages)
    cols = {}
    for c in SAME_COLUMNS:
        v = pages[c].to_numpy()
        cols[f"same_{c}"] = (v[:, None] == v[None, :]).ravel()
    df = pd.DataFrame({"dist": d, "w": w, **cols})
    g = df.groupby("dist")
    out = pd.DataFrame({k: g.apply(lambda x, k=k: float((x[k] * x["w"]).sum() / x["w"].sum())) for k in cols})
    return out


def far_pairs(pairs: pd.DataFrame, expect: pd.DataFrame, lo: int = FAR, hi: int | None = None) -> dict:
    """Share of pairs at distance ``lo..hi`` on the same section/language/hand/quire/bifolio and its expectation.

    The expectation is distance matched: the mean over the selected pairs of
    ``expect.loc[dist]`` (:func:`same_by_distance`), i.e. the share random page pairs at exactly the
    same distances would give.  Returns ``n`` and, per column, ``<col>`` and ``<col>_exp``.
    """
    m = pairs["dist"].to_numpy() >= lo
    if hi is not None:
        m &= pairs["dist"].to_numpy() <= hi
    sel = pairs[m]
    out = {"n": int(m.sum())}
    e = expect.reindex(sel["dist"].to_numpy())
    for c in SAME_COLUMNS:
        k = f"same_{c}"
        out[c] = float(sel[k].mean()) if len(sel) else np.nan
        out[f"{c}_exp"] = float(e[k].mean()) if len(sel) else np.nan
    return out


def form_by_distance(occ: pd.DataFrame, pairs: pd.DataFrame) -> dict:
    """Share of identical-form pairs per form class and its expectation from the stem's paradigm.

    For a stem with ``k`` occurrences and form counts ``n_f`` the probability that two occurrences
    drawn without replacement are identical is ``sum n_f (n_f - 1) / (k (k - 1))``; the expectation
    of a class is the mean of this over its pairs (form chosen independently of the distance).
    Returns ``n_<class>``, ``same_<class>`` and ``exp_<class>``.
    """
    counts = occ.groupby(["stem", "token"], sort=False).size()
    k = occ.groupby("stem", sort=False).size()
    num = (counts * (counts - 1)).groupby(level="stem").sum()
    p_same = (num / (k * (k - 1))).reindex(pairs["stem"].to_numpy()).to_numpy()
    cls = form_bin(pairs["dist"].to_numpy())
    same = pairs["same_form"].to_numpy()
    out = {}
    for c in FORM_LABELS:
        m = cls == c
        out[f"n_{c}"] = int(m.sum())
        out[f"same_{c}"] = float(same[m].mean()) if m.any() else np.nan
        out[f"exp_{c}"] = float(p_same[m].mean()) if m.any() else np.nan
    return out


def _prob_at_most_two(weights: np.ndarray, k: int) -> float:
    """P(k independent draws from ``weights`` hit at most two distinct categories)."""
    w = np.asarray(weights, dtype=float)
    w = w[w > 0]
    one = float((w ** k).sum())
    s = (w[:, None] + w[None, :]) ** k - w[:, None] ** k - w[None, :] ** k
    two = float(np.triu(s, 1).sum())
    return one + two


def far_clusters(occ: pd.DataFrame, pages: pd.DataFrame, by: str = "quire", min_groups: int = 3) -> dict:
    """Stems whose occurrences spread over at least ``min_groups`` distinct ``by`` groups (quires, sections).

    ``n`` observed, ``share`` of the stems with >= ``min_groups`` occurrences, ``exp`` (exact
    expectation under the size-weighted page permutation, ``min_groups = 3`` only) and ``n_stems``.
    """
    k = occ.groupby("stem", sort=False).size()
    k = k[k >= min_groups]
    groups = occ.groupby("stem", sort=False)[by].nunique().reindex(k.index)
    n = int((groups >= min_groups).sum())
    out = {"n_stems": int(len(k)), "n": n, "share": float(n / len(k)) if len(k) else np.nan}
    if min_groups == 3 and len(k):
        w = pages.groupby(by)["weight"].sum().to_numpy()
        exp = sum((1.0 - _prob_at_most_two(w, int(kk))) * cnt for kk, cnt in k.value_counts().items())
        out["exp"] = float(exp)
        out["share_exp"] = float(exp / len(k))
    else:
        out["exp"], out["share_exp"] = np.nan, np.nan
    return out


# --------------------------------------------------------------------------- #
# labels
# --------------------------------------------------------------------------- #


def label_stems(dfd_labels: pd.DataFrame, subtypes: Sequence[str] = LABEL_SUBTYPES) -> set:
    """Non-empty stems of the label tokens (``generic_type == 'L'``) with ``subtype`` in ``subtypes``."""
    m = dfd_labels["generic_type"].astype(str) == "L"
    if "subtype" in dfd_labels.columns and subtypes:
        m &= dfd_labels["subtype"].astype(str).isin(list(subtypes))
    s = dfd_labels.loc[m, "stem"]
    return set(s[s != ""].tolist())


# --------------------------------------------------------------------------- #
# all metrics of one occurrence table
# --------------------------------------------------------------------------- #


def recurrence_metrics(occ: pd.DataFrame, pages: pd.DataFrame, expect_same: pd.DataFrame | None = None,
                       expect_profile: dict | None = None) -> dict:
    """Flat dictionary of the 015 metrics of one occurrence table.

    Keys: ``n_stems``, ``n_tokens``, ``n_pairs``; ``p_<bin>`` / ``e_<bin>`` observed and expected share
    per distance bin, ``p_far``/``e_far``/``x_far`` (excess), ``p_far_nl``/``e_far_nl``/``x_far_nl``;
    ``far_<col>``/``far_<col>_exp`` and ``mid_<col>``/``mid_<col>_exp`` (distance 2..15) for the
    same-page-attribute shares; ``form_*`` (:func:`form_by_distance`); ``clq_*`` / ``cls_*``
    (:func:`far_clusters` by quire / section).  ``expect_same`` and ``expect_profile`` can be
    precomputed (they depend on ``pages`` only).
    """
    pairs = occurrence_pairs(occ)
    exp_prof = expect_profile or expected_profile(pages)
    expect_same = same_by_distance(pages) if expect_same is None else expect_same
    out = {"n_stems": int(occ["stem"].nunique()), "n_tokens": int(len(occ)), "n_pairs": int(len(pairs))}
    obs = observed_profile(pairs)
    for b in DISTANCE_LABELS:
        out[f"p_{b}"], out[f"e_{b}"] = obs[b], exp_prof[b]
    out["p_far"], out["e_far"] = obs["far"], exp_prof["far"]
    out["x_far"] = obs["far"] - exp_prof["far"]
    out["p_far_nl"], out["e_far_nl"] = obs["far_nonlocal"], exp_prof["far_nonlocal"]
    out["x_far_nl"] = obs["far_nonlocal"] - exp_prof["far_nonlocal"]
    for tag, lo, hi in (("far", FAR, None), ("mid", 2, FAR - 1)):
        for k, v in far_pairs(pairs, expect_same, lo, hi).items():
            out[f"{tag}_{k}"] = v
    for k, v in form_by_distance(occ, pairs).items():
        out[f"form_{k}"] = v
    for tag, by in (("clq", "quire"), ("cls", "section")):
        for k, v in far_clusters(occ, pages, by).items():
            out[f"{tag}_{k}"] = v
    return out
