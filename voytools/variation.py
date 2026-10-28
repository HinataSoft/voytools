"""Within-word variation: how much of the choice of an affix given the stem is explained by covariates.

Experiment 011.  For a decomposed token table (columns ``prefix``, ``stem``, ``suffix`` from
:func:`voytools.morphology.decompose_table`) the functions here measure

* the decomposition of the conditional entropy H(S | Z) of a target S (suffix or prefix) given the
  stem Z into the part explained by covariates X (position in the line, Currier language, scribe,
  section, page ...) and the rest (:func:`conditional_entropy_decomposition`).  The plug-in
  reduction H(S|Z) - H(S|Z,X) is biased upwards for covariates with many values; the bias is
  estimated by permuting X among the tokens of the same stem (``n_perm`` times) and subtracted;
* the residual structure of the choice: dependence on a neighbour's affix given the stem and the
  conditioning covariates (:func:`neighbour_conditional_mi`), and the consistency of a stem's
  affix inside a page or paragraph (:func:`mode_consistency`).

Covariate columns in the classes of experiment 007 are added by :func:`add_covariates`.  All
estimators work on integer codes (:func:`codes`) and plug-in entropies in bits.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from ._common import line_ids, page_ids

# --------------------------------------------------------------------------- #
# codes and entropies
# --------------------------------------------------------------------------- #


def codes(values) -> np.ndarray:
    """Integer codes 0..k-1 of a sequence of hashable values (missing values form one code)."""
    s = pd.Series(list(values) if not isinstance(values, (pd.Series, np.ndarray)) else values, dtype=object)
    s = s.where(s.notna(), "<NA>")
    return pd.factorize(s)[0].astype(np.int64)


def joint_codes(*arrays: np.ndarray) -> np.ndarray:
    """Integer code of the joint value of several integer code arrays (re-factorised, stays small).

    >>> joint_codes(np.array([0, 0, 1, 1]), np.array([0, 1, 0, 0])).tolist()
    [0, 1, 2, 2]
    """
    if not arrays:
        raise ValueError("need at least one array")
    key = np.asarray(arrays[0], dtype=np.int64)
    for a in arrays[1:]:
        a = np.asarray(a, dtype=np.int64)
        key = pd.factorize(key * (int(a.max()) + 1 if len(a) else 1) + a)[0].astype(np.int64)
    return pd.factorize(key)[0].astype(np.int64)


def entropy(c: np.ndarray) -> float:
    """Plug-in entropy (bits) of an integer code array."""
    n = len(c)
    if n == 0:
        return float("nan")
    p = np.bincount(c) / n
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def _pair_entropy(a: np.ndarray, ka: int, b: np.ndarray) -> float:
    """H(A, B) for codes ``a`` (``ka`` values) and ``b``; avoids re-factorising."""
    key = a * (int(b.max()) + 1 if len(b) else 1) + b
    if ka * (int(b.max()) + 1) < 50_000_000:
        n = len(key)
        p = np.bincount(key) / n
        p = p[p > 0]
        return float(-(p * np.log2(p)).sum())
    return entropy(pd.factorize(key)[0])


def conditional_entropy(y: np.ndarray, *given: np.ndarray) -> float:
    """Plug-in H(Y | G1, G2, ...) in bits for integer code arrays (no ``given``: H(Y)).

    >>> conditional_entropy(np.array([0, 1, 0, 1]), np.array([0, 1, 0, 1]))
    0.0
    """
    y = np.asarray(y, dtype=np.int64)
    if not given:
        return entropy(y)
    g = joint_codes(*given)
    return _pair_entropy(g, int(g.max()) + 1, y) - entropy(g)


def permute_within(groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Index array that permutes rows among the rows of the same group (``x[idx]`` is the shuffled ``x``)."""
    groups = np.asarray(groups)
    n = len(groups)
    base = np.argsort(groups, kind="stable")
    order = np.lexsort((rng.random(n), groups))
    out = np.empty(n, dtype=np.int64)
    out[base] = order
    return out


# --------------------------------------------------------------------------- #
# covariates (classes of experiment 007)
# --------------------------------------------------------------------------- #


def line_place(df: pd.DataFrame) -> np.ndarray:
    """Place of the token in its line: ``'first'``, ``'mid'``, ``'last'`` or ``'single'``."""
    f = df["is_first"].to_numpy(dtype=bool)
    last = df["is_last"].to_numpy(dtype=bool)
    return np.where(f & last, "single", np.where(f, "first", np.where(last, "last", "mid"))).astype(object)


def position_class(df: pd.DataFrame) -> np.ndarray:
    """Position class of experiment 007: line place x ``rel_pos`` in 5 bins x ``line_in_par`` 0 / 1 / 2+.

    Identical to ``position_classes`` in ``experiments/007_position_vs_content/run.py`` (22 occupied
    classes in the ZL paragraph text).
    """
    place = line_place(df)
    rel = np.clip((df["rel_pos"].to_numpy(dtype=float) * 5).astype(np.int64), 0, 4)
    lip = pd.to_numeric(df["line_in_par"], errors="coerce").to_numpy(dtype=float)
    lip = np.where(np.isnan(lip), 2, np.minimum(lip, 2)).astype(np.int64)
    return np.asarray([f"{a}|{b}|{c}" for a, b, c in zip(place, rel, lip)], dtype=object)


def add_covariates(df: pd.DataFrame) -> pd.DataFrame:
    """Copy of a token table with the covariate columns used by experiment 011.

    ``cv_pos`` (position class, :func:`position_class`), ``cv_relpos`` (``rel_pos`` in 5 bins),
    ``cv_place`` (:func:`line_place`), ``cv_lip`` (``line_in_par`` 0 / 1 / 2+), ``cv_lang``
    (``language``), ``cv_hand`` (``hand``), ``cv_sec`` (``illustration``), ``cv_page`` (page id),
    ``cv_par`` (paragraph id: page + ``par_index``, a line without index is its own unit),
    ``cv_line`` (line id).
    """
    out = df.copy()
    out["cv_pos"] = position_class(df)
    out["cv_relpos"] = np.clip((df["rel_pos"].to_numpy(dtype=float) * 5).astype(np.int64), 0, 4)
    out["cv_place"] = line_place(df)
    lip = pd.to_numeric(df["line_in_par"], errors="coerce").to_numpy(dtype=float)
    out["cv_lip"] = np.where(np.isnan(lip), 2, np.minimum(lip, 2)).astype(np.int64)
    for new, col in (("cv_lang", "language"), ("cv_hand", "hand"), ("cv_sec", "illustration")):
        out[new] = df[col].astype(object).fillna("?").astype(str).to_numpy() if col in df.columns else "?"
    pg = page_ids(df)
    ln = line_ids(df)
    out["cv_page"] = pg
    out["cv_line"] = ln
    if "par_index" in df.columns:
        par = df["par_index"].tolist()
        keys = [f"p{a}|{b}" if not pd.isna(b) else f"l{c}" for a, b, c in zip(pg, par, ln)]
    else:
        keys = [f"p{a}" for a in pg]
    out["cv_par"] = codes(keys)
    return out


def _columns(spec) -> list[str]:
    return [spec] if isinstance(spec, str) else list(spec)


# --------------------------------------------------------------------------- #
# entropy decomposition
# --------------------------------------------------------------------------- #


def conditional_entropy_decomposition(
    dfd: pd.DataFrame,
    target: str = "suffix",
    given: str | Sequence[str] = "stem",
    by: Mapping[str, str | Sequence[str]] | Sequence[str] | None = None,
    n_perm: int = 100,
    min_count: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Share of H(target | given) explained by each covariate set in ``by``, bias-corrected by permutation.

    Parameters
    ----------
    dfd:
        Token table with the columns ``target``, ``given`` and the covariates.
    given:
        Column (or list of columns) conditioned on in both terms.  Rows are restricted to values of
        the *first* given column (the stem) with at least ``min_count`` tokens.
    by:
        Mapping name -> column or list of columns (a joint covariate X), or a list of column names
        (each its own X).
    n_perm:
        Permutations of X among the rows with the same ``given`` value; their mean
        H(target | given, X_perm) estimates the conditional entropy expected without association
        (the plug-in reduction of an unrelated X with many values is not zero).

    Returns one row per X with ``H_given`` = H(S|Z), ``H_cond`` = H(S|Z,X), ``H_cond_perm`` (mean
    over permutations), ``reduction_raw`` = H(S|Z) - H(S|Z,X), ``bias`` = H(S|Z) - H_cond_perm,
    ``reduction`` = H_cond_perm - H(S|Z,X) (corrected), ``explained`` = reduction / H(S|Z),
    ``explained_raw``, ``perm_sd`` (SD of H_cond_perm), ``z`` = reduction / perm_sd, the numbers of
    rows, stems and X values.
    """
    gcols = _columns(given)
    if by is None:
        raise ValueError("by= is required")
    items = list(by.items()) if isinstance(by, Mapping) else [(c, c) for c in by]
    z0 = codes(dfd[gcols[0]].to_numpy(dtype=object))
    cnt = np.bincount(z0)
    keep = cnt[z0] >= min_count
    sub = dfd.loc[keep]
    y = codes(sub[target].to_numpy(dtype=object))
    z = joint_codes(*[codes(sub[c].to_numpy(dtype=object)) for c in gcols]) if len(sub) else np.zeros(0, np.int64)
    rows = []
    if len(sub) == 0:
        return pd.DataFrame(rows)
    kz = int(z.max()) + 1
    h_z = entropy(z)
    h_given = _pair_entropy(z, kz, y) - h_z
    rng = np.random.default_rng(seed)
    perms = [permute_within(z, rng) for _ in range(n_perm)]
    for name, cols in items:
        x = joint_codes(*[codes(sub[c].to_numpy(dtype=object)) for c in _columns(cols)])
        zx = joint_codes(z, x)
        kzx = int(zx.max()) + 1
        h_cond = _pair_entropy(zx, kzx, y) - entropy(zx)
        sims = []
        for idx in perms:
            zxp = joint_codes(z, x[idx])
            sims.append(_pair_entropy(zxp, int(zxp.max()) + 1, y) - entropy(zxp))
        sims = np.asarray(sims, dtype=float)
        h_perm = float(sims.mean()) if len(sims) else float("nan")
        sd = float(sims.std(ddof=1)) if len(sims) > 1 else float("nan")
        red = h_perm - h_cond
        rows.append({"X": name, "n_tokens": int(len(sub)), "n_stems": int(len(np.unique(z0[keep]))),
                     "n_x_values": int(x.max()) + 1, "H_given": h_given, "H_cond": h_cond, "H_cond_perm": h_perm,
                     "reduction_raw": h_given - h_cond, "bias": h_given - h_perm, "reduction": red,
                     "explained": red / h_given if h_given > 0 else float("nan"),
                     "explained_raw": (h_given - h_cond) / h_given if h_given > 0 else float("nan"),
                     "perm_sd": sd, "z": red / sd if sd and sd > 0 else float("nan")})
    return pd.DataFrame(rows)


def crossval_explained(
    dfd: pd.DataFrame,
    target: str = "suffix",
    given: str = "stem",
    by: Mapping[str, str | Sequence[str]] | Sequence[str] | None = None,
    k: int = 5,
    alpha: float | Sequence[float] = (1.0, 3.0, 10.0, 30.0, 100.0),
    min_count: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Held-out version of :func:`conditional_entropy_decomposition` (sensitivity check).

    Tokens are split into ``k`` random folds.  For every held-out token the cross-entropy (bits) of
    its target is measured under p0(s | z) = (n(z, s) + alpha p(s)) / (n(z) + alpha) and under
    p1(s | z, x) = (n(z, x, s) + alpha p0(s | z)) / (n(z, x) + alpha), both estimated on the other
    folds (``p(s)``: add-one smoothed marginal).  ``explained_cv`` = (CE0 - CE1) / CE0.  Unlike the
    permutation correction this is not conservative for covariates with many values, but it charges
    sparse cells through the smoothing; a covariate without information gives a value <= 0.
    With several ``alpha`` values each cross-entropy uses the value that minimises it (one
    parameter chosen on the held-out folds: slightly optimistic for both terms).
    Returns one row per X with ``CE_given``, ``CE_cond``, ``explained_cv``, ``alpha_given``, ``alpha_cond``.
    """
    alphas = [float(alpha)] if np.isscalar(alpha) else [float(a) for a in alpha]
    gcol = given if isinstance(given, str) else _columns(given)[0]
    items = list(by.items()) if isinstance(by, Mapping) else [(c, c) for c in (by or [])]
    z0 = codes(dfd[gcol].to_numpy(dtype=object))
    keep = np.bincount(z0)[z0] >= min_count
    sub = dfd.loc[keep]
    y = codes(sub[target].to_numpy(dtype=object))
    z = codes(sub[gcol].to_numpy(dtype=object))
    n, ky = len(y), int(y.max()) + 1 if len(y) else 1
    fold = np.random.default_rng(seed).integers(0, k, size=n)
    xs = {name: joint_codes(z, *[codes(sub[c].to_numpy(dtype=object)) for c in _columns(cols)])
          for name, cols in items}
    na = len(alphas)
    ce0 = np.zeros((na, n))
    ce1 = {name: np.zeros((na, na, n)) for name in xs}
    kz = int(z.max()) + 1 if n else 1
    zof = {}
    for name, zx in xs.items():
        zof[name] = np.zeros(int(zx.max()) + 1, dtype=np.int64)
        zof[name][zx] = z  # z of every (z, x) cell (zx is a refinement of z)
    for f in range(k):
        te, tr = fold == f, fold != f
        ps = (np.bincount(y[tr], minlength=ky) + 1.0) / (tr.sum() + ky)
        nzs = np.bincount(z[tr] * ky + y[tr], minlength=kz * ky).reshape(kz, ky).astype(float)
        nz = nzs.sum(axis=1)
        p0s = [(nzs + a * ps) / (nz[:, None] + a) for a in alphas]
        for i, p0 in enumerate(p0s):
            ce0[i, te] = -np.log2(p0[z[te], y[te]])
        for name, zx in xs.items():
            kzx = int(zx.max()) + 1
            nzxs = np.bincount(zx[tr] * ky + y[tr], minlength=kzx * ky).reshape(kzx, ky).astype(float)
            nzx = nzxs.sum(axis=1)
            cell, yy = zx[te], y[te]
            for i, p0 in enumerate(p0s):
                back = p0[zof[name][cell], yy]
                for j, a in enumerate(alphas):
                    ce1[name][i, j, te] = -np.log2((nzxs[cell, yy] + a * back) / (nzx[cell] + a))
    rows = []
    if n == 0:
        return pd.DataFrame(rows)
    m0 = ce0.mean(axis=1)
    i0 = int(np.argmin(m0))
    c0 = float(m0[i0])
    for name in xs:
        m1 = ce1[name].mean(axis=2)
        i, j = np.unravel_index(int(np.argmin(m1)), m1.shape)
        c1 = float(m1[i, j])
        rows.append({"X": name, "CE_given": c0, "CE_cond": c1, "explained_cv": (c0 - c1) / c0 if c0 > 0 else float("nan"),
                     "alpha_given": alphas[i0], "alpha_cond": alphas[j]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# residual structure
# --------------------------------------------------------------------------- #


def neighbour_conditional_mi(
    dfd: pd.DataFrame,
    target: str = "suffix",
    neighbour: str = "suffix",
    offset: int = -1,
    given: str | Sequence[str] = ("stem",),
    unit: str = "cv_line",
    perm_by: str | Sequence[str] | None = None,
    n_perm: int = 100,
    min_count: int = 20,
    seed: int = 0,
) -> dict:
    """I(target_t ; neighbour_{t+offset} | given_t), pairs inside the same ``unit``, permutation-corrected.

    Rows whose first ``given`` column (the stem) has fewer than ``min_count`` tokens are dropped as
    targets (the neighbour may be any token).  The correction permutes the neighbour's value among
    the pairs with the same value of the target row's ``given`` columns (a conditional permutation:
    it keeps the distribution of the neighbour inside each conditioning cell and breaks only its
    association with the target), ``n_perm`` times; ``mi`` = raw - mean(permuted).  ``perm_by``
    (columns of the neighbour's own row) replaces the cells by coarser groups; that null also breaks
    the neighbour's dependence on the target's stem and is biased downwards.  Returns raw, permutation mean and SD, corrected value,
    ``p`` (share of permutations >= raw, +1 smoothing), number of pairs.
    """
    gcols = _columns(given)
    pcols = _columns(perm_by) if perm_by is not None else None
    n = len(dfd)
    if n <= abs(offset):
        return {"mi_raw": float("nan"), "mi_perm": float("nan"), "mi": float("nan"), "perm_sd": float("nan"),
                "p": float("nan"), "n_pairs": 0}
    u = dfd[unit].to_numpy()
    t_idx = np.arange(n)
    nb_idx = t_idx + offset
    ok = (nb_idx >= 0) & (nb_idx < n)
    t_idx, nb_idx = t_idx[ok], nb_idx[ok]
    ok = u[t_idx] == u[nb_idx]
    z0 = codes(dfd[gcols[0]].to_numpy(dtype=object))
    ok &= np.bincount(z0)[z0[t_idx]] >= min_count
    t_idx, nb_idx = t_idx[ok], nb_idx[ok]
    y = codes(dfd[target].to_numpy(dtype=object))[t_idx]
    nb = codes(dfd[neighbour].to_numpy(dtype=object))[nb_idx]
    g = joint_codes(*[codes(dfd[c].to_numpy(dtype=object))[t_idx] for c in gcols])
    pg = g if pcols is None else joint_codes(*[codes(dfd[c].to_numpy(dtype=object))[nb_idx] for c in pcols])
    h_g = conditional_entropy(y, g)

    def cmi(v: np.ndarray) -> float:
        return h_g - conditional_entropy(y, g, v)

    raw = cmi(nb)
    rng = np.random.default_rng(seed)
    sims = np.asarray([cmi(nb[permute_within(pg, rng)]) for _ in range(n_perm)], dtype=float)
    return {"mi_raw": raw, "mi_perm": float(sims.mean()), "mi": raw - float(sims.mean()),
            "perm_sd": float(sims.std(ddof=1)) if n_perm > 1 else float("nan"),
            "p": float((1 + (sims >= raw - 1e-12).sum()) / (1 + n_perm)), "n_pairs": int(len(t_idx))}


def mode_consistency(
    dfd: pd.DataFrame,
    target: str = "suffix",
    stem: str = "stem",
    unit: str = "cv_page",
    groups: str | Sequence[str] | None = None,
    min_count: int = 20,
    min_unit: int = 2,
    n_perm: int = 100,
    seed: int = 0,
) -> dict:
    """Is a stem's affix more uniform inside a page (``unit``) than its tokens spread over the manuscript?

    For stems with >= ``min_count`` tokens and every (stem, unit) cell with >= ``min_unit`` tokens,
    ``within`` = share of tokens that carry the most frequent ``target`` value of their cell
    (token-weighted over cells); ``global`` = share carrying the stem's overall most frequent value
    (same tokens).  Expectation: the unit labels are permuted among the tokens of the same stem (and
    the same ``groups`` value, e.g. position class, language and scribe, so that cells keep their
    composition in those covariates), ``n_perm`` times, and ``within`` recomputed on the same cell
    sizes.  Returns ``within``, ``expected`` (permutation mean), ``excess`` = within - expected,
    ``perm_sd``, ``z``, ``p``, ``global``, numbers of stems, cells and tokens.
    """
    z0 = codes(dfd[stem].to_numpy(dtype=object))
    keep = np.bincount(z0)[z0] >= min_count
    sub = dfd.loc[keep]
    z = codes(sub[stem].to_numpy(dtype=object))
    y = codes(sub[target].to_numpy(dtype=object))
    un = codes(sub[unit].to_numpy(dtype=object))
    gcols = _columns(groups) if groups is not None else []
    perm_groups = joint_codes(z, *[codes(sub[c].to_numpy(dtype=object)) for c in gcols])
    ky = int(y.max()) + 1 if len(y) else 1

    def within(u: np.ndarray) -> tuple[float, int, int]:
        cell = joint_codes(z, u)
        size = np.bincount(cell)
        ok = size[cell] >= min_unit
        if not ok.any():
            return float("nan"), 0, 0
        c2, y2 = pd.factorize(cell[ok])[0], y[ok]
        tab = np.bincount(c2 * ky + y2, minlength=(int(c2.max()) + 1) * ky).reshape(-1, ky)
        return float(tab.max(axis=1).sum() / ok.sum()), int(tab.shape[0]), int(ok.sum())

    w, n_cells, n_tok = within(un)
    gtab = np.bincount(z * ky + y, minlength=(int(z.max()) + 1) * ky).reshape(-1, ky) if len(z) else np.zeros((0, ky))
    mode = gtab.argmax(axis=1) if len(gtab) else np.zeros(0, dtype=np.int64)
    cell = joint_codes(z, un) if len(z) else z
    ok = np.bincount(cell)[cell] >= min_unit if len(z) else np.zeros(0, bool)
    glob = float((y[ok] == mode[z[ok]]).mean()) if ok.any() else float("nan")
    rng = np.random.default_rng(seed)
    sims = np.asarray([within(un[permute_within(perm_groups, rng)])[0] for _ in range(n_perm)], dtype=float)
    exp = float(np.nanmean(sims)) if len(sims) else float("nan")
    sd = float(np.nanstd(sims, ddof=1)) if len(sims) > 1 else float("nan")
    return {"within": w, "expected": exp, "excess": w - exp, "perm_sd": sd,
            "z": (w - exp) / sd if sd and sd > 0 else float("nan"),
            "p": float((1 + (sims >= w - 1e-12).sum()) / (1 + len(sims))), "global": glob,
            "n_stems": int(len(np.unique(z))), "n_cells": n_cells, "n_tokens": n_tok}


def affix_by_covariate(
    dfd: pd.DataFrame, target: str = "suffix", stem: str = "stem", covariate: str = "cv_place", top: int = 10
) -> pd.DataFrame:
    """Long table of p(target | stem, covariate) for the ``top`` most frequent stems.

    Columns: stem, stem rank, covariate value, target value, count, share inside (stem, covariate).
    """
    vc = dfd[stem].value_counts()
    stems = vc.index[:top].tolist()
    sub = dfd[dfd[stem].isin(stems)]
    t = sub.groupby([stem, covariate, target], sort=False).size().rename("count").reset_index()
    t["share"] = t["count"] / t.groupby([stem, covariate])["count"].transform("sum")
    t.insert(1, "stem_rank", t[stem].map({s: i + 1 for i, s in enumerate(stems)}))
    return t.sort_values(["stem_rank", covariate, "count"], ascending=[True, True, False]).reset_index(drop=True)
