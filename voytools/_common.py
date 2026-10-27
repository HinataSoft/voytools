"""Internal helpers shared by the analysis modules (not part of the public API)."""

from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from .glyphs import collapse_eva, default_preset

NA_LABEL = "NA"


# --------------------------------------------------------------------------- #
# presets
# --------------------------------------------------------------------------- #


def resolve_preset(df: pd.DataFrame | None, preset: Any) -> Any:
    """Resolve ``preset='auto'`` from ``df.attrs`` (``glyph_preset`` or ``alphabet``).

    Falls back to ``'basic'`` when the table carries no alphabet information (for
    example after a round trip through CSV).  Pass ``preset=None`` explicitly for
    non-Eva data in that case.
    """
    if not (isinstance(preset, str) and preset == "auto"):
        return preset
    attrs = getattr(df, "attrs", {}) or {}
    if "glyph_preset" in attrs:
        return attrs["glyph_preset"]
    if attrs.get("alphabet"):
        return default_preset(attrs["alphabet"])
    return "basic"


def glyph_lists(df: pd.DataFrame, preset: Any = "auto", column: str = "token") -> list[tuple[str, ...]]:
    """Glyph tuple for every row of ``df[column]`` (cached per distinct token)."""
    preset = resolve_preset(df, preset)
    tokens = df[column].tolist()
    cache: dict[str, tuple[str, ...]] = {}
    out: list[tuple[str, ...]] = []
    for t in tokens:
        g = cache.get(t)
        if g is None:
            g = collapse_eva(t, preset)
            cache[t] = g
        out.append(g)
    return out


# --------------------------------------------------------------------------- #
# lines, pages, groups
# --------------------------------------------------------------------------- #


def _key_strings(df: pd.DataFrame, cols: Sequence[str]) -> list[str]:
    parts = []
    for c in cols:
        if c in df.columns:
            parts.append(["" if pd.isna(v) else str(v) for v in df[c].tolist()])
    if not parts:
        raise KeyError(f"none of the columns {list(cols)} found in the table")
    if len(parts) == 1:
        return parts[0]
    return ["|".join(items) for items in zip(*parts)]


def _run_ids(keys: list[str]) -> np.ndarray:
    n = len(keys)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    arr = np.asarray(keys, dtype=object)
    change = np.empty(n, dtype=bool)
    change[0] = True
    change[1:] = arr[1:] != arr[:-1]
    return np.cumsum(change) - 1


def line_ids(df: pd.DataFrame) -> np.ndarray:
    """Integer id per row; consecutive rows of the same locus (and transcriber) share an id."""
    if "locus_id" in df.columns:
        return _run_ids(_key_strings(df, ["folio", "locus_id", "transcriber"]))
    if "line_id" in df.columns:
        return _run_ids([str(v) for v in df["line_id"].tolist()])
    raise KeyError("token table needs a 'locus_id' (or 'line_id') column to recognise lines")


def page_ids(df: pd.DataFrame) -> np.ndarray:
    """Integer id per row; consecutive rows of the same page (and transcriber) share an id."""
    if "folio" in df.columns:
        return _run_ids(_key_strings(df, ["folio", "transcriber"]))
    raise KeyError("token table needs a 'folio' column to recognise pages")


def by_columns(by: str | Sequence[str] | None) -> list[str]:
    """Normalise a ``by=`` argument to a list of column names."""
    if by is None:
        return []
    if isinstance(by, str):
        return [by]
    return list(by)


def iter_groups(df: pd.DataFrame, by: str | Sequence[str] | None):
    """Yield ``(key_tuple, sub_dataframe)`` for every group of ``by`` (missing values form group 'NA').

    Groups are ordered by their key strings (``'NA'`` last).  With ``by=None`` a single
    group with the empty key is yielded.  The sub-frames are copies.
    """
    cols = by_columns(by)
    if not cols:
        yield (), df
        return
    for c in cols:
        if c not in df.columns:
            raise KeyError(f"grouping column {c!r} not in table")
    keys = list(
        zip(*[[NA_LABEL if pd.isna(v) else str(v) for v in df[c].tolist()] for c in cols])
    )
    buckets: dict[tuple[str, ...], list[int]] = {}
    for i, k in enumerate(keys):
        buckets.setdefault(k, []).append(i)
    for k in sorted(buckets, key=lambda kk: tuple((x == NA_LABEL, x) for x in kk)):
        yield k, df.iloc[buckets[k]].copy()


def assemble(rows: list[tuple[tuple[str, ...], dict[str, Any]]], by: str | Sequence[str] | None) -> pd.DataFrame:
    """Build the result frame of a by-capable statistic from ``(key, metrics)`` rows."""
    cols = by_columns(by)
    if not cols:
        return pd.DataFrame([m for _, m in rows], index=["all"] * len(rows))
    records = []
    for key, metrics in rows:
        rec = dict(zip(cols, key))
        rec.update(metrics)
        records.append(rec)
    return pd.DataFrame(records)


def is_nested(seq: Any) -> bool:
    """True if ``seq`` is a non-empty sequence of lists/tuples (several segments)."""
    if isinstance(seq, (str, bytes)):
        return False
    for first in seq:
        return isinstance(first, (list, tuple))
    return False


def segments(seq: Any) -> list[Sequence[Any]]:
    """Return ``seq`` as a list of segments (a flat sequence becomes one segment)."""
    if is_nested(seq):
        return [s for s in seq]
    return [list(seq) if not isinstance(seq, (list, tuple, str)) else seq]


def as_str_array(values: Iterable[Any]) -> np.ndarray:
    """Object array of strings (missing -> '')."""
    return np.asarray(["" if (v is None or (isinstance(v, float) and v != v)) else str(v) for v in values], dtype=object)
