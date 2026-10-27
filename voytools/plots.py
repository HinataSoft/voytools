"""Optional plotting helpers (need matplotlib; imported lazily).

The functions build :class:`matplotlib.figure.Figure` objects without ``pyplot``, so they work
on machines without a display; pass ``path`` to save a PNG/PDF/SVG, or show the returned
figure yourself.  Install the extra with ``pip install voytools[plot]``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

__all__ = ["plot_zipf", "plot_positional", "plot_bpe_curve", "have_matplotlib"]


def have_matplotlib() -> bool:
    """True if matplotlib can be imported."""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


def _figure(figsize: tuple[float, float]):
    try:
        from matplotlib.figure import Figure
    except ImportError as exc:  # pragma: no cover - exercised only without matplotlib
        raise ImportError("matplotlib is required for plotting: pip install matplotlib") from exc
    try:
        return Figure(figsize=figsize, layout="constrained")
    except TypeError:  # matplotlib < 3.5
        return Figure(figsize=figsize, constrained_layout=True)


def _finish(fig, path: str | Path | None, dpi: int = 150):
    if path is not None:
        fig.savefig(str(path), dpi=dpi)
    return fig


def plot_zipf(freqs: pd.DataFrame | Mapping[str, pd.DataFrame], path: str | Path | None = None):
    """Rank-frequency (Zipf) plot on log-log axes.

    ``freqs`` is the output of :func:`voytools.stats.token_frequency` or a dict
    ``{label: that output}`` to compare corpora.
    """
    fig = _figure((6.0, 4.2))
    ax = fig.subplots()
    data = freqs if isinstance(freqs, Mapping) else {"tokens": freqs}
    for label, f in data.items():
        counts = np.sort(f["count"].to_numpy(dtype=float))[::-1]
        ax.loglog(np.arange(1, len(counts) + 1), counts, label=label, linewidth=1.2)
    ax.set_xlabel("rank")
    ax.set_ylabel("frequency")
    ax.set_title("Rank-frequency distribution")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    return _finish(fig, path)


def plot_positional(profile: Any, top: int = 8, path: str | Path | None = None):
    """Grouped bars of the most frequent token-initial/final glyphs per position class.

    ``profile`` is a :class:`voytools.stats.PositionalProfile`.
    """
    classes = [c for c in profile.shares.columns if profile.n[c] > 0]
    glyphs = (
        profile.shares[classes].max(axis=1).sort_values(ascending=False, kind="stable").head(top).index.tolist()
    )
    fig = _figure((7.5, 4.2))
    ax = fig.subplots()
    width = 0.8 / max(len(classes), 1)
    x = np.arange(len(glyphs))
    for i, c in enumerate(classes):
        ax.bar(x + i * width, profile.shares.loc[glyphs, c].to_numpy(dtype=float), width, label=c)
    ax.set_xticks(x + width * (len(classes) - 1) / 2)
    ax.set_xticklabels(glyphs)
    ax.set_ylabel("share of tokens in class")
    ax.set_title(f"{profile.which.capitalize()} {profile.unit} of the token by position in the line")
    ax.legend(fontsize=7)
    return _finish(fig, path)


def plot_bpe_curve(curves: pd.DataFrame | Mapping[str, pd.DataFrame], path: str | Path | None = None):
    """Dependence gap ``H(U) - H(U_{t+1}|U_t)`` against the number of BPE merges.

    ``curves`` is the output of :func:`voytools.morphology.bpe_dependence_curve` or a dict of them.
    """
    fig = _figure((6.0, 4.2))
    ax = fig.subplots()
    data = curves if isinstance(curves, Mapping) else {"text": curves}
    for label, c in data.items():
        ax.plot(c["merges"], c["gap"], marker="o", label=label)
    ax.set_xlabel("BPE merges")
    ax.set_ylabel("H(U) - H(U_t+1 | U_t)  [bits]")
    ax.set_title("Dependence gap of within-token BPE units")
    ax.legend()
    ax.grid(True, alpha=0.25)
    return _finish(fig, path)
