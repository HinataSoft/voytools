"""voytools: statistics for Voynich manuscript transliterations in IVTFF 2.0 format.

Quick start::

    from voytools import load_tables, filter_tokens
    from voytools.stats import vocabulary_summary, glyph_entropy_report

    tables = load_tables("ZL3b-n.txt")
    text = filter_tokens(tables.tokens, generic_type="P")      # paragraph text, clean tokens
    print(vocabulary_summary(text, by="language"))
    print(glyph_entropy_report(text, preset="extended"))

Modules: :mod:`~voytools.ivtff` (parser), :mod:`~voytools.table` (DataFrames),
:mod:`~voytools.glyphs` (Eva composites), :mod:`~voytools.stats`, :mod:`~voytools.morphology`,
:mod:`~voytools.nulls`, :mod:`~voytools.selfcitation` (Timm & Schinner generator), :mod:`~voytools.wordgrammar`,
:mod:`~voytools.experiment` (checkpointed many-seed runs),
:mod:`~voytools.plots` (optional, needs matplotlib), :mod:`~voytools.cli`.
"""

__version__ = "0.1.0"

from .experiment import load_runs, merge_runs, run_seeds, seed_range, summarize_runs
from .glyphs import collapse_eva, glyph_sequence, to_symbol_string
from .ivtff import (
    IVTFFError,
    Transliteration,
    parse_ivtff,
    parse_ivtff_text,
    render_high_ascii,
    tokenize_locus,
)
from .table import Tables, build_tables, export, filter_tokens, load_tables

__all__ = [
    "__version__",
    "IVTFFError",
    "Tables",
    "Transliteration",
    "build_tables",
    "collapse_eva",
    "export",
    "filter_tokens",
    "glyph_sequence",
    "load_runs",
    "load_tables",
    "merge_runs",
    "parse_ivtff",
    "parse_ivtff_text",
    "render_high_ascii",
    "run_seeds",
    "seed_range",
    "summarize_runs",
    "to_symbol_string",
    "tokenize_locus",
]
