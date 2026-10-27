"""Command line interface: ``python -m voytools <command> FILE [options]``.

Commands
--------
``summary``    pages, loci by type, tokens, types, hapax, lengths, glyph entropy, breakdowns
``export``     write the token table to CSV (or Parquet if pyarrow is installed)
``freq``       most frequent tokens
``entropy``    glyph entropies H(X), H(X_{t+1}|X_t) for the Eva presets
``positions``  token-initial/final glyphs by position in the line, Jensen-Shannon divergences
``affixes``    productive prefixes / suffixes
``mi``         shuffle-corrected adjacent-token and edge-glyph mutual information, neighbour similarity

Common options select tokens (``--generic-type``, ``--language``, ``--illustration``,
``--hand``, ``--transcriber``, ``--clean/--no-clean``) and control the tokenizer
(``--space-policy``, ``--alternatives``).  The statistical commands default to paragraph
text (``--generic-type P``) and clean tokens; ``summary`` and ``export`` default to all loci.
Output is plain text tables; ``--csv`` prints CSV instead.
"""

from __future__ import annotations

import argparse
import sys
from typing import Any, Sequence

import pandas as pd

from . import __version__
from .glyphs import is_eva_alphabet
from .ivtff import render_high_ascii
from .table import Tables, build_tables, export as export_table, filter_tokens, most_complete_transcriber
from .ivtff import parse_ivtff

_STAT_COMMANDS = {"freq", "entropy", "positions", "affixes", "mi"}


# --------------------------------------------------------------------------- #
# output helpers
# --------------------------------------------------------------------------- #


def _render_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        vals = out[c].tolist()
        if any(isinstance(v, str) for v in vals):
            out[c] = [render_high_ascii(v) if isinstance(v, str) else v for v in vals]
    ix = out.index.tolist()
    if any(isinstance(v, str) for v in ix):
        out.index = [render_high_ascii(v) if isinstance(v, str) else v for v in ix]
    return out


def _emit(df: pd.DataFrame, args: argparse.Namespace, title: str | None = None, index: bool = False) -> None:
    df = _render_frame(df)
    if args.csv:
        if title:
            print(f"# {title}")
        print(df.to_csv(index=index, float_format="%.6g"), end="")
        print()
        return
    if title:
        print(f"\n== {title} ==")
    print(df.to_string(index=index, float_format=lambda x: f"{x:.4g}" if abs(x) < 1e-3 and x != 0 else f"{x:.4f}"))


def _note(args: argparse.Namespace, text: str) -> None:
    if args.csv:
        print(f"# {text}")
    else:
        print(text)


def _warn(text: str) -> None:
    print(f"voytools: {text}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# argument parsing
# --------------------------------------------------------------------------- #


def _common_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group("token selection")
    g.add_argument("--generic-type", default=None, metavar="T",
                   help="locus generic type(s), e.g. P, 'P,L' or 'all' (default: P for statistics, all for summary/export)")
    g.add_argument("--language", default=None, metavar="A|B", help="Currier language ($L) of the page")
    g.add_argument("--illustration", default=None, metavar="H", help="illustration type ($I): A B C H P S T Z")
    g.add_argument("--hand", default=None, metavar="1", help="scribe hand ($H): 1-5")
    g.add_argument("--transcriber", default=None, metavar="X",
                   help="transcriber ID for interlinear files (default: the most complete one, with a notice)")
    g.add_argument("--clean", action=argparse.BooleanOptionalAction, default=None,
                   help="drop tokens with unreadable characters, alternative readings or high-ASCII characters "
                        "(default: on for statistics, off for export)")
    t = p.add_argument_group("tokenizer")
    t.add_argument("--space-policy", choices=["split", "join", "strict"], default="split",
                   help="how the uncertain space ',' is treated (default: split)")
    t.add_argument("--alternatives", choices=["first", "mark"], default="first",
                   help="how [a:b] readings are treated (default: first)")
    t.add_argument("--infer-paragraphs", action="store_true",
                   help="infer paragraphs for files without <%%>/<$> marks")
    o = p.add_argument_group("output")
    o.add_argument("--csv", action="store_true", help="print CSV instead of aligned text")
    return p


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser (exposed for tests)."""
    common = _common_parser()
    ap = argparse.ArgumentParser(
        prog="python -m voytools",
        description="Statistics for Voynich manuscript transliterations in IVTFF format.",
    )
    ap.add_argument("--version", action="version", version=f"voytools {__version__}")
    sub = ap.add_subparsers(dest="command", required=True, metavar="command")

    s = sub.add_parser("summary", parents=[common], help="overview of a transliteration file")
    s.add_argument("file")

    s = sub.add_parser("export", parents=[common], help="write the token table to CSV/Parquet")
    s.add_argument("file")
    s.add_argument("out", help="output path (.csv or .parquet)")
    s.add_argument("--tables", action="store_true", help="also write OUT_loci.csv and OUT_pages.csv")

    s = sub.add_parser("freq", parents=[common], help="token frequencies")
    s.add_argument("file")
    s.add_argument("--top", type=int, default=50, help="number of tokens per group (default 50)")
    s.add_argument("--by", default=None, help="group by this column, e.g. language")

    s = sub.add_parser("entropy", parents=[common], help="glyph entropies")
    s.add_argument("file")
    s.add_argument("--preset", choices=["basic", "extended", "none", "auto"], default=None,
                   help="Eva composite preset (default: basic and extended for Eva, none otherwise)")
    s.add_argument("--by", default=None, help="group by this column, e.g. language")
    s.add_argument("--boundary", default=None, metavar="SYM", help="insert this symbol between tokens")
    s.add_argument("--scope", choices=["line", "all"], default="line",
                   help="'line': bigrams do not cross line breaks; 'all': one stream for the whole text")

    s = sub.add_parser("positions", parents=[common], help="token-initial/final glyphs by line position")
    s.add_argument("file")
    s.add_argument("--which", choices=["first", "last", "both"], default="both")
    s.add_argument("--preset", choices=["basic", "extended", "none", "auto"], default="auto")
    s.add_argument("--top", type=int, default=8, help="glyphs shown per class (default 8)")
    s.add_argument("--permutations", type=int, default=0, help="permutation test for the JSD matrix (default 0 = off)")
    s.add_argument("--seed", type=int, default=0)

    s = sub.add_parser("affixes", parents=[common], help="productive prefixes and suffixes")
    s.add_argument("file")
    s.add_argument("--side", choices=["prefix", "suffix", "both"], default="both")
    s.add_argument("--max-len", type=int, default=3, help="maximum affix length in glyphs (default 3)")
    s.add_argument("--min-count", type=int, default=20)
    s.add_argument("--top", type=int, default=30, help="rows shown (default 30)")
    s.add_argument("--preset", choices=["basic", "extended", "none", "auto"], default="auto")
    s.add_argument("--by", default=None)

    s = sub.add_parser("mi", parents=[common], help="adjacent-token / edge-glyph MI and neighbour similarity")
    s.add_argument("file")
    s.add_argument("--shuffles", type=int, default=20)
    s.add_argument("--cap", type=int, default=2000, help="token types kept for the token MI (default 2000)")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--preset", choices=["basic", "extended", "none", "auto"], default="auto")
    s.add_argument("--window", type=int, default=5, help="neighbour-similarity window (default 5)")
    s.add_argument("--no-neighbor", action="store_true", help="skip the neighbour-similarity test")
    s.add_argument("--by", default=None)
    return ap


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #


def _preset_arg(name: str | None) -> Any:
    if name is None or name == "auto":
        return "auto"
    return None if name == "none" else name


def _load(args: argparse.Namespace) -> tuple[Tables, pd.DataFrame]:
    translit = parse_ivtff(args.file)
    stat = args.command in _STAT_COMMANDS
    transcriber = args.transcriber
    if transcriber is None and translit.is_interlinear and args.command != "export":
        transcriber = most_complete_transcriber(translit)
        _warn(
            f"interlinear file with {len(translit.transcribers)} transcribers "
            f"({''.join(translit.transcribers)}); using '{transcriber}' (most loci). Select with --transcriber."
        )
    tables = build_tables(
        translit,
        transcriber=transcriber,
        space_policy=args.space_policy,
        alternatives=args.alternatives,
        infer_paragraphs=args.infer_paragraphs,
    )
    generic = args.generic_type
    if generic is None and stat:
        generic = "P"
    clean = args.clean if args.clean is not None else stat or args.command == "summary"
    tokens = filter_tokens(
        tables.tokens,
        generic_type=generic,
        language=args.language,
        illustration=args.illustration,
        hand=args.hand,
        clean=clean,
    )
    if not len(tokens):
        _warn("no tokens left after filtering")
    return tables, tokens


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def _cmd_summary(args: argparse.Namespace) -> int:
    from .stats import glyph_entropy_report, length_distribution, vocabulary_summary

    tables, tokens = _load(args)
    h = tables.header
    _note(args, f"File: {args.file}")
    _note(args, f"Alphabet: {h.alphabet}   IVTFF version: {h.version}   origin: {h.origin}")
    if tables.transcribers:
        _note(args, f"Transcribers: {' '.join(tables.transcribers)}")
    n_loci = len(tables.loci)
    _note(args, f"Pages: {len(tables.pages)}   Loci: {n_loci}   Tokens in loci: {int(tables.loci['n_tokens'].sum())}")
    by_type = tables.loci.groupby("generic_type", sort=True).size().rename("loci").reset_index()
    _emit(by_type, args, "Loci by generic type")
    by_sub = tables.loci.groupby("locus_type", sort=True).size().rename("loci").reset_index()
    _emit(by_sub, args, "Loci by type and subtype")
    sel = filter_desc(args)
    _note(args, f"\nSelected tokens ({sel}): {len(tokens)}")
    if not len(tokens):
        return 0
    vs = vocabulary_summary(tokens)
    _emit(vs.reset_index(drop=True), args, "Vocabulary")
    for col, title in (("generic_type", "by generic type"), ("language", "by language"),
                       ("illustration", "by illustration"), ("hand", "by hand")):
        if tokens[col].nunique() > 1:
            _emit(vocabulary_summary(tokens, by=col), args, f"Vocabulary {title}")
    ld = length_distribution(tokens, unit="chars")
    lens = tokens["len_chars"]
    lg = tokens["len_glyphs"]
    stats = pd.DataFrame(
        {"mean": [lens.mean(), lg.mean()], "median": [lens.median(), lg.median()],
         "sd": [lens.std(), lg.std()], "min": [lens.min(), lg.min()], "max": [lens.max(), lg.max()]},
        index=["chars", "glyphs"],
    )
    _emit(stats.reset_index().rename(columns={"index": "unit"}), args, "Token length")
    _emit(ld.head(15), args, "Length distribution (characters, first rows)")
    presets = ["basic", "extended"] if is_eva_alphabet(h.alphabet) else [None]
    frames = []
    for pr in presets:
        r = glyph_entropy_report(tokens, pr)
        r.insert(0, "preset", pr or "none")
        frames.append(r.reset_index(drop=True))
    _emit(pd.concat(frames, ignore_index=True), args, "Glyph entropy (spaces removed, bits)")
    if tables.warnings:
        _note(args, f"\nParser warnings: {len(tables.warnings)} (first: {tables.warnings[0]})")
    return 0


def filter_desc(args: argparse.Namespace) -> str:
    """Describe the active token selection for the summary header."""
    parts = []
    for name in ("generic_type", "language", "illustration", "hand", "transcriber"):
        v = getattr(args, name, None)
        if v:
            parts.append(f"{name}={v}")
    clean = args.clean if args.clean is not None else True
    parts.append("clean" if clean else "unclean")
    return ", ".join(parts)


def _cmd_export(args: argparse.Namespace) -> int:
    tables, tokens = _load(args)
    path = export_table(tokens, args.out)
    print(f"wrote {len(tokens)} tokens to {path}")
    if args.tables:
        stem = str(path.with_suffix(""))
        for name in ("loci", "pages"):
            p = export_table(getattr(tables, name), f"{stem}_{name}.csv")
            print(f"wrote {name} to {p}")
    return 0


def _cmd_freq(args: argparse.Namespace) -> int:
    from .stats import token_frequency

    _, tokens = _load(args)
    f = token_frequency(tokens, by=args.by)
    if args.by:
        f = f.groupby(args.by, sort=False).head(args.top)
    else:
        f = f.head(args.top)
    _emit(f, args, f"Top {args.top} tokens")
    return 0


def _cmd_entropy(args: argparse.Namespace) -> int:
    from .stats import glyph_entropy_report

    tables, tokens = _load(args)
    if args.preset is None:
        presets: list[Any] = ["basic", "extended"] if is_eva_alphabet(tables.header.alphabet) else [None]
    else:
        presets = [_preset_arg(args.preset)]
    frames = []
    for pr in presets:
        r = glyph_entropy_report(tokens, pr, by=args.by, boundary=args.boundary, scope=args.scope)
        r.insert(0, "preset", pr if isinstance(pr, str) and pr != "auto" else ("none" if pr is None else pr))
        frames.append(r.reset_index(drop=True))
    _emit(pd.concat(frames, ignore_index=True), args, "Glyph entropy (bits): H1 = H(X), H2 = H(X_t+1 | X_t)")
    return 0


def _cmd_positions(args: argparse.Namespace) -> int:
    from .stats import positional_profile

    _, tokens = _load(args)
    preset = _preset_arg(args.preset)
    for which in (["first", "last"] if args.which == "both" else [args.which]):
        prof = positional_profile(tokens, which=which, preset=preset, permutations=args.permutations, seed=args.seed)
        _emit(prof.n.rename("tokens").reset_index().rename(columns={"index": "class"}), args,
              f"Position classes ({which} glyph of the token)")
        rows = []
        for cls in prof.shares.columns:
            if prof.n[cls] == 0:
                continue
            top = prof.top(cls, args.top)
            rows.append({"class": cls, **{f"#{i + 1}": f"{g} {v * 100:.1f}%" for i, (g, v) in enumerate(top.items())}})
        _emit(pd.DataFrame(rows), args, f"Most frequent {which} glyphs per class")
        _emit(prof.jsd.round(4).reset_index().rename(columns={"index": "class"}), args,
              f"Jensen-Shannon divergence between classes ({which} glyph, bits)")
        if prof.pvalues is not None:
            _emit(prof.pvalues.round(4).reset_index().rename(columns={"index": "class"}), args,
                  f"Permutation p-values ({args.permutations} permutations)")
    return 0


def _cmd_affixes(args: argparse.Namespace) -> int:
    from .morphology import affix_inventory

    _, tokens = _load(args)
    preset = _preset_arg(args.preset)
    for side in (["prefix", "suffix"] if args.side == "both" else [args.side]):
        inv = affix_inventory(tokens, side, args.max_len, args.min_count, preset, by=args.by)
        if args.by:
            inv = inv.groupby(args.by, sort=False).head(args.top)
        else:
            inv = inv.head(args.top)
        _emit(inv, args, f"{'Prefixes' if side == 'prefix' else 'Suffixes'} (top {args.top} by count)")
    return 0


def _cmd_mi(args: argparse.Namespace) -> int:
    from .stats import adjacent_token_mi, edge_glyph_mi, neighbor_similarity

    _, tokens = _load(args)
    preset = _preset_arg(args.preset)
    _emit(adjacent_token_mi(tokens, cap=args.cap, shuffles=args.shuffles, seed=args.seed, by=args.by),
          args, f"Adjacent-token MI (bits; within lines; {args.shuffles} within-line shuffles)", index=args.by is None)
    _emit(edge_glyph_mi(tokens, shuffles=args.shuffles, seed=args.seed, preset=preset, by=args.by),
          args, "Edge-glyph MI: last glyph of token t vs first glyph of token t+1 (bits)", index=args.by is None)
    if not args.no_neighbor:
        _emit(neighbor_similarity(tokens, window=args.window, shuffles=args.shuffles, seed=args.seed, by=args.by),
              args, f"Neighbour similarity (window {args.window}, within-page shuffles)", index=args.by is None)
    return 0


_COMMANDS = {
    "summary": _cmd_summary,
    "export": _cmd_export,
    "freq": _cmd_freq,
    "entropy": _cmd_entropy,
    "positions": _cmd_positions,
    "affixes": _cmd_affixes,
    "mi": _cmd_mi,
}


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point; returns the process exit code."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")  # type: ignore[attr-defined]
        except (AttributeError, ValueError):
            pass
    args = build_parser().parse_args(argv)
    try:
        return _COMMANDS[args.command](args)
    except (OSError, ValueError, KeyError) as exc:
        print(f"voytools: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
