"""Experiment 020: paragraph and page as units of content (stem sharing between lines).

Usage (from the repository root; see ``spec.md``)::

    python experiments/020_paragraph_as_unit/run.py run --combo ZL_split                  # all sources
    python experiments/020_paragraph_as_unit/run.py run --combo GC_join --sources real,perm_lines
    python experiments/020_paragraph_as_unit/run.py run --combo ctl:celsus_chap             # a control
    python experiments/020_paragraph_as_unit/run.py summarize

Combinations: ``{ZL, GC}_{split, join}`` (real text, ``P`` clean) and ``ctl:<name>`` controls in the ZL split
layout (lines and pages of ZL):

* ``celsus_chap`` / ``pliny_chap`` – Latin prose with its **own paragraphs** = chapters of the source
  (``controls.tess_paragraphs``, ``controls.paragraph_control``): a chapter's end closes the line, the next
  chapter starts a new line; sections = books;
* ``celsus_sent40`` / ``pliny_sent40`` – the same prose with weak paragraphs: sentence groups of >= 40 tokens
  inside a chapter (the size of a Voynich paragraph);
* ``linnaeus_par`` – Linnaeus, *Species Plantarum*, item = paragraph (items < 3 tokens merged; positive control);
* ``pliny_prose`` – Pliny as running prose in the ZL layout **with ZL's paragraph boundaries** (009; the
  boundaries mean nothing for the text: negative control of the expectation);
* ``gibberish_pooled`` – human gibberish (013) streamed into the 60-page ZL crop with ZL's boundaries (the
  transcriptions carry no paragraph marks; negative control); ``gibberish_native`` – the writers' own lines and
  pages, one paragraph per page (page-boundary metric only).

Sources (``seed`` = repetition): ``real`` (seed 0); ``real_boot`` – bootstrap over pages (200); ``shuffle_bounds``
– null (a): paragraph lengths of every page permuted (200); ``perm_lines`` – null (b): lines permuted on the page
(200); ``selfcitation_21`` – ts2020 config 21 (004), 30 samples, Eva only; ``template_005`` – the 005 template
with topic, 30 samples.  Controls: ``real``, ``real_boot``, ``shuffle_bounds``, ``perm_lines``.

Variants (column prefix): ``st`` stems with >= 2 occurrences (decomposition (a), learned affixes), ``st50``
the same without the 50 most frequent stems, ``tok`` whole tokens with >= 2 occurrences, ``ng`` stems with
leading gallows glyphs stripped (artefact check of the first-line metric).  Metrics:
``voytools.paragraphs.all_metrics`` (``<variant>__<metric>__<key>``).
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_tables, run_seeds  # noqa: E402
from voytools import controls as C  # noqa: E402
from voytools import paragraphs as P  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.selfcitation import selfcitation_like  # noqa: E402
from voytools.templates import ParadigmTemplate  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M6 = _load("exp006", "006_topic_at_distance")
M9 = _load("exp009", "009_register_profile")

DATA = ROOT / "data"
OUT = ROOT / "results" / "020_paragraph_as_unit"
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
COMBOS = ["ZL_split", "ZL_join", "GC_split", "GC_join"]
CONTROLS = ("celsus_chap", "pliny_chap", "celsus_sent40", "pliny_sent40", "linnaeus_par", "pliny_prose",
            "gibberish_pooled", "gibberish_native")
SENT_TARGET = 40
LINNAEUS_MIN_ITEM = 3
VARIANTS = {"st": dict(unit="stem", min_count=2), "st50": dict(unit="stem", min_count=2, drop_top=50),
            "tok": dict(unit="token", min_count=2), "ng": dict(unit="stem", min_count=2, gallows=True)}
SOURCES = {"real": range(1), "real_boot": range(200), "shuffle_bounds": range(200), "perm_lines": range(200),
           "selfcitation_21": range(30), "template_005": range(30)}
CONTROL_SOURCES = ("real", "real_boot", "shuffle_bounds", "perm_lines")
TS_SEED0, TPL_SEED0 = 2100, 100_000


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


def _tess_divisions(pattern: str, target: int | None) -> list:
    paths = sorted(M9._vu("").glob(pattern), key=lambda p: (len(p.name), p.name))
    return [d for p in paths for d in C.tess_paragraphs(C.read_text(p), level=1, target=target)], paths


def control_table(name: str) -> pd.DataFrame:
    """Token table of a 020 control (``attrs['control']`` describes the layout)."""
    tpl = M9.real_table("ZL_split")[0]
    if name in ("celsus_chap", "pliny_chap", "celsus_sent40", "pliny_sent40"):
        pattern = "latin/celsus*.tess" if name.startswith("celsus") else "herbal/pliny*.tess"
        divs, paths = _tess_divisions(pattern, SENT_TARGET if name.endswith("sent40") else None)
        t = C.paragraph_control(divs, tpl, name)
    elif name == "linnaeus_par":
        divs, paths = M9._linnaeus()
        t = C.paragraph_control(divs, tpl, name, min_item=LINNAEUS_MIN_ITEM)
    elif name == "pliny_prose":
        return M9.control_table(name)
    elif name in ("gibberish_pooled", "gibberish_native"):
        M13 = _load("exp013", "013_gibberish_and_cipher_controls")
        return M13.build(name)
    else:
        raise ValueError(name)
    t.attrs["control"]["files"] = [M9.rel_source(p) for p in paths]
    return t


class Context:
    """One corpus: token table with stems, page order, the end-mark flag of pages, line tables of the real text."""

    def __init__(self, combo: str):
        self.combo = combo
        self.is_control = combo.startswith("ctl:")
        if self.is_control:
            name = combo[4:]
            self.key, self.policy, self.file = "ZL", "split", name
            self.real = control_table(name).reset_index(drop=True)
            self.rank = {p: i for i, p in enumerate(pd.unique(self.real["folio"].to_numpy()))}
            if "par_continues" in self.real.columns:
                self.open_end = self.real.groupby("folio", sort=False)["par_continues"].first().astype(bool).to_dict()
            else:
                self.open_end = None
        else:
            self.key, self.policy = combo.split("_")
            self.file = FILES[self.key]
            tables = load_tables(DATA / self.file, space_policy=self.policy)
            self.rank = {p: i for i, p in enumerate(pd.unique(tables.tokens["folio"].to_numpy()))}
            self.real = filter_tokens(tables.tokens, generic_type="P", clean=True).reset_index(drop=True)
            loci = tables.loci[tables.loci["generic_type"] == "P"]
            last = loci.groupby("folio", sort=False).tail(1)
            self.open_end = dict(zip(last["folio"], ~last["paragraph_end"].astype(bool)))
        self.real_d = decompose_table(self.real)
        self.prefixes, self.suffixes = self.real_d.attrs["prefixes"], self.real_d.attrs["suffixes"]
        self.pages = list(pd.unique(self.real["folio"].to_numpy()))
        self._tpl = None

    def columns(self) -> dict:
        return {"combo": self.combo, "file": self.file, "space_policy": self.policy, "n_pages": len(self.pages),
                "n_lines": int(self.real["locus_id"].nunique())}

    def lines(self, df_d: pd.DataFrame) -> dict:
        """``{variant: (line table, content pairs)}`` of a decomposed sample."""
        out = {}
        for var, kw in VARIANTS.items():
            L = P.line_table(df_d, rank=self.rank, **kw)
            out[var] = (L, P.pair_content(L))
        return out

    def decompose(self, sample: pd.DataFrame) -> pd.DataFrame:
        return decompose_table(sample, prefixes=self.prefixes, suffixes=self.suffixes)

    def template(self) -> ParadigmTemplate:
        if self._tpl is None:
            self._tpl = ParadigmTemplate.learn(self.real_d, **M6.template_config(self.key))
        return self._tpl


def metrics(ctx: Context, lines: dict, weights: dict | None = None, relabel=None) -> dict:
    """All metrics of all variants; ``relabel(L, seed)`` -> new line table (layout nulls)."""
    out = {}
    for var, (L, pc) in lines.items():
        Lx = L if relabel is None else relabel(L)
        pairs = P.pair_layout(pc, Lx)
        out.update(P.all_metrics(Lx, pairs, weights=weights, open_end=ctx.open_end, prefix=f"{var}__"))
        out[f"{var}__n_units"] = L.attrs["n_units"]
        out[f"{var}__n_kept"] = L.attrs["n_kept"]
    return out


def page_weights(ctx: Context, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    draw = rng.integers(0, len(ctx.pages), size=len(ctx.pages))
    cnt = np.bincount(draw, minlength=len(ctx.pages))
    return dict(zip(ctx.pages, cnt.astype(float)))


def make_runner(ctx: Context, source: str, real_lines: dict):
    def one(seed: int) -> dict:
        base = {**ctx.columns(), "source": source}
        if source == "real":
            return {**base, **metrics(ctx, real_lines)}
        if source == "real_boot":
            return {**base, **metrics(ctx, real_lines, weights=page_weights(ctx, seed))}
        if source == "shuffle_bounds":
            return {**base, **metrics(ctx, real_lines, relabel=lambda L: P.shuffle_boundaries(L, seed))}
        if source == "perm_lines":
            return {**base, **metrics(ctx, real_lines, relabel=lambda L: P.permute_lines(L, seed))}
        if source == "selfcitation_21":
            if ctx.key != "ZL" or ctx.is_control:
                raise ValueError("selfcitation_like generates basic Eva; no null for a non-Eva file")
            sample = selfcitation_like(ctx.real, seed=TS_SEED0 + seed, **M6.selfcitation_params())
        elif source == "template_005":
            sample = ctx.template().sample(ctx.real, seed=TPL_SEED0 + seed, pool_seed=seed)
        else:
            raise ValueError(source)
        return {**base, **metrics(ctx, ctx.lines(ctx.decompose(sample)))}

    return one


def progress(i: int, n: int, s: int) -> None:
    if i % 25 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_run(args) -> None:
    default = CONTROL_SOURCES if args.combo.startswith("ctl:") else list(SOURCES)
    sources = list(default) if args.sources == "all" else args.sources.split(",")
    ctx = Context(args.combo)
    real_lines = ctx.lines(ctx.real_d)
    parts = OUT / "runs" / args.combo.replace(":", "_")
    parts.mkdir(parents=True, exist_ok=True)
    L = real_lines["st"][0]
    L[P.LINE_COLS + ["n", "locus_id"]].to_csv(parts / "lines.csv", index=False)
    for source in sources:
        seeds = M6.parse_seeds(args.seeds) if args.seeds else list(SOURCES[source])
        print(f"{args.combo} {source}: {len(seeds)} seeds", flush=True)
        params = {"combo": args.combo, "source": source, "variants": VARIANTS, "n_pages": len(ctx.pages),
                  "n_lines": int(ctx.real["locus_id"].nunique()),
                  "n_pars": int(L.groupby("folio")["par"].nunique().sum()),
                  "n_pages_2pars": int((L.groupby("folio")["n_pars"].first() >= 2).sum()),
                  "n_open_end": int(sum(ctx.open_end.values())) if ctx.open_end else 0}
        if source == "selfcitation_21":
            params["selfcitation_004"] = {"row": M6.TS_ROW, **{k: (list(v) if isinstance(v, tuple) else v)
                                                               for k, v in M6.selfcitation_params().items()}}
            params["seed_offset"] = TS_SEED0
        if source == "template_005":
            params["template_005"] = M6.template_config(ctx.key)
            params["seed_offset"], params["pool_seed"] = TPL_SEED0, "repetition"
        if ctx.is_control:
            params["control"] = {k: v for k, v in ctx.real.attrs.get("control", {}).items()}
        run_seeds(f"020_{args.combo}_{source}", seeds, make_runner(ctx, source, real_lines),
                  parts / source, params=params, progress=progress)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--combo", required=True)
    p.add_argument("--sources", default="all")
    p.add_argument("--seeds", default=None, help="'0:30', '0:200/4:1' or '3,5'")
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("summarize")
    p.set_defaults(func=lambda a: __import__("summarize").summarize())
    args = ap.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    args.func(args)


if __name__ == "__main__":
    main()
