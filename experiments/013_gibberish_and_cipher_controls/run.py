"""Experiment 013: missing reference points - human gibberish and the cipher of a real text.

Sources (built with ``voytools.controls``):

* ``gibberish_pooled`` - the 38 writers of Gaskell & Bowern (2022), cleaned like their ``stats.py``,
  concatenated in a random writer order (seed 13) and set as running text into the ZL split layout
  cropped to the size of the gibberish (evenly spaced pages, ``controls.crop_template``); writer =
  section = scribe (majority per page);
* ``gibberish_native`` - the same text in the writers' own lines and pages (for line metrics);
* ``voynich_crop`` - ZL split restricted to the same cropped layout (Voynich at the gibberish size);
  ``caesar_crop``, ``celsus_crop`` - the prose controls of 009 in the cropped layout;
* ``naibbe_caesar`` (Naibbe ciphertext of Caesar, R&T driver), ``naibbe_pliny`` (Greshko's ciphertext of
  Pliny XVI) and ``naibbe_plaintext`` (Greshko's respaced plaintext, aligned chunk-for-chunk) - cropped to
  the ZL split size, running text in the full ZL layout, sections = 7 contiguous blocks (books are not
  marked in the files).  Ciphertext is Eva (glyph units as for Voynich); plaintext is letters.

Metrics: the register profile of 009 (``M9.profile``), the line metrics of 010 (``M10.line_profile``)
and the decomposition of 011 (``M11.measure``), with the same code; per writer: vocabulary, H2, word
length, suffixes per stem, line-end length, compared with Voynich chunks of the same sizes.

Usage (from the repository root)::

    python experiments/013_gibberish_and_cipher_controls/run.py sources
    python experiments/013_gibberish_and_cipher_controls/run.py profile --source gibberish_pooled   # 009
    python experiments/013_gibberish_and_cipher_controls/run.py lines --source gibberish_native     # 010
    python experiments/013_gibberish_and_cipher_controls/run.py variation --source naibbe_pliny     # 011
    python experiments/013_gibberish_and_cipher_controls/run.py writers
    python experiments/013_gibberish_and_cipher_controls/run.py summarize

Seeds as in 009: ``-1`` = full text, ``0..99`` = subsamples of 80 % of the pages (011: ``0..29``).
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import controls as C  # noqa: E402
from voytools import run_seeds  # noqa: E402
from voytools._common import glyph_lists  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.stats import conditional_entropy, entropy, vocabulary_summary  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M11 = _load("exp011", "011_suffix_variation")
M9 = M11.M9
M10 = _load("exp010", "010_line_as_item")

OUT = ROOT / "results" / "013_gibberish_and_cipher_controls"
GASKELL = Path(os.environ.get("VOYNICH_GIBBERISH", ROOT.parent / "ext" / "gaskell"))
GIB_ZIP = GASKELL / "data" / "gibberish_transcriptions.zip"
NAIBBE = M9.UNITS / "data" / "controls" / "naibbe"
NAIBBE_FILES = {"naibbe_caesar": "naibbe_caesar_ciphertext_lines.txt",
                "naibbe_pliny": "greshko_nathist_output_ciphertext.txt",
                "naibbe_plaintext": "greshko_nathist_pre_encryption_respaced_plaintext.txt"}
WRITER_SEED = 13
N_BLOCKS = 7
SOURCES = ("gibberish_pooled", "gibberish_native", "voynich_crop", "caesar_crop", "celsus_crop",
           "naibbe_caesar", "naibbe_pliny", "naibbe_plaintext")
CROP_CONTROLS = {"caesar_crop": "caesar_prose", "celsus_crop": "celsus_prose"}
N_SUB = 100
N_SUB_011 = 30
MIN_STEM_WRITER = 5


# --------------------------------------------------------------------------- #
# sources
# --------------------------------------------------------------------------- #

_CACHE: dict = {}


def writers() -> dict:
    if "writers" not in _CACHE:
        _CACHE["writers"] = C.read_gibberish(GIB_ZIP)
    return _CACHE["writers"]


def gib_divisions() -> list:
    w = writers()
    ids = list(w)
    order = [ids[i] for i in np.random.default_rng(WRITER_SEED).permutation(len(ids))]
    return [(k, [ln for page in w[k] for ln in page]) for k in order]


def zl() -> pd.DataFrame:
    return M9.real_table("ZL_split")[0]


def crop() -> pd.DataFrame:
    if "crop" not in _CACHE:
        _CACHE["crop"] = C.crop_template(zl(), C.n_tokens(gib_divisions()))
    return _CACHE["crop"]


def build(name: str) -> pd.DataFrame:
    """Token table of a 013 source (``attrs['control']`` describes the layout)."""
    if name == "gibberish_pooled":
        return C.stream_control(gib_divisions(), crop(), name, writer_as_hand=True)
    if name == "gibberish_native":
        pages = [(f"{w}.{i + 1}", w, w, page) for w, pp in writers().items() for i, page in enumerate(pp)]
        t = C.lines_table(pages, name=name)
        t.attrs["control"] = {"name": name, "mode": "native", "strategy": "own lines and pages",
                              "sections": list(writers()), "source_tokens": len(t)}
        return t
    if name == "voynich_crop":
        t = crop().copy()
        t.attrs["control"] = {"name": name, "mode": "real", "strategy": "evenly spaced pages",
                              "sections": sorted(t["illustration"].astype(str).unique()), "source_tokens": len(zl())}
        return t
    if name in CROP_CONTROLS:
        ctl = CROP_CONTROLS[name]
        _, _, loader, mode, keep = M9.CONTROLS[ctl]
        divs, paths = loader()
        t = C.build_control(divs, crop(), M9.K_SECTIONS, mode, name, keep=keep, min_k=M9.MIN_K)
        t.attrs["control"]["files"] = [M9.rel_source(p) for p in paths]
        return t
    if name in NAIBBE_FILES:
        path = NAIBBE / NAIBBE_FILES[name]
        toks = [t for ln in C.read_lines(path) for t in ln]
        tpl = zl()
        divs = C.blocks(toks, N_BLOCKS, limit=len(tpl))
        eva = name != "naibbe_plaintext"
        t = C.stream_control(divs, tpl, name, alphabet="Eva-" if eva else "latin",
                             glyph_preset=tpl.attrs.get("glyph_preset") if eva else None)
        t.attrs["control"]["source_tokens"] = len(toks)
        t.attrs["control"]["files"] = [M9.rel_source(path)]
        return t
    raise ValueError(name)


class Ctx(M9.Ctx):
    """009 context for a 013 source (affixes learned on the source itself, page ranks of ZL)."""

    def __init__(self, name: str):
        self.source = f"ctl:{name}"  # 'ctl:' -> 011 measures decomposition (a) only
        self.base_tag = M9.TEMPLATE_TAG
        if name == "gibberish_native":
            pooled = build("gibberish_pooled")
            d = decompose_table(pooled)
            self.real = build(name).reset_index(drop=True)
        else:
            self.real = build(name).reset_index(drop=True)
            d = decompose_table(self.real)
        pages = list(pd.unique(self.real["folio"]))
        rank0 = M9.real_table(M9.TEMPLATE_TAG)[1]
        self.rank = {p: rank0.get(p, i) for i, p in enumerate(pages)}
        eva = str(self.real.attrs.get("alphabet", "")).lower().startswith("eva")
        self.gpreset = "extended" if eva else None
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]
        self.enc_ns = SimpleNamespace(prefixes=self.prefixes, suffixes=self.suffixes, gpreset=self.gpreset, merges=[])
        self._tmpl = None
        self._m7 = None


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def _run(kind: str, source: str, seeds_arg: str | None) -> None:
    ctx = Ctx(source)
    n_sub = N_SUB_011 if kind == "p011" else N_SUB

    def one(seed: int) -> dict:
        df = ctx.real if seed < 0 else M9.page_subsample(ctx.real, seed)
        level = "point" if seed < 0 else "sub"
        if kind == "p009":
            m = M9.profile(df, ctx, seed, level)
        elif kind == "p010":
            m = M10.line_profile(df, ctx, seed, level)
        else:
            m = M11.measure(df, ctx, seed, level)
        return {"source": f"ctl:{source}", "kind": level, **m}

    seeds = M9.parse_seeds(seeds_arg, [-1] + list(range(n_sub)))
    params = {"source": source, "kind": kind, "prefixes": ctx.prefixes, "suffixes": ctx.suffixes,
              "glyph_units": ctx.gpreset, "control": ctx.real.attrs.get("control"), "writer_seed": WRITER_SEED}
    run_seeds(f"013_{kind}", seeds, one, OUT / "runs" / kind / source, params=params, progress=progress)


def cmd_sources(args) -> None:
    rows = []
    for name in SOURCES:
        t = build(name)
        info = t.attrs.get("control", {})
        v = vocabulary_summary(t).iloc[0]
        n_lines = t.groupby(["folio", "locus_id"], sort=False).ngroups
        files = info.get("files") or ([f"danielgaskell/voynich/data/{GIB_ZIP.name}"] if "gibberish" in name else
                                      ["data/ZL3b-n.txt"])
        rows.append({"source": name, "files": ";".join(files), "layout": info.get("mode"),
                     "strategy": info.get("strategy"), "source_tokens": info.get("source_tokens"), "tokens": len(t),
                     "lines": n_lines, "pages": t["folio"].nunique(), "sections": t["illustration"].nunique(),
                     "writers": len(writers()) if "gibberish" in name else np.nan,
                     "n_types": int(v["n_types"]), "hapax_share": round(float(v["hapax_share"]), 4),
                     "alphabet": t.attrs.get("alphabet"), "first_tokens": " ".join(t["token"].head(10))})
        print(name, len(t), flush=True)
    w = writers()
    per = pd.DataFrame([{"writer": k, "pages": len(p), "lines": sum(len(pg) for pg in p),
                         "tokens": sum(len(ln) for pg in p for ln in pg)} for k, p in w.items()])
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "sources_used.csv", index=False)
    per.to_csv(OUT / "writers_inventory.csv", index=False)
    print(per.describe().round(1).to_string())


# --------------------------------------------------------------------------- #
# per writer
# --------------------------------------------------------------------------- #


def writer_metrics(t: pd.DataFrame, prefixes, suffixes, preset) -> dict:
    """Short-text metrics of one writer (or one Voynich chunk of the same size)."""
    v = vocabulary_summary(t).iloc[0]
    gl = glyph_lists(t, preset)
    lens = np.asarray([len(g) for g in gl], dtype=float)
    d = decompose_table(t, prefixes=prefixes, suffixes=suffixes)
    sc = d["stem"].value_counts()
    stems = sc.index[sc >= MIN_STEM_WRITER]
    sub = d[d["stem"].isin(stems)]
    nsuf = sub.groupby("stem")["suffix"].nunique()
    hsg = np.nan
    if len(stems) >= 3:
        hsg = float(sum(len(g) / len(sub) * entropy(g["suffix"].value_counts().to_numpy())
                        for _, g in sub.groupby("stem")))
    last = t["is_last"].to_numpy(dtype=bool)
    return {"n_tokens": len(t), "n_lines": t.groupby(["folio", "locus_id"], sort=False).ngroups,
            "n_types": int(v["n_types"]), "hapax_share": float(v["hapax_share"]),
            "H1_glyph": entropy(pd.Series([x for g in gl for x in g]).value_counts().to_numpy()), "H2_glyph": conditional_entropy(gl, order=1),
            "mean_len": float(lens.mean()), "n_stems_5": int(len(stems)),
            "mean_suffixes_per_stem": float(nsuf.mean()) if len(stems) >= 3 else np.nan,
            "H_suffix_given_stem": hsg,
            "end_len_diff": float(lens[last].mean() - lens[~last].mean()) if last.any() and (~last).any() else np.nan}


def cmd_writers(args) -> None:
    ctx = Ctx("gibberish_native")
    t = ctx.real
    rows = []
    for w in writers():
        tw = t[t["hand"] == w]
        rows.append({"group": "gibberish", "unit": w, **writer_metrics(tw, ctx.prefixes, ctx.suffixes, None)})
    # Voynich chunks of the same sizes: consecutive whole lines of ZL split from evenly spaced starting lines
    vctx = M9.Ctx("ZL_split")
    real = vctx.real
    lid = real.groupby(["folio", "locus_id"], sort=False).ngroup().to_numpy()
    starts = np.flatnonzero(np.r_[True, lid[1:] != lid[:-1]])
    sizes = [r["n_tokens"] for r in rows]
    for i, n in enumerate(sizes):
        s0 = starts[int(i * len(starts) / len(sizes))]
        end = s0 + n
        cut = starts[starts >= end]
        stop = int(cut[0]) if len(cut) else len(real)
        chunk = real.iloc[s0:stop]
        rows.append({"group": "voynich_ZL_chunk", "unit": f"chunk{i + 1:02d}",
                     **writer_metrics(chunk, vctx.prefixes, vctx.suffixes, vctx.gpreset)})
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.round(5).to_csv(OUT / "writers.csv", index=False)
    print(out.drop(columns="unit").groupby("group").agg(["mean", "std"]).T.round(3).to_string())


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("sources")
    for name in ("profile", "lines", "variation"):
        p = sub.add_parser(name)
        p.add_argument("--source", required=True, choices=SOURCES)
        p.add_argument("--seeds", default=None)
    sub.add_parser("writers")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from summarize import summarize  # noqa: E402  (experiments/013_gibberish_and_cipher_controls/summarize.py)

        summarize()
        return
    if args.cmd == "sources":
        cmd_sources(args)
    elif args.cmd == "writers":
        cmd_writers(args)
    else:
        _run({"profile": "p009", "lines": "p010", "variation": "p011"}[args.cmd], args.source, args.seeds)


if __name__ == "__main__":
    main()
