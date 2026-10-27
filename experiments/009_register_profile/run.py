"""Experiment 009: register profile - where the Voynich text lies between narrative, catalogue and generator.

One profile of ~37 metrics in 7 groups (vocabulary, token order, position vs section, paradigm,
locality, sequence, line) is measured with the same code on

* the Voynich paragraph text (``ZL_split``, ``ZL_join``, ``GC_split``),
* real-language controls from ``lrozanova/voynich-units`` (and optionally a second Pliny edition from
  ``alexanderboxer/voynich-attack``) set into the line/page layout of ZL split (``voytools.controls``),
* null models: ts2020 calB (004 config 21), structured template (005 full model), auto-copy 86 (003),
  within-line shuffle.

Usage (from the repository root)::

    python experiments/009_register_profile/run.py controls            # controls_used.csv
    python experiments/009_register_profile/run.py profile --source ZL_split            # point + 100 subsamples
    python experiments/009_register_profile/run.py profile --source ctl:celsus_prose
    python experiments/009_register_profile/run.py profile --source ts2020:ZL_split     # 30 seeds
    python experiments/009_register_profile/run.py summarize

Seeds: real text and controls ``-1`` = full text (point estimate, more shuffles), ``0..99`` = subsamples
of 80 % of the pages without replacement; null sources ``0..29``.  Every command is resumable
(``run_seeds``); outputs go to ``results/009_register_profile/runs/<source>/``.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import controls as C  # noqa: E402
from voytools import filter_tokens, load_runs, load_tables, run_seeds  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.nulls import shuffle_within_lines  # noqa: E402
from voytools.selfcitation import selfcitation_like  # noqa: E402
from voytools.stats import adjacent_token_mi, edge_glyph_mi, jsd, page_similarity, page_term_counts, \
    vocabulary_summary, zipf_fit  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M4 = _load("exp004", "004_selfcitation_port")
M5 = _load("exp005", "005_structured_template_null")
M6 = _load("exp006", "006_topic_at_distance")
M7 = _load("exp007", "007_position_vs_content")

DATA = ROOT / "data"
OUT = ROOT / "results" / "009_register_profile"
UNITS = Path(os.environ.get("VOYNICH_UNITS", ROOT.parent / "lrozanova" / "voynich-units"))
BOXER = Path(os.environ.get("VOYNICH_ATTACK", ROOT.parent / "ext" / "boxer"))
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
REAL_TAGS = ("ZL_split", "ZL_join", "GC_split")
TEMPLATE_TAG = "ZL_split"   # layout of every control
K_SECTIONS = 7              # ZL3b paragraph text has 7 illustration codes
MIN_K = 4
N_SUB = 100
SUB_FRAC = 0.8
N_NULL = 30
TS2020_CONFIG = 21          # 004 best_configs.csv row ts2020_calB
AC_SEED0 = M7.AC_SEED0      # 003 phase-B seeds of config 86
TMPL_SEED0 = 9000
NULLS = {"ts2020": ("ZL_split", "ZL_join"), "tmpl": ("ZL_split", "GC_split"),
         "ac86": ("ZL_split", "GC_split"), "line": ("ZL_split", "GC_split")}
SHUF = {"point": {"adj": 50, "edge": 50, "bias": 50, "mi": 50, "line": 50},
        "sub": {"adj": 50, "edge": 50, "bias": 10, "mi": 10, "line": 20}}
TOP_CLASSES = 5
MAIN_K = 30                 # 006: function stems excluded from page similarity
DMAX = 10
POS_TERMS = (("I_pos", "pos", None), ("I_sec", "sec", None), ("I_sec_g_pos", "sec", "pos"))

GROUPS = {
    "vocab": ["hapax_share", "n_types", "top50_coverage", "zipf_slope"],
    "order": ["token_mi_excess_share", "edge_mi_excess"],
    "position": [f"{y}__{t}" for y in M7.TARGETS for t, _, _ in POS_TERMS],
    "paradigm": ["mean_suffixes_per_stem", "frac_stems_multi_suffix", "H_suffix_given_stem", "mi_prefix_suffix",
                 "cmi_excess", "stem_dispersion", "jsd_stem_sections"],
    "locality": ["loc_log_dist", "loc_same_section", "loc_contrast_adj"],
    "sequence": ["mi_stem_d1", "mi_stem_area", "mi_affix_d1", "mi_affix_area"],
    "line": ["first_top5_excess", "last_top5_excess", "jsd_first_mid", "jsd_last_mid", "I_first_rest"],
}
METRICS = [m for g in GROUPS.values() for m in g]


# --------------------------------------------------------------------------- #
# controls
# --------------------------------------------------------------------------- #

def _bundle(name: str) -> Path:
    return UNITS / "voynich_decipherment_repro_bundle" / "decipherment_attack_v5" / "latin_controls" / name


def _vu(rel: str) -> Path:
    return UNITS / "data" / "voynich-units" / rel


def _tess(pattern: str) -> tuple[list, list[str]]:
    paths = sorted(_vu("").glob(pattern), key=lambda p: (len(p.name), p.name))
    return [d for p in paths for d in C.parse_tess(C.read_text(p))], paths


def _caesar():
    paths = [_bundle("caesar_la.txt"), _bundle("caesar_la2.txt")]
    a = C.parse_headed(C.gutenberg_body(C.read_text(paths[0])), r"C\. IULI CAESARIS DE BELLO GALLICO COMMENTARIUS (\w+)",
                       "{1}")
    b = C.parse_headed(C.gutenberg_body(C.read_text(paths[1])), r"Liber ([IVX]+)", "{1}")
    return a + b, paths


def _linnaeus():
    paths = [_bundle("species1.txt"), _bundle("species2.txt")]
    return [d for p in paths for d in C.parse_linnaeus(C.read_text(p))], paths


def _culpeper(item: str):
    path = _vu("herbal/culpeper_en.txt")
    return C.parse_culpeper(C.read_text(path), item=item), [path]


def _italian():
    path = _vu("italian/27608.txt")
    return C.blocks(C.words(C.gutenberg_body(C.read_text(path))), 14), [path]


def _pliny_perseus():
    path = BOXER / "corpora" / "latin" / "legacy" / "pliny" / "pliny_lat0.csv"
    df = pd.read_csv(path)
    df.columns = ["ref", "text"]
    pairs = [(str(r).split(".")[0], C.words(str(t))) for r, t in zip(df["ref"], df["text"])]
    divs = []
    for label, item in pairs:
        if divs and divs[-1][0] == label:
            divs[-1][1].append(item)
        else:
            divs.append((label, [item]))
    return divs, [path]


PLINY_BOOKS = [str(b) for b in range(16, 31)]  # the books of the voynich-units copy
# name -> (genre, language, loader, layout mode, keep)
CONTROLS = {
    "caesar_prose": ("narrative", "la", _caesar, "prose", None),
    "celsus_prose": ("medical", "la", lambda: _tess("latin/celsus*.tess"), "prose", None),
    "pliny_prose": ("botanical", "la", lambda: _tess("herbal/pliny*.tess"), "prose", None),
    "plinyP_prose": ("botanical", "la", _pliny_perseus, "prose", PLINY_BOOKS),
    "culpeper_prose": ("herbal", "en", lambda: _culpeper("paragraph"), "prose", None),
    "culpeper_items": ("herbal", "en", lambda: _culpeper("sentence"), "items", None),
    "italian_prose": ("prose", "it", _italian, "prose", [f"b{i}" for i in range(2, 15)]),  # b1: front matter
    "linnaeus_prose": ("catalogue", "la", _linnaeus, "prose", None),
    "linnaeus_items": ("catalogue", "la", _linnaeus, "items", None),
}


def control_available(name: str) -> bool:
    return name != "plinyP_prose" or (BOXER / "corpora" / "latin" / "legacy" / "pliny" / "pliny_lat0.csv").exists()


def rel_source(p: Path) -> str:
    """Path of a source file inside its repository, prefixed with the repository name."""
    for name, base in (("voynich-units", UNITS), ("voynich-attack", BOXER)):
        try:
            return f"{name}/{Path(p).resolve().relative_to(base.resolve()).as_posix()}"
        except ValueError:
            pass
    return str(p)


_REAL: dict = {}


def real_table(tag: str) -> pd.DataFrame:
    if tag not in _REAL:
        key, policy = tag.split("_")
        t = load_tables(DATA / FILES[key], space_policy=policy)
        _REAL[tag] = (filter_tokens(t.tokens, generic_type="P", clean=True).reset_index(drop=True),
                      {p: i for i, p in enumerate(pd.unique(t.tokens["folio"].to_numpy()))})
    return _REAL[tag]


def control_table(name: str) -> pd.DataFrame:
    _, _, loader, mode, keep = CONTROLS[name]
    divs, paths = loader()
    tpl = real_table(TEMPLATE_TAG)[0]
    out = C.build_control(divs, tpl, K_SECTIONS, mode, name, keep=keep, min_k=MIN_K)
    out.attrs["control"]["files"] = [rel_source(p) for p in paths]
    return out


# --------------------------------------------------------------------------- #
# corpus context
# --------------------------------------------------------------------------- #


class Ctx:
    """A full corpus (real text or control) with its affix lists, glyph units and page ranks."""

    def __init__(self, source: str):
        self.source = source
        if source.startswith("ctl:"):
            self.real = control_table(source[4:]).reset_index(drop=True)
            _, self.rank = real_table(TEMPLATE_TAG)
            self.base_tag = TEMPLATE_TAG
        else:
            self.base_tag = source.split(":")[-1]
            self.real, self.rank = real_table(self.base_tag)
        eva = bool(str(self.real.attrs.get("alphabet", "")).lower().startswith("eva"))
        self.gpreset = "extended" if eva else None
        d = decompose_table(self.real)
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]
        self.enc_ns = SimpleNamespace(prefixes=self.prefixes, suffixes=self.suffixes, gpreset=self.gpreset, merges=[])
        self._tmpl = None
        self._m7 = None

    def corpus(self, kind: str, seed: int) -> pd.DataFrame:
        if kind == "ts2020":
            return selfcitation_like(self.real, seed=seed, **M4.config_params(TS2020_CONFIG))
        if kind == "tmpl":
            if self._tmpl is None:
                _, cfg, _ = M5.ablation_base(self.base_tag)
                self._tmpl = M5.template(self.base_tag).with_params(**{**cfg, "suffix_by": None})
            return self._tmpl.sample(self.real, seed=TMPL_SEED0 + seed, pool_seed=seed)
        if kind == "ac86":
            if self._m7 is None:
                key, policy = self.base_tag.split("_")
                self._m7 = M7.Context(FILES[key], policy)
            return self._m7.autocopy(AC_SEED0 + seed)
        if kind == "line":
            return shuffle_within_lines(self.real, seed=seed)
        raise ValueError(kind)


def page_subsample(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    pages = pd.unique(df["folio"].to_numpy())
    rng = np.random.default_rng(seed)
    keep = rng.choice(pages, size=int(round(SUB_FRAC * len(pages))), replace=False)
    return df[df["folio"].isin(keep)].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #


def m_vocab(df: pd.DataFrame) -> dict:
    v = vocabulary_summary(df).iloc[0]
    counts = df["token"].value_counts().to_numpy()
    return {"hapax_share": float(v["hapax_share"]), "n_types": float(v["n_types"]),
            "top50_coverage": float(v["top50_coverage"]), "zipf_slope": float(zipf_fit(counts)["slope"])}


def m_order(df: pd.DataFrame, ctx: Ctx, k: dict, seed: int) -> dict:
    a = adjacent_token_mi(df, cap=2000, shuffles=k["adj"], seed=seed).iloc[0]
    e = edge_glyph_mi(df, shuffles=k["edge"], seed=seed, preset=ctx.gpreset).iloc[0]
    return {"token_mi_excess_share": float(a["excess_share"]), "edge_mi_excess": float(e["excess_bits"])}


def m_position(enc: dict, k: dict, seed: int) -> dict:
    """007 decomposition restricted to position and section: raw minus shuffle bias (Y permuted within Z)."""
    rng = np.random.default_rng(seed)
    X = {kk: v[0] for kk, v in enc["X"].items()}
    zero = np.zeros(enc["n"], dtype=np.int64)
    sh = {None: M7.Shuffler(zero), "pos": M7.Shuffler(X["pos"])}
    out = {}
    for yname in M7.TARGETS:
        y = enc["Y"][yname][0]
        for name, xn, zn in POS_TERMS:
            z = None if zn is None else X[zn]
            raw = M7.cmi(y, X[xn], z)
            bias = float(np.mean([M7.cmi(y[sh[zn].src(rng)], X[xn], z) for _ in range(k["bias"])]))
            out[f"{yname}__{name}"] = raw - bias
    return out


def m_paradigm(dfd: pd.DataFrame, seed: int) -> dict:
    m = M5.metrics_decomposed(dfd, seed)
    return {c: float(m[c]) for c in GROUPS["paradigm"]}


def m_locality(dfd: pd.DataFrame, rank: dict) -> dict:
    """006 model ``spec+size`` (cosine, K = 30 function stems excluded) and size-adjusted contrast."""
    size = dfd.groupby("folio", sort=False).size()
    pages = [p for p in size.index if size[p] >= M6.MIN_TOKENS]
    meta = dfd.groupby("folio", sort=False)[["illustration", "language", "hand"]].first().loc[pages]
    vc = dfd.loc[dfd["stem"] != "", "stem"].value_counts()
    counts = page_term_counts(dfd, "stem", exclude=vc.index[:MAIN_K].tolist()).reindex(pages, fill_value=0)
    S = page_similarity(counts.to_numpy(), "cosine")
    ns = SimpleNamespace(rank=np.asarray([rank[p] for p in pages]), sec=meta["illustration"].astype(str).to_numpy(),
                         lang=meta["language"].astype(str).to_numpy(), hand=meta["hand"].astype(str).to_numpy())
    iu, ju = np.triu_indices(len(pages), 1)
    sz = size.loc[pages].to_numpy(dtype=float)
    r = M6.pair_metrics(S[iu, ju], M6.pair_design(ns, iu, ju), sz[iu], sz[ju])
    return {"loc_log_dist": r["s_log_dist"], "loc_same_section": r["s_same_section"],
            "loc_contrast_adj": r["contrast_adj"]}


def m_sequence(enc: dict, k: dict, seed: int) -> dict:
    """007 MI(d) of stems and affix frames with the ``pc`` correction (position class in the paragraph)."""
    para = enc["para"]
    same = {d: para[:-d] == para[d:] for d in range(1, DMAX + 1)}
    sub = dict(enc, U={u: enc["U"][u] for u in ("stem", "affix")})
    units = M7.UNITS
    M7.UNITS = ("stem", "affix")
    try:
        raw = M7.mi_curve(sub, same)
        g = M7.correction_groups(enc)["pc"]
        sh = M7.Shuffler(np.asarray(g, dtype=np.int64))
        rng = np.random.default_rng(seed)
        acc = {key: 0.0 for key in raw}
        for _ in range(k["mi"]):
            for key, v in M7.mi_curve(sub, same, sh.src(rng)).items():
                acc[key] += v / k["mi"]
    finally:
        M7.UNITS = units
    out = {}
    for u in ("stem", "affix"):
        corr = {d: raw[(u, d)] - acc[(u, d)] for d in range(1, DMAX + 1)}
        out[f"mi_{u}_d1"] = corr[1]
        out[f"mi_{u}_area"] = float(sum(corr[d] for d in range(2, DMAX + 1)))
    return out


def _top_share(c: np.ndarray, kk: int) -> float:
    n = np.bincount(c)
    return float(np.sort(n)[::-1][:kk].sum() / n.sum()) if n.sum() else float("nan")


def m_line(df: pd.DataFrame, enc: dict, k: dict, seed: int) -> dict:
    """Line conventions on affix-frame classes (prefix|suffix, >= 20 occurrences, else <other>)."""
    cls, kc = enc["U"]["affix"][0]
    first = df["is_first"].to_numpy(dtype=bool)
    last = df["is_last"].to_numpy(dtype=bool)
    mid = ~(first | last)
    out = {"first_top5_excess": _top_share(cls[first], TOP_CLASSES) - _top_share(cls[mid], TOP_CLASSES),
           "last_top5_excess": _top_share(cls[last], TOP_CLASSES) - _top_share(cls[mid], TOP_CLASSES),
           "jsd_first_mid": jsd(np.bincount(cls[first], minlength=kc), np.bincount(cls[mid], minlength=kc)),
           "jsd_last_mid": jsd(np.bincount(cls[last], minlength=kc), np.bincount(cls[mid], minlength=kc))}
    # I(class of the first token; class of each other token of the line), lines' first classes permuted
    line = enc["line"]
    first_idx = np.flatnonzero(first)
    lid_first = line[first_idx]
    fmap = np.full(int(line.max()) + 1, -1, dtype=np.int64)
    rest = ~first
    rng = np.random.default_rng(seed)

    def mi(fc: np.ndarray) -> float:
        fmap[:] = -1
        fmap[lid_first] = fc
        f = fmap[line[rest]]
        ok = f >= 0
        return M7.mi_pairs(f[ok], cls[rest][ok], kc, kc)

    fc0 = cls[first_idx]
    raw = mi(fc0)
    sims = [mi(fc0[rng.permutation(len(fc0))]) for _ in range(k["line"])]
    out["I_first_rest"] = raw - float(np.mean(sims))
    return out


def profile(df: pd.DataFrame, ctx: Ctx, seed: int, level: str) -> dict:
    """All metrics of one corpus (``level`` 'point' or 'sub' sets the number of shuffles)."""
    k = SHUF[level]
    s = max(seed, 0)
    df = df.reset_index(drop=True)
    dfd = decompose_table(df, prefixes=ctx.prefixes, suffixes=ctx.suffixes)
    enc = M7.Context.encode(ctx.enc_ns, df)
    out = {"n_tokens": len(df), "n_pages": int(df["folio"].nunique()), "n_sections": int(df["illustration"].nunique())}
    out.update(m_vocab(df))
    out.update(m_order(df, ctx, k, s))
    out.update(m_position(enc, k, s))
    out.update(m_paradigm(dfd, s))
    out.update(m_locality(dfd, ctx.rank))
    out.update(m_sequence(enc, k, s))
    out.update(m_line(df, enc, k, s))
    return out


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def safe(source: str) -> str:
    return source.replace(":", "__")


def parse_seeds(spec: str | None, default: list[int]) -> list[int]:
    if not spec:
        return default
    if ":" in spec:
        a, b = spec.split(":")
        return list(range(int(a), int(b)))
    return [int(x) for x in spec.split(",") if x.strip()]


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_profile(args) -> None:
    source = args.source
    kind = source.split(":")[0] if ":" in source and not source.startswith("ctl:") else "real"
    ctx = Ctx(source if kind == "real" else f"real:{source.split(':')[1]}")
    if kind == "real":
        def one(seed: int) -> dict:
            df = ctx.real if seed < 0 else page_subsample(ctx.real, seed)
            return {"source": source, "kind": "point" if seed < 0 else "sub",
                    **profile(df, ctx, seed, "point" if seed < 0 else "sub")}
        default = [-1] + list(range(N_SUB))
    else:
        def one(seed: int) -> dict:
            return {"source": source, "kind": "null", **profile(ctx.corpus(kind, seed), ctx, seed, "sub")}
        default = list(range(N_NULL))
    seeds = parse_seeds(args.seeds, default)
    params = {"source": source, "prefixes": ctx.prefixes, "suffixes": ctx.suffixes, "glyph_units": ctx.gpreset,
              "shuffles": SHUF, "sub_frac": SUB_FRAC, "n_sub": N_SUB, "n_null": N_NULL,
              "ts2020_config": TS2020_CONFIG, "ac_seed0": AC_SEED0, "tmpl_seed0": TMPL_SEED0,
              "control": ctx.real.attrs.get("control")}
    run_seeds("009_register_profile", seeds, one, OUT / "runs" / safe(source), params=params, progress=progress)


def cmd_controls(args) -> None:
    rows = []
    for name, (genre, lang, loader, mode, keep) in CONTROLS.items():
        if not control_available(name):
            print(f"{name}: source not available", flush=True)
            continue
        t = control_table(name)
        info = t.attrs["control"]
        sec = t.groupby("illustration", sort=False)["folio"].nunique()
        v = vocabulary_summary(t).iloc[0]
        rows.append({"control": name, "genre": genre, "language": lang, "layout": mode, "strategy": info["strategy"],
                     "files": ";".join(info["files"]), "source_tokens": info["source_tokens"], "tokens": len(t),
                     "n_sections": len(info["sections"]), "sections": " | ".join(info["sections"]),
                     "pages_per_section": " | ".join(f"{s}:{n}" for s, n in sec.items()),
                     "n_types": int(v["n_types"]), "hapax_share": round(float(v["hapax_share"]), 4),
                     "first_tokens": " ".join(t["token"].head(12))})
        print(name, info["strategy"], len(info["sections"]), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "controls_used.csv", index=False)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profile")
    p.add_argument("--source", required=True)
    p.add_argument("--seeds", default=None)
    sub.add_parser("controls")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        from summarize import summarize  # noqa: E402  (experiments/009_register_profile/summarize.py)

        summarize()
        return
    {"profile": cmd_profile, "controls": cmd_controls}[args.cmd](args)


if __name__ == "__main__":
    main()
