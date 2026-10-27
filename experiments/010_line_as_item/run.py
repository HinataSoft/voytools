"""Experiment 010: the line as an item - line independence, line templates, vertical copying.

Six line-level metric groups, measured with the same code on the real text (ZL/GC x split/join), the
null models (ts2020 calB, template 005, auto-copy 86, within-line shuffle, position-class shuffle,
permutation of whole lines of equal length) and the real-language controls of experiment 009
(laid out in the ZL split layout with ``voytools.controls``; prose broken mechanically and
catalogue/herbal laid out item per line):

1. boundary - corrected MI at d = 1 inside the line / across a line / across a paragraph boundary
   (``stats.line_boundary_mi``) for edge glyphs, affix classes and stems;
2. header - I(class of the first token; class of the other tokens / 2nd / 3rd / last token),
   null: first tokens permuted among lines of equal length;
3. templates - entropies of the class by position, line signatures (lines of 4..8 tokens),
   zlib compression of the class sequence vs within-line and same-length-column shuffles;
4. same length - Jaccard similarity of class sets of line pairs, same vs different length, minus the
   expectation after a global token shuffle that keeps line lengths;
5. vertical copying - ``stats.vertical_similarity`` (edit distance <= 1 to the token above);
6. line end - ``-m``/``-g`` endings (Eva only), length of the last token, end-glyph enrichment.

Usage (from the repository root)::

    python experiments/010_line_as_item/run.py profile --source ZL_split        # point + 100 subsamples
    python experiments/010_line_as_item/run.py profile --source ctl:linnaeus_items
    python experiments/010_line_as_item/run.py profile --source ts2020:ZL_split  # 30 seeds
    python experiments/010_line_as_item/run.py summarize

Seeds as in 009: ``-1`` = full text, ``0..99`` page subsamples (80 %), nulls ``0..29``.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import run_seeds  # noqa: E402
from voytools._common import glyph_lists, line_ids  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.nulls import shuffle_within_position_class  # noqa: E402
from voytools.stats import line_boundary_mi, vertical_similarity  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R9 = _load("exp009", "009_register_profile")  # controls, contexts and null corpora of 009
M7 = R9.M7

OUT = ROOT / "results" / "010_line_as_item"
REAL_TAGS = ("ZL_split", "ZL_join", "GC_split", "GC_join")
NULLS = {"ts2020": ("ZL_split", "ZL_join"), "tmpl": ("ZL_split", "GC_split"), "ac86": ("ZL_split", "GC_split"),
         "line": ("ZL_split", "GC_split"), "posclass": ("ZL_split", "GC_split"), "samelen": ("ZL_split", "GC_split")}
CONTROLS = ("caesar_prose", "pliny_prose", "celsus_prose", "italian_prose", "culpeper_prose", "culpeper_items",
            "linnaeus_prose", "linnaeus_items")
N_SUB, N_NULL = 100, 30
N_CLASSES = 20       # affix classes: 20 most frequent (prefix, suffix) pairs + other
TOP_STEMS = 200
SIG_LENGTHS = range(4, 9)
N_JAC_PAIRS = 20000
K = {"point": {"bnd": 50, "head": 200, "comp": 5, "jac": 5, "vert": 200},
     "sub": {"bnd": 10, "head": 20, "comp": 2, "jac": 2, "vert": 20}}
END_GLYPHS = ("m", "g")


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


def permute_equal_length_lines(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Null: whole lines (token sequences) permuted among the lines of equal length; layout kept."""
    rng = np.random.default_rng(seed)
    toks = df["token"].to_numpy(dtype=object)
    lid = line_ids(df)
    starts = np.flatnonzero(np.r_[True, lid[1:] != lid[:-1]])
    lens = np.diff(np.r_[starts, len(df)])
    src = np.arange(len(starts))
    for n in np.unique(lens):
        idx = np.flatnonzero(lens == n)
        src[idx] = idx[rng.permutation(len(idx))]
    order = np.concatenate([np.arange(starts[s], starts[s] + lens[s]) for s in src])
    out = df.copy()
    for c in ("token", "token_raw", "len_chars", "len_glyphs"):
        if c in out.columns:
            out[c] = out[c].to_numpy()[order]
    out.attrs.update(df.attrs)
    return out


def corpus(ctx, kind: str, seed: int) -> pd.DataFrame:
    if kind == "posclass":
        return shuffle_within_position_class(ctx.real, seed=seed)
    if kind == "samelen":
        return permute_equal_length_lines(ctx.real, seed)
    return ctx.corpus(kind, seed)


# --------------------------------------------------------------------------- #
# encoding
# --------------------------------------------------------------------------- #


def codes(values) -> np.ndarray:
    return pd.factorize(pd.Series(list(values), dtype=object))[0].astype(np.int64)


def encode(df: pd.DataFrame, ctx) -> dict:
    """Affix class (20 most frequent prefix|suffix + other), capped stem, edge glyphs, line structure."""
    dfd = decompose_table(df, prefixes=ctx.prefixes, suffixes=ctx.suffixes)
    frame = (dfd["prefix"].astype(str) + "|" + dfd["suffix"].astype(str)).to_numpy(dtype=object)
    vc = pd.Series(frame).value_counts()
    top = set(vc.index[:N_CLASSES])
    cls = np.asarray([f if f in top else "<other>" for f in frame], dtype=object)
    stem = dfd["stem"].to_numpy(dtype=object)
    sc = pd.Series(stem).value_counts()
    tops = set(sc.index[:TOP_STEMS])
    gl = glyph_lists(df, ctx.gpreset)
    lid = line_ids(df)
    starts = np.flatnonzero(np.r_[True, lid[1:] != lid[:-1]])
    lens = np.diff(np.r_[starts, len(df)])
    return {
        "cls_label": cls, "cls": codes(cls),
        "stem": np.asarray([s if s in tops else "<other>" for s in stem], dtype=object),
        "first_glyph": [g[0] if g else "" for g in gl], "last_glyph": [g[-1] if g else "" for g in gl],
        "len": np.asarray([len(g) for g in gl]),
        "starts": starts, "lens": lens, "line": lid,
        "tokpos": np.arange(len(df)) - np.repeat(starts, lens),
    }


def H(c: np.ndarray) -> float:
    if len(c) == 0:
        return float("nan")
    p = np.bincount(c) / len(c)
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def joint(*arrays: np.ndarray) -> np.ndarray:
    key = arrays[0].astype(np.int64)
    for a in arrays[1:]:
        key = key * (int(a.max()) + 1 if len(a) else 1) + a
    return pd.factorize(key)[0].astype(np.int64)


def mi(x: np.ndarray, y: np.ndarray) -> float:
    return H(x) + H(y) - H(joint(x, y))


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #


def m_boundary(df: pd.DataFrame, e: dict, k: dict, seed: int) -> dict:
    out = {}
    units = {"edge": (e["last_glyph"], e["first_glyph"]), "affix": (e["cls_label"], e["cls_label"]),
             "stem": (e["stem"], e["stem"])}
    for u, (left, right) in units.items():
        r = line_boundary_mi(df, left=left, right=right, shuffles=k["bnd"], seed=seed).set_index("kind")
        for kind, name in (("inside", "in"), ("line", "line"), ("paragraph", "par")):
            out[f"bnd_{u}_{name}"] = float(r.loc[kind, "excess_bits"])
        out[f"bnd_{u}_ratio_line"] = float(r.loc["line", "ratio_to_inside"])
        out[f"bnd_{u}_ratio_par"] = float(r.loc["paragraph", "ratio_to_inside"])
        if u == "edge":
            out["n_pairs_line"] = int(r.loc["line", "n_pairs"])
            out["n_pairs_par"] = int(r.loc["paragraph", "n_pairs"])
    return out


def m_header(e: dict, k: dict, seed: int) -> dict:
    """I(first class; X) - mean of the same with first classes permuted among lines of equal length."""
    cls, starts, lens = e["cls"], e["starts"], e["lens"]
    ok = lens >= 2
    st, ln = starts[ok], lens[ok]
    first = cls[st]
    rest_line = np.repeat(np.arange(len(st)), ln - 1)
    rest_tok = np.concatenate([np.arange(s + 1, s + n) for s, n in zip(st, ln)])
    targets = {"rest": (rest_line, cls[rest_tok])}
    for name, off in (("second", 1), ("third", 2)):
        m = ln > off
        targets[name] = (np.flatnonzero(m), cls[st[m] + off])
    m3 = ln >= 3
    targets["last"] = (np.flatnonzero(m3), cls[st[m3] + ln[m3] - 1])
    groups = {n: np.flatnonzero(ln == n) for n in np.unique(ln)}
    rng = np.random.default_rng(seed)
    obs = {t: mi(first[li], y) for t, (li, y) in targets.items()}
    sims = {t: [] for t in targets}
    for _ in range(k["head"]):
        perm = np.arange(len(st))
        for idx in groups.values():
            perm[idx] = idx[rng.permutation(len(idx))]
        f = first[perm]
        for t, (li, y) in targets.items():
            sims[t].append(mi(f[li], y))
    out = {}
    for t in targets:
        s = np.asarray(sims[t])
        out[f"head_{t}"] = obs[t] - float(s.mean())
        out[f"head_{t}_p"] = float((1 + (s >= obs[t]).sum()) / (1 + len(s)))
    return out


def _signatures(cls: np.ndarray, starts: np.ndarray, lens: np.ndarray) -> dict[int, pd.Series]:
    out = {}
    for n in SIG_LENGTHS:
        st = starts[lens == n]
        if len(st):
            sig = pd.Series([" ".join(map(str, cls[s: s + n])) for s in st])
            out[n] = sig.value_counts()
    return out


def _zratio(cls: np.ndarray, starts: np.ndarray, lens: np.ndarray, order: np.ndarray | None = None) -> float:
    c = cls if order is None else cls[order]
    sep = int(c.max()) + 1 if len(c) else 0
    seq = np.insert(c, starts[1:], sep).astype(np.uint8) if sep < 255 else None
    raw = seq.tobytes()
    return len(zlib.compress(raw, 9)) / len(raw)


def _column_shuffle(tokpos: np.ndarray, lens_tok: np.ndarray, rng) -> np.ndarray:
    """Order that permutes tokens with the same (line length, position) among lines."""
    key = lens_tok * 1000 + tokpos
    order = np.arange(len(key))
    for g in np.unique(key):
        idx = np.flatnonzero(key == g)
        order[idx] = idx[rng.permutation(len(idx))]
    return order


def _within_line_shuffle(line: np.ndarray, rng) -> np.ndarray:
    return np.lexsort((rng.random(len(line)), line))


def m_templates(e: dict, k: dict, seed: int) -> dict:
    cls, starts, lens, tokpos = e["cls"], e["starts"], e["lens"], e["tokpos"]
    lens_tok = np.repeat(lens, lens)
    sel = (lens_tok >= SIG_LENGTHS.start) & (lens_tok < SIG_LENGTHS.stop)
    pos_key = joint(lens_tok[sel], tokpos[sel])
    inside = e["line"][:-1] == e["line"][1:]
    out = {"tpl_H_class": H(cls),
           "tpl_H_given_pos": H(joint(pos_key, cls[sel])) - H(pos_key),
           "tpl_H_given_prev": H(joint(cls[:-1][inside], cls[1:][inside])) - H(cls[:-1][inside])}
    out["tpl_gain_pos"] = H(cls[sel]) - out["tpl_H_given_pos"]
    out["tpl_gain_prev"] = H(cls[1:][inside]) - out["tpl_H_given_prev"]
    sigs = _signatures(cls, starts, lens)
    n_lines = sum(int(s.sum()) for s in sigs.values())
    out["tpl_sig_per1000"] = 1000 * sum(len(s) for s in sigs.values()) / n_lines
    out["tpl_top20_cov"] = sum(int(s.iloc[:20].sum()) for s in sigs.values()) / n_lines
    rng = np.random.default_rng(seed)
    z = _zratio(cls, starts, lens)
    zl = np.mean([_zratio(cls, starts, lens, _within_line_shuffle(e["line"], rng)) for _ in range(k["comp"])])
    zc = np.mean([_zratio(cls, starts, lens, _column_shuffle(tokpos, lens_tok, rng)) for _ in range(k["comp"])])
    # the same signature statistics after the column shuffle (expectation from position-wise marginals)
    o = _column_shuffle(tokpos, lens_tok, rng)
    sigs0 = _signatures(cls[o], starts, lens)
    out["tpl_sig_per1000_colsh"] = 1000 * sum(len(s) for s in sigs0.values()) / n_lines
    out["tpl_top20_cov_colsh"] = sum(int(s.iloc[:20].sum()) for s in sigs0.values()) / n_lines
    out.update({"tpl_zratio": z, "tpl_zratio_linesh": float(zl), "tpl_zratio_colsh": float(zc),
                "tpl_zgain_linesh": float(zl - z), "tpl_zgain_colsh": float(zc - z)})
    return out


def _jaccard_pairs(cls: np.ndarray, starts: np.ndarray, lens: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    sets = [frozenset(cls[s: s + n].tolist()) for s, n in zip(starts, lens)]
    return np.asarray([len(sets[i] & sets[j]) / len(sets[i] | sets[j]) for i, j in zip(a, b)])


def m_same_length(e: dict, k: dict, seed: int) -> dict:
    cls, starts, lens = e["cls"], e["starts"], e["lens"]
    rng = np.random.default_rng(10_000 + seed)
    n = len(starts)
    a, b = rng.integers(n, size=N_JAC_PAIRS), rng.integers(n, size=N_JAC_PAIRS)
    keep = a != b
    a, b = a[keep], b[keep]
    same = lens[a] == lens[b]
    j = _jaccard_pairs(cls, starts, lens, a, b)
    ex_same, ex_diff = [], []
    for _ in range(k["jac"]):
        j0 = _jaccard_pairs(cls[rng.permutation(len(cls))], starts, lens, a, b)
        ex_same.append(j0[same].mean())
        ex_diff.append(j0[~same].mean())
    out = {"jac_same": float(j[same].mean()), "jac_diff": float(j[~same].mean()),
           "jac_same_exp": float(np.mean(ex_same)), "jac_diff_exp": float(np.mean(ex_diff))}
    out["jac_same_excess"] = out["jac_same"] - out["jac_same_exp"]
    out["jac_diff_excess"] = out["jac_diff"] - out["jac_diff_exp"]
    out["jac_samelen_effect"] = out["jac_same_excess"] - out["jac_diff_excess"]
    return out


def m_vertical(df: pd.DataFrame, k: dict, seed: int) -> dict:
    v = vertical_similarity(df, k=1, shuffles=k["vert"], seed=seed)
    return {f"vert_{c}": float(v[c]) for c in ("above", "above_left", "above_right", "else_above", "expected_above",
                                               "left_same", "excess_above", "ratio_above", "p_perm")}


def m_line_end(df: pd.DataFrame, e: dict, ctx) -> dict:
    last = df["is_last"].to_numpy(dtype=bool) & (np.repeat(e["lens"], e["lens"]) >= 2)
    other = ~df["is_last"].to_numpy(dtype=bool)
    lg = np.asarray(e["last_glyph"], dtype=object)
    out = {"end_len_last": float(e["len"][last].mean()), "end_len_other": float(e["len"][other].mean())}
    out["end_len_diff"] = out["end_len_last"] - out["end_len_other"]
    if ctx.gpreset is not None:  # Eva
        mg = np.isin(lg, END_GLYPHS)
        out["end_mg_last"] = float(mg[last].mean())
        out["end_mg_other"] = float(mg[other].mean())
    else:
        out["end_mg_last"] = out["end_mg_other"] = float("nan")
    vc = pd.Series(lg[other]).value_counts()
    cand = vc.index[vc >= 50]
    pl = pd.Series(lg[last]).value_counts(normalize=True)
    po = pd.Series(lg[other]).value_counts(normalize=True)
    ratios = pd.Series({g: pl.get(g, 0.0) / po[g] for g in cand})
    out["end_enrich_max"] = float(ratios.max()) if len(ratios) else float("nan")
    out["end_enrich_glyph"] = str(ratios.idxmax()) if len(ratios) else ""
    return out


def line_profile(df: pd.DataFrame, ctx, seed: int, level: str) -> dict:
    k = K[level]
    s = max(seed, 0)
    df = df.reset_index(drop=True)
    e = encode(df, ctx)
    out = {"n_tokens": len(df), "n_lines": len(e["starts"])}
    out.update(m_boundary(df, e, k, s))
    out.update(m_header(e, k, s))
    out.update(m_templates(e, k, s))
    out.update(m_same_length(e, k, s))
    out.update(m_vertical(df, k, s))
    out.update(m_line_end(df, e, ctx))
    return out


def signature_table(df: pd.DataFrame, ctx, source: str) -> pd.DataFrame:
    """20 most frequent class signatures per line length 4..8 of a full corpus (labels decoded)."""
    e = encode(df.reset_index(drop=True), ctx)
    labels = pd.Series(e["cls_label"]).groupby(e["cls"]).first()
    rows = []
    for n, s in _signatures(e["cls"], e["starts"], e["lens"]).items():
        tot = int(s.sum())
        for r, (sig, c) in enumerate(s.iloc[:20].items(), 1):
            rows.append({"source": source, "length": n, "rank": r, "count": int(c), "share": c / tot,
                         "n_lines": tot, "signature": " ".join(labels[int(x)] for x in sig.split())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def safe(source: str) -> str:
    return source.replace(":", "__")


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_profile(args) -> None:
    source = args.source
    kind = source.split(":")[0] if ":" in source and not source.startswith("ctl:") else "real"
    ctx = R9.Ctx(source if kind == "real" else f"real:{source.split(':')[1]}")
    if kind == "real":
        sig = OUT / "signatures" / f"{safe(source)}.csv"
        if not sig.exists():
            sig.parent.mkdir(parents=True, exist_ok=True)
            signature_table(ctx.real, ctx, source).to_csv(sig, index=False, float_format="%.5g")

        def one(seed: int) -> dict:
            df = ctx.real if seed < 0 else R9.page_subsample(ctx.real, seed)
            return {"source": source, "kind": "point" if seed < 0 else "sub",
                    **line_profile(df, ctx, seed, "point" if seed < 0 else "sub")}
        default = [-1] + list(range(N_SUB))
    else:
        def one(seed: int) -> dict:
            return {"source": source, "kind": "null", **line_profile(corpus(ctx, kind, seed), ctx, seed, "sub")}
        default = list(range(N_NULL))
    seeds = R9.parse_seeds(args.seeds, default)
    params = {"source": source, "prefixes": ctx.prefixes, "suffixes": ctx.suffixes, "glyph_units": ctx.gpreset,
              "n_classes": N_CLASSES, "top_stems": TOP_STEMS, "shuffles": K, "n_sub": N_SUB, "n_null": N_NULL,
              "sig_lengths": list(SIG_LENGTHS), "n_jaccard_pairs": N_JAC_PAIRS,
              "control": ctx.real.attrs.get("control")}
    run_seeds("010_line_as_item", seeds, one, OUT / "runs" / safe(source), params=params, progress=progress)


def all_sources() -> list[str]:
    out = list(REAL_TAGS) + [f"ctl:{c}" for c in CONTROLS if R9.control_available(c)]
    out += [f"{k}:{t}" for k, tags in NULLS.items() for t in tags]
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profile")
    p.add_argument("--source", required=True)
    p.add_argument("--seeds", default=None)
    sub.add_parser("sources")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        from summarize import summarize  # noqa: E402  (experiments/010_line_as_item/summarize.py)

        summarize()
        return
    if args.cmd == "sources":
        print("\n".join(all_sources()))
        return
    cmd_profile(args)


if __name__ == "__main__":
    main()
