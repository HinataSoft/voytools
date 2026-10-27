"""Experiment 007: position versus content.

(1) How is the information about affixes, stems and token length divided between position in the
line/paragraph, section (``illustration``), scribe (``hand``) and Currier language?  (2) Does a
dependence between tokens at distance d = 1..10 remain after positional effects are removed?

Usage (from the repository root; see ``spec.md`` and ``notes.md``)::

    # metric 1: information decomposition (real: point + 200 page subsamples; autocopy 86: 30 seeds)
    python experiments/007_position_vs_content/run.py decomp --file ZL3b-n.txt --space-policy split --source real
    python experiments/007_position_vs_content/run.py decomp --file ZL3b-n.txt --space-policy split --source autocopy
    # metric 2: residual MI(d); sources real, posclass, line, autocopy
    python experiments/007_position_vs_content/run.py mi --file ZL3b-n.txt --space-policy split --source real
    # metric 3: glyph-class drift along the line and down the paragraph (real, autocopy)
    python experiments/007_position_vs_content/run.py drift --file ZL3b-n.txt --space-policy split --source real
    # tables, plots
    python experiments/007_position_vs_content/run.py summarize
    # experiment 008: the canonical copying null (004, config 21) and the 005 template, written elsewhere
    python experiments/007_position_vs_content/run.py decomp --source selfcitation --out-root results/008_ts2020_recheck/007

Seeds: real subsamples 0..199 (80 % of pages without replacement; the full-text point estimate is
written separately); null sources 30 seeds.  The auto-copy null is config 86 of experiment 003
(trigram edits) with the seeds 8600..8629 of 003 phase B, so the samples are the ones of 003; they
are cached as token lists under ``scratch/`` (gitignored).  Every command is resumable.

Sources added for experiment 008 (``--out-root`` sends their runs elsewhere; the outputs of the
sources above do not change): ``selfcitation`` = ``selfcitation_like`` with config 21 of 004
(row ``ts2020_calB`` of ``results/004_selfcitation_port/best_configs.csv``), seeds 2100..2129 (the 004
phase-B seeds, so the samples are the ones of 004; Eva only, cached like the auto-copy samples);
``template`` = the full model of 005 (``ParadigmTemplate``, configuration of the ablation ``best``
of ``{ZL,GC}_split`` in ``results/005_structured_template_null/ablation/``), seeds 100000..100029 with
``pool_seed = seed - 100000`` (the 005 ablation seeds), learned on the decomposed real text of the set.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_tables, run_seeds, seed_range  # noqa: E402
from voytools._common import glyph_lists, line_ids, page_ids  # noqa: E402
from voytools.morphology import apply_bpe, decompose_table, learn_bpe  # noqa: E402
from voytools.nulls import as_table, shuffle_within_lines, shuffle_within_position_class  # noqa: E402
from voytools.selfcitation import selfcitation_like  # noqa: E402
from voytools.templates import ParadigmTemplate  # noqa: E402

DATA = ROOT / "data"
OUT = ROOT / "results" / "007_position_vs_content"
CACHE = ROOT / "scratch" / "007_cache"
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
POLICIES = ("split", "join")

N_SUB = 200            # real: page subsamples
SUB_FRAC = 0.8         # share of pages per subsample (without replacement)
N_NULL = 30            # seeds per null source
AC_CONFIG = 86         # experiment 003, trigram edits
AC_SEED0 = 8600        # 003 phase-B seeds of config 86: 8600..8649
TS_BEST = ROOT / "results" / "004_selfcitation_port" / "best_configs.csv"
TS_ROW = "ts2020_calB"  # config 21, the canonical copying null since 004
TS_SEED0 = 2100         # 004 phase-B seeds of config 21: 2100..2149
TPL_ABLATION = ROOT / "results" / "005_structured_template_null" / "ablation"
TPL_SEED0 = 100_000     # 005 ablation 'best': seed = 1000 * 100 + rep, pool_seed = rep
NEW_SOURCES = ("selfcitation", "template")  # experiment 008
DMAX = 10
TOP_STEMS = 200        # decomposition: stem class = top 200 stems + <other>
MIN_UNIT = 20          # MI: stems / affix frames with >= 20 occurrences + <other>
BPE_MERGES = 64
LEN_CAP = 10
K_BIAS = {"point": 50, "sub": 10, "null": 10}     # shuffles for the plug-in bias correction
K_PERM = 200                                       # permutation nulls of the point estimate
K_SHUF = {"point": 100, "sub": 20, "null": 20}     # shuffles per correction in MI(d)
TARGETS = ("prefix", "suffix", "stem", "length")
TERMS = (  # (name, X, Z) for I(Y; X | Z)
    ("I_pos", "pos", None), ("I_sec", "sec", None), ("I_hand", "hand", None), ("I_lang", "lang", None),
    ("I_sec_g_pos", "sec", "pos"), ("I_pos_g_sec", "pos", "sec"),
    ("I_sec_g_hand", "sec", "hand"), ("I_hand_g_sec", "hand", "sec"),
)
PAGE_LEVEL = {"sec", "hand", "lang"}
UNITS = ("edge", "bpe", "affix", "stem")
CORRECTIONS = {  # name -> (class, scope)
    "pc": ("spec", "paragraph"),      # the spec's null: position class inside the paragraph
    "pcpage": ("spec", "page"),       # sensitivity: same classes, whole page
    "pcfine": ("fine", "page"),       # metric 3 covariate: fine bins, whole page
    "line": (None, "line"),           # within-line shuffle
}
GLYPH_TOP = 12
DRIFT_AXES = ("rel_pos", "rel_pos_mid", "line_in_par")


def parse_seeds(spec: str) -> list[int]:
    """``'0:20'`` -> seeds 0..19; ``'0:100/4:1'`` -> chunk 1 of 4 of 0..99; ``'3,5,8'`` -> those seeds."""
    if "/" in spec:
        rng, part = spec.split("/")
        chunks, k = part.split(":")
        a, b = rng.split(":")
        return seed_range(int(a), int(b), int(chunks))[int(k)]
    if ":" in spec:
        a, b = spec.split(":")
        return list(range(int(a), int(b)))
    return [int(s) for s in spec.split(",") if s.strip()]


def tag_of(file: str, policy: str) -> str:
    return f"{'GC' if file.startswith('GC') else 'ZL'}_{policy}"


def _load_003():
    spec = importlib.util.spec_from_file_location("exp003", ROOT / "experiments" / "003_autocopy_strict" / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def selfcitation_params() -> dict:
    """``selfcitation_like`` keyword arguments of config 21 (row ``ts2020_calB`` of 004's best_configs.csv)."""
    r = pd.read_csv(TS_BEST).set_index("source").loc[TS_ROW]
    return {"p_copy": float(r["p_copy"]), "replace_weights": tuple(float(x) for x in str(r["replace_weights"]).split("/")),
            "p_same_position": float(r["p_same_position"])}


def template_config(key: str) -> dict:
    """Full 005 template model of the file ``key`` ('ZL'/'GC'): configuration of the split ablation ``best``."""
    meta = json.loads((TPL_ABLATION / f"{key}_split" / "best" / "part0" / "meta.json").read_text())
    c = meta["params"]["config"]
    return {"n_stems": int(c["inventory_n_stems"]), "zipf_s": float(c["zipf_s"]), "paradigm": c["paradigm"],
            "prefix_dependence": c["prefix_dependence"], "topic_concentration": float(c["topic_concentration"]),
            "page_locality": float(c["page_locality"]), "suffix_by": c["suffix_by"]}


# --------------------------------------------------------------------------- #
# position classes and groups
# --------------------------------------------------------------------------- #


def position_classes(df: pd.DataFrame, fine: bool = False) -> np.ndarray:
    """Spec classes (= default of ``shuffle_within_position_class``) or the fine variant.

    spec: line place (first/mid/last/single) x rel_pos in 5 bins x line_in_par 0 / 1 / 2+;
    fine: line place x rel_pos in 10 bins x line_in_par 0..5 / 6+.
    """
    is_first = df["is_first"].to_numpy(dtype=bool)
    is_last = df["is_last"].to_numpy(dtype=bool)
    place = np.where(is_first & is_last, "single", np.where(is_first, "first", np.where(is_last, "last", "mid")))
    bins = 10 if fine else 5
    rel = np.clip((df["rel_pos"].to_numpy(dtype=float) * bins).astype(np.int64), 0, bins - 1)
    lip = pd.to_numeric(df["line_in_par"], errors="coerce").to_numpy(dtype=float)
    cap = 6 if fine else 2
    lip = np.where(np.isnan(lip), cap, np.minimum(lip, cap)).astype(np.int64)
    return np.asarray([f"{a}|{b}|{c}" for a, b, c in zip(place, rel, lip)], dtype=object)


def paragraph_units(df: pd.DataFrame) -> np.ndarray:
    """Paragraph id (page + par_index; a line without paragraph index is its own unit)."""
    par = df["par_index"].tolist()
    keys = [f"p{pg}|{pi}" if not pd.isna(pi) else f"l{li}" for pg, pi, li in zip(page_ids(df), par, line_ids(df))]
    return pd.factorize(pd.Series(keys, dtype=object))[0].astype(np.int64)


def codes(values) -> tuple[np.ndarray, int]:
    c, u = pd.factorize(pd.Series(list(values), dtype=object))
    return c.astype(np.int64), len(u)


def combine(*arrays: np.ndarray) -> np.ndarray:
    """Joint code of several integer arrays (re-factorised, so it stays small)."""
    key = arrays[0].astype(np.int64)
    for a in arrays[1:]:
        key = key * (int(a.max()) + 1 if len(a) else 1) + a
    return pd.factorize(key)[0].astype(np.int64)


class Shuffler:
    """Permutations inside groups, identical to ``voytools.nulls._permute_within``."""

    def __init__(self, groups: np.ndarray):
        self.groups = groups
        self.base = np.argsort(groups, kind="stable")

    def src(self, rng: np.random.Generator) -> np.ndarray:
        n = len(self.groups)
        order = np.lexsort((rng.random(n), self.groups))
        out = np.empty(n, dtype=np.int64)
        out[self.base] = order
        return out


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #


class Context:
    """Real text of one file / space policy with the inventories that null corpora reuse."""

    def __init__(self, file: str, policy: str):
        self.file, self.policy, self.tag = file, policy, tag_of(file, policy)
        self.real = filter_tokens(load_tables(DATA / file, space_policy=policy).tokens, generic_type="P", clean=True)
        self.real = self.real.reset_index(drop=True)
        eva = self.real.attrs.get("glyph_preset") is not None
        self.gpreset = "extended" if eva else None  # glyph-level units (CLAUDE.md)
        d = decompose_table(self.real)
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]
        self.merges = learn_bpe(glyph_lists(self.real, self.gpreset), BPE_MERGES)
        self._ac = None
        self._tpl = None

    # ---- null corpora
    def autocopy(self, seed: int) -> pd.DataFrame:
        path = CACHE / self.tag / f"ac_{seed}.txt"
        if path.exists():
            flat = path.read_text(encoding="utf-8").split("\n")
        else:
            if self._ac is None:
                mod = _load_003()
                self._ac = (mod, mod.Context(self.file, self.policy))
            mod, ctx3 = self._ac
            flat = ctx3.generate(mod.config_params(AC_CONFIG), seed)["token"].tolist()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("\n".join(flat), encoding="utf-8")
        lengths = np.bincount(line_ids(self.real))
        lines, pos = [], 0
        for n in lengths:
            lines.append(flat[pos: pos + n])
            pos += n
        return as_table(lines, self.real)

    def selfcitation(self, seed: int) -> pd.DataFrame:
        if self.gpreset is None:
            raise ValueError("selfcitation_like generates basic Eva; no null for a non-Eva file")
        path = CACHE / self.tag / f"ts_{seed}.txt"
        if path.exists():
            lines = [ln.split(" ") if ln else [] for ln in path.read_text(encoding="utf-8").split("\n")]
            return as_table(lines, self.real, glyph_preset=self.real.attrs.get("glyph_preset"))
        x = selfcitation_like(self.real, seed=seed, **selfcitation_params())
        lines = x.groupby(line_ids(x), sort=True)["token"].apply(list).tolist()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(" ".join(ln) for ln in lines), encoding="utf-8")
        return x

    def template(self, seed: int) -> pd.DataFrame:
        path = CACHE / self.tag / f"tpl_{seed}.txt"
        if path.exists():
            lines = [ln.split(" ") if ln else [] for ln in path.read_text(encoding="utf-8").split("\n")]
            return as_table(lines, self.real, glyph_preset=self.real.attrs.get("glyph_preset"))
        if self._tpl is None:
            dfd = decompose_table(self.real, prefixes=self.prefixes, suffixes=self.suffixes)
            self._tpl = ParadigmTemplate.learn(dfd, **template_config(self.tag.split("_")[0]))
        x = self._tpl.sample(self.real, seed=seed, pool_seed=seed - TPL_SEED0)
        lines = x.groupby(line_ids(x), sort=True)["token"].apply(list).tolist()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(" ".join(ln) for ln in lines), encoding="utf-8")
        return x

    def corpus(self, source: str, seed: int) -> pd.DataFrame:
        if source == "real":
            return self.real
        if source == "autocopy":
            return self.autocopy(seed)
        if source == "selfcitation":
            return self.selfcitation(seed)
        if source == "template":
            return self.template(seed)
        if source == "posclass":
            return shuffle_within_position_class(self.real, seed=seed)
        if source == "line":
            return shuffle_within_lines(self.real, seed=seed)
        raise ValueError(source)

    # ---- encoding
    def encode(self, df: pd.DataFrame) -> dict:
        """Integer codes of targets, units and covariates of one corpus (categories from the corpus)."""
        dfd = decompose_table(df, prefixes=self.prefixes, suffixes=self.suffixes)
        pre = dfd["prefix"].to_numpy(dtype=object)
        stem = dfd["stem"].to_numpy(dtype=object)
        suf = dfd["suffix"].to_numpy(dtype=object)
        sc = pd.Series(stem).value_counts()
        top = set(sc.index[:TOP_STEMS])
        frequent = set(sc.index[sc >= MIN_UNIT])
        frame = pd.Series(pre + "|" + suf)
        fc = frame.value_counts()
        frame_ok = set(fc.index[fc >= MIN_UNIT])
        gl = glyph_lists(df, self.gpreset)
        bpe = apply_bpe(gl, self.merges)
        enc = {
            "n": len(df),
            "Y": {
                "prefix": codes(pre),
                "suffix": codes(suf),
                "stem": codes([s if s in top else "<other>" for s in stem]),
                "length": codes(np.minimum(df["len_glyphs"].to_numpy(dtype=np.int64), LEN_CAP)),
            },
            "U": {
                "edge": (codes([g[-1] if g else "" for g in gl]), codes([g[0] if g else "" for g in gl])),
                "bpe": (codes([u[-1] if u else "" for u in bpe]), codes([u[0] if u else "" for u in bpe])),
                "affix": (codes([f if f in frame_ok else "<other>" for f in frame]),) * 2,
                "stem": (codes([s if s in frequent else "<other>" for s in stem]),) * 2,
            },
            "X": {
                "pos": codes(position_classes(df)),
                "sec": codes(df["illustration"].fillna("?")),
                "hand": codes(df["hand"].fillna("?")),
                "lang": codes(df["language"].fillna("?")),
            },
            "page": page_ids(df),
            "para": paragraph_units(df),
            "line": line_ids(df),
            "pc_spec": position_classes(df),
            "pc_fine": position_classes(df, fine=True),
            "glyphs": gl,
            "drift": {
                "rel_pos": np.clip((df["rel_pos"].to_numpy(dtype=float) * 5).astype(np.int64), 0, 4),
                "line_in_par": np.minimum(pd.to_numeric(df["line_in_par"], errors="coerce").fillna(5)
                                          .to_numpy(dtype=np.int64), 5),
                "mid": ~(df["is_first"].to_numpy(dtype=bool) | df["is_last"].to_numpy(dtype=bool)),
            },
        }
        return enc


def subset(enc: dict, mask: np.ndarray) -> dict:
    """Restrict an encoding to the rows in ``mask`` (whole pages, so paragraphs stay intact)."""
    def re(cn):
        return cn[0][mask], cn[1]

    out = {
        "n": int(mask.sum()),
        "Y": {k: re(v) for k, v in enc["Y"].items()},
        "U": {k: (re(v[0]), re(v[1])) for k, v in enc["U"].items()},
        "X": {k: re(v) for k, v in enc["X"].items()},
        "glyphs": [g for g, m in zip(enc["glyphs"], mask) if m],
        "drift": {k: v[mask] for k, v in enc["drift"].items()},
    }
    for k in ("page", "para", "line", "pc_spec", "pc_fine"):
        out[k] = enc[k][mask]
    return out


def page_subsample(enc: dict, seed: int) -> np.ndarray:
    pages = np.unique(enc["page"])
    rng = np.random.default_rng(seed)
    keep = rng.choice(pages, size=int(round(SUB_FRAC * len(pages))), replace=False)
    return np.isin(enc["page"], keep)


# --------------------------------------------------------------------------- #
# metric 1: information decomposition
# --------------------------------------------------------------------------- #


def H(c: np.ndarray) -> float:
    n = len(c)
    if n == 0:
        return float("nan")
    p = np.bincount(c) / n
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def cmi(y: np.ndarray, x: np.ndarray, z: np.ndarray | None) -> float:
    """Plug-in I(Y; X | Z) in bits (Z = None: unconditional)."""
    if z is None:
        return H(y) + H(x) - H(combine(x, y))
    return H(combine(y, z)) + H(combine(x, z)) - H(combine(x, y, z)) - H(z)


def decomposition(enc: dict, k_bias: int, seed: int) -> dict:
    """H(Y) and the TERMS for every target: raw, shuffle bias (Y permuted within Z), corrected."""
    rng = np.random.default_rng(seed)
    out = {}
    X = {k: v[0] for k, v in enc["X"].items()}
    zero = np.zeros(enc["n"], dtype=np.int64)
    shufflers = {z: Shuffler(zero if z is None else X[z]) for z in {t[2] for t in TERMS}}
    for yname in TARGETS:
        y = enc["Y"][yname][0]
        out[f"{yname}__H"] = H(y)
        for name, xn, zn in TERMS:
            z = None if zn is None else X[zn]
            raw = cmi(y, X[xn], z)
            sims = [cmi(y[shufflers[zn].src(rng)], X[xn], z) for _ in range(k_bias)]
            bias = float(np.mean(sims))
            out[f"{yname}__{name}__raw"] = raw
            out[f"{yname}__{name}__bias"] = bias
            out[f"{yname}__{name}__corr"] = raw - bias
    return out


def permutation_nulls(enc: dict, k: int, seed: int) -> pd.DataFrame:
    """Raw terms under two nulls: Y permuted inside pages; page-level labels permuted among pages."""
    rng = np.random.default_rng(seed)
    X = {kk: v[0] for kk, v in enc["X"].items()}
    page = enc["page"]
    within_page = Shuffler(page)
    upages, first = np.unique(page, return_index=True)
    pidx = np.searchsorted(upages, page)
    rows = []
    for r in range(k):
        src = within_page.src(rng)
        perm = rng.permutation(len(upages))
        Xp = dict(X)
        for v in PAGE_LEVEL:  # one page's (section, hand, language) record moves to another page
            Xp[v] = X[v][first][perm][pidx]
        for yname in TARGETS:
            y = enc["Y"][yname][0]
            for name, xn, zn in TERMS:
                rows.append({"rep": r, "Y": yname, "term": name, "null": "within_page",
                             "value": cmi(y[src], X[xn], None if zn is None else X[zn])})
                if xn in PAGE_LEVEL or zn in PAGE_LEVEL:
                    rows.append({"rep": r, "Y": yname, "term": name, "null": "page_labels",
                                 "value": cmi(y, Xp[xn], None if zn is None else Xp[zn])})
    return pd.DataFrame(rows)


def cmd_decomp(args) -> None:
    ctx = Context(args.file, args.space_policy)
    base = OUT / "decomp" / ctx.tag
    if args.source == "real":
        enc = ctx.encode(ctx.real)
        point_path = base / "point.csv"
        if not point_path.exists():
            t0 = time.time()
            pt = decomposition(enc, K_BIAS["point"], seed=0)
            perm = permutation_nulls(enc, K_PERM, seed=1)
            rows = []
            for yname in TARGETS:
                for name, _, _ in TERMS:
                    row = {"Y": yname, "term": name, "H": pt[f"{yname}__H"],
                           **{c: pt[f"{yname}__{name}__{c}"] for c in ("raw", "bias", "corr")}}
                    for nl in ("within_page", "page_labels"):
                        v = perm.query("Y == @yname and term == @name and null == @nl")["value"].to_numpy()
                        if len(v):
                            row[f"{nl}_mean"] = float(v.mean())
                            row[f"{nl}_q975"] = float(np.quantile(v, 0.975))
                            row[f"{nl}_p"] = float((1 + (v >= row["raw"] - 1e-12).sum()) / (1 + len(v)))
                    rows.append(row)
            base.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(rows).to_csv(point_path, index=False)
            print(f"point estimate {ctx.tag}: {time.time() - t0:.0f} s", flush=True)

        def one(seed: int) -> dict:
            return {"tag": ctx.tag, "source": "real_sub", **decomposition(subset(enc, page_subsample(enc, seed)),
                                                                       K_BIAS["sub"], seed)}

        seeds = parse_seeds(args.seeds) if args.seeds else list(range(N_SUB))
        out = base / "real_sub"
    else:
        def one(seed: int) -> dict:
            return {"tag": ctx.tag, "source": args.source,
                    **decomposition(ctx.encode(ctx.corpus(args.source, seed)), K_BIAS["null"], seed)}

        seeds = parse_seeds(args.seeds) if args.seeds else null_seeds(args.source)
        out = base / args.source
    run_seeds("007_decomp", seeds, one, out, params=params_of(args, ctx), progress=progress)


# --------------------------------------------------------------------------- #
# metric 2: MI(d)
# --------------------------------------------------------------------------- #


def mi_pairs(x: np.ndarray, y: np.ndarray, kx: int, ky: int) -> float:
    n = len(x)
    if n == 0:
        return float("nan")
    joint = np.bincount(x * ky + y, minlength=kx * ky)
    nx = np.bincount(x, minlength=kx).astype(float)
    ny = np.bincount(y, minlength=ky).astype(float)
    nz = np.nonzero(joint)[0]
    c = joint[nz].astype(float)
    return float((c / n * np.log2(c * n / (nx[nz // ky] * ny[nz % ky]))).sum())


def mi_curve(enc: dict, same: dict, src: np.ndarray | None = None) -> dict:
    out = {}
    for u in UNITS:
        (L, kl), (R, kr) = enc["U"][u]
        if src is not None:
            L, R = L[src], R[src]
        for d in range(1, DMAX + 1):
            m = same[d]
            out[(u, d)] = mi_pairs(L[:-d][m], R[d:][m], kl, kr)
    return out


def correction_groups(enc: dict) -> dict:
    groups = {}
    for name, (cls, scope) in CORRECTIONS.items():
        if scope == "line":
            groups[name] = enc["line"]
            continue
        unit = enc["para"] if scope == "paragraph" else enc["page"]
        labels = enc["pc_spec"] if cls == "spec" else enc["pc_fine"]
        groups[name] = pd.factorize(pd.Series([f"{a}#{b}" for a, b in zip(unit, labels)], dtype=object))[0]
    return groups


def mi_stats(enc: dict, k: int, seed: int) -> dict:
    para = enc["para"]
    same = {d: para[:-d] == para[d:] for d in range(1, DMAX + 1)}
    rng = np.random.default_rng(seed)
    out = {f"n_pairs_d{d}": int(same[d].sum()) for d in (1, DMAX)}
    raw = mi_curve(enc, same)
    for (u, d), v in raw.items():
        out[f"{u}_d{d}_raw"] = v
    for name, g in correction_groups(enc).items():
        sh = Shuffler(np.asarray(g, dtype=np.int64))
        acc = {key: 0.0 for key in raw}
        for _ in range(k):
            for key, v in mi_curve(enc, same, sh.src(rng)).items():
                acc[key] += v / k
        for (u, d), v in acc.items():
            out[f"{u}_d{d}_sh_{name}"] = v
        out[f"frozen_{name}"] = float(np.mean(sh.src(rng) == np.arange(len(g))))
    return out


def cmd_mi(args) -> None:
    ctx = Context(args.file, args.space_policy)
    base = OUT / "mi" / ctx.tag
    if args.source == "real":
        enc = ctx.encode(ctx.real)
        # the run.py groups must reproduce the package null exactly
        sh = Shuffler(correction_groups(enc)["pc"]).src(np.random.default_rng(5))
        assert (ctx.real["token"].to_numpy()[sh] == shuffle_within_position_class(ctx.real, seed=5)["token"]
                .to_numpy()).all()
        point_path = base / "real_point"
        if not (point_path / "runs.csv").exists():
            run_seeds("007_mi_point", [-1], lambda s: {"tag": ctx.tag, "source": "real",
                                                       **mi_stats(enc, K_SHUF["point"], 0)},
                      point_path, params=params_of(args, ctx), progress=progress)

        def one(seed: int) -> dict:
            return {"tag": ctx.tag, "source": "real_sub",
                    **mi_stats(subset(enc, page_subsample(enc, seed)), K_SHUF["sub"], seed)}

        seeds = parse_seeds(args.seeds) if args.seeds else list(range(N_SUB))
        out = base / "real_sub"
    else:
        def one(seed: int) -> dict:
            return {"tag": ctx.tag, "source": args.source,
                    **mi_stats(ctx.encode(ctx.corpus(args.source, seed)), K_SHUF["null"], seed)}

        seeds = parse_seeds(args.seeds) if args.seeds else null_seeds(args.source)
        out = base / args.source
    run_seeds("007_mi", seeds, one, out, params=params_of(args, ctx), progress=progress)


# --------------------------------------------------------------------------- #
# metric 3: drift
# --------------------------------------------------------------------------- #


def drift_stats(enc: dict, classes: list[str], seed: int, k_bias: int) -> dict:
    """Glyph-class shares per bin and I(glyph class; bin) (token-level shuffle-corrected)."""
    cl = {g: i for i, g in enumerate(classes)}
    other = len(classes)
    tok_of = np.repeat(np.arange(len(enc["glyphs"])), [len(g) for g in enc["glyphs"]])
    gcode = np.fromiter((cl.get(x, other) for g in enc["glyphs"] for x in g), dtype=np.int64, count=len(tok_of))
    names = classes + ["<other>"]
    rng = np.random.default_rng(seed)
    out = {}
    for axis in DRIFT_AXES:
        b = enc["drift"]["line_in_par" if axis == "line_in_par" else "rel_pos"]
        keep_tok = enc["drift"]["mid"] if axis == "rel_pos_mid" else np.ones(len(b), dtype=bool)
        sel = keep_tok[tok_of]
        g, bt = gcode[sel], b[tok_of][sel]
        nb = int(b.max()) + 1
        tab = np.zeros((nb, len(names)))
        np.add.at(tab, (bt, g), 1)
        for i in range(nb):
            tot = tab[i].sum()
            for j, nm in enumerate(names):
                out[f"{axis}__{i}__{nm}"] = tab[i, j] / tot if tot else float("nan")
        raw = H(g) + H(bt) - H(combine(bt, g))
        # bias: bins permuted among tokens (all glyphs of a token keep one bin)
        tok_sel = np.nonzero(keep_tok)[0]
        sims = []
        for _ in range(k_bias):
            pb = b.copy()
            pb[tok_sel] = b[tok_sel][rng.permutation(len(tok_sel))]
            sims.append(H(g) + H(pb[tok_of][sel]) - H(combine(pb[tok_of][sel], g)))
        out[f"{axis}__I_raw"] = raw
        out[f"{axis}__I_corr"] = raw - float(np.mean(sims))
    return out


def cmd_drift(args) -> None:
    ctx = Context(args.file, args.space_policy)
    enc_real = ctx.encode(ctx.real)
    flat = pd.Series([x for g in enc_real["glyphs"] for x in g]).value_counts()
    classes = flat.index[:GLYPH_TOP].tolist()
    base = OUT / "drift" / ctx.tag
    if args.source == "real":
        point = base / "real_point"
        if not (point / "runs.csv").exists():
            run_seeds("007_drift_point", [-1], lambda s: {"tag": ctx.tag, "source": "real",
                                                          **drift_stats(enc_real, classes, 0, K_BIAS["point"])},
                      point, params=params_of(args, ctx, glyph_classes=classes))

        def one(seed: int) -> dict:
            return {"tag": ctx.tag, "source": "real_sub",
                    **drift_stats(subset(enc_real, page_subsample(enc_real, seed)), classes, seed, K_BIAS["sub"])}

        seeds = parse_seeds(args.seeds) if args.seeds else list(range(N_SUB))
        out = base / "real_sub"
    else:
        def one(seed: int) -> dict:
            return {"tag": ctx.tag, "source": args.source,
                    **drift_stats(ctx.encode(ctx.corpus(args.source, seed)), classes, seed, K_BIAS["null"])}

        seeds = parse_seeds(args.seeds) if args.seeds else null_seeds(args.source)
        out = base / args.source
    run_seeds("007_drift", seeds, one, out, params=params_of(args, ctx, glyph_classes=classes), progress=progress)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def null_seeds(source: str) -> list[int]:
    seed0 = {"autocopy": AC_SEED0, "selfcitation": TS_SEED0, "template": TPL_SEED0}.get(source, 0)
    return list(range(seed0, seed0 + N_NULL))


def params_of(args, ctx: Context, **extra) -> dict:
    if args.source == "selfcitation":
        extra = {**extra, "selfcitation_004": {"row": TS_ROW, **{k: (list(v) if isinstance(v, tuple) else v)
                                                                  for k, v in selfcitation_params().items()}}}
    if args.source == "template":
        extra = {**extra, "template_005": template_config(ctx.tag.split("_")[0]), "pool_seed": f"seed - {TPL_SEED0}"}
    return {"command": args.cmd, "file": ctx.file, "space_policy": ctx.policy, "source": args.source,
            "prefixes": ctx.prefixes, "suffixes": ctx.suffixes, "glyph_preset_units": ctx.gpreset,
            "n_sub": N_SUB, "sub_frac": SUB_FRAC, "n_null": N_NULL, "autocopy_config_003": AC_CONFIG,
            "dmax": DMAX, "top_stems": TOP_STEMS, "min_unit": MIN_UNIT, "bpe_merges": BPE_MERGES,
            "k_bias": K_BIAS, "k_perm": K_PERM, "k_shuf": K_SHUF,
            "corrections": {k: list(v) for k, v in CORRECTIONS.items()}, **extra}


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("decomp", "mi", "drift"):
        p = sub.add_parser(name)
        p.add_argument("--file", default="ZL3b-n.txt")
        p.add_argument("--space-policy", default="split", choices=POLICIES)
        p.add_argument("--source", default="real", choices=("real", "autocopy", "posclass", "line") + NEW_SOURCES)
        p.add_argument("--seeds", default=None)
        p.add_argument("--out-root", default=None, help="write under this folder instead of results/007_...")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    global OUT
    if getattr(args, "out_root", None):
        OUT = Path(args.out_root).resolve()
    if args.cmd == "summarize":
        from summarize import summarize  # noqa: E402  (experiments/007_position_vs_content/summarize.py)

        summarize()
        return
    {"decomp": cmd_decomp, "mi": cmd_mi, "drift": cmd_drift}[args.cmd](args)


if __name__ == "__main__":
    main()
