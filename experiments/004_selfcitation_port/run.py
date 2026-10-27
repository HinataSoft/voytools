"""Experiment 004: port of the Timm & Schinner self-citation generator as the canonical copying null.

Usage (from the repository root; see ``spec.md``)::

    # reference repo of Rozanova & Temerev (MIT), cloned once into a temporary folder
    git clone --depth 1 https://github.com/lrozanova/voynich-units /tmp/voynich-units
    # real text: point estimates of all metrics (seed 0), ZL split and join
    python experiments/004_selfcitation_port/run.py real --space-policy split
    # step 1: validation of the port (identity with their samples, metrics vs their samples and tab. 7)
    python experiments/004_selfcitation_port/run.py validate --rt-repo /tmp/voynich-units
    # step 2: calibration grids A and B (one run set, 4 parallel chunks), then ranking
    python experiments/004_selfcitation_port/run.py grid --seeds 0:6000/4:0 --out results/004_selfcitation_port/grid/part0
    python experiments/004_selfcitation_port/run.py rank      # -> grid_A.csv, grid_B.csv, phaseB_choice.json
    # step 2, phase B: chosen configs x 50 seeds, all metrics, split and join
    python experiments/004_selfcitation_port/run.py phaseB --space-policy split
    # 003 reference (trigram config 86) with the paradigm metrics that 003 did not compute
    python experiments/004_selfcitation_port/run.py trigram86 --space-policy split
    # step 3: vocabulary size as a confound
    python experiments/004_selfcitation_port/run.py vocab
    # step 4: tables, profile, plot
    python experiments/004_selfcitation_port/run.py summarize

Grid: ``config_id`` indexes ``itertools.product`` over :data:`AXES` (60 configurations); a run has
``seed = config_id * 100 + rep``.  Calibration A uses the configurations with
``p_same_position == 28`` and ``rep < 5``; calibration B all configurations with ``rep < 3`` (both read
from the same run set).  The replacement weights ``(0, 0, 100)`` are not in the spec grid: they are the
values with which Rozanova & Temerev produced their published samples (found in the validation step),
added as an out-of-spec row (``in_spec_grid = 0``) and never chosen as ``calA`` / ``calB``.
The generator works on Eva only; every metric is computed on ZL3b (split and join), GC2a is not used.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_runs, load_tables, merge_runs, run_seeds, seed_range  # noqa: E402
from voytools.morphology import (  # noqa: E402
    bpe_dependence_curve,
    decompose_table,
    inventory_jaccard,
    paradigm_metrics,
    productivity_slope,
    stem_dispersion,
)
from voytools.nulls import as_table, autocopy_like  # noqa: E402
from voytools.selfcitation import selfcitation_generate, selfcitation_like, selfcitation_template  # noqa: E402
from voytools.stats import (  # noqa: E402
    adjacent_token_mi,
    edge_glyph_mi,
    glyph_entropy_report,
    neighbor_similarity,
    positional_profile,
    vocabulary_summary,
)
from voytools.wordgrammar import learn_word_grammar  # noqa: E402

DATA = ROOT / "data"
FILE = "ZL3b-n.txt"
OUT = ROOT / "results" / "004_selfcitation_port"
REF001 = ROOT / "results" / "001_autocopy_calibration"
REF002 = ROOT / "results" / "002_paradigm_coherence"
REF003 = ROOT / "results" / "003_autocopy_strict"

RT_WEIGHTS = (0, 0, 100)  # weights of the published R&T samples (validation step)
AXES = {
    "p_copy": (0.0, 0.1, 0.2, 0.3),
    "replace_weights": ((30, 50, 20), (50, 35, 15), (20, 50, 30), (60, 30, 10), RT_WEIGHTS),
    "p_same_position": (14, 28, 42),
}
CONFIGS = list(itertools.product(*AXES.values()))  # 60
REPS_A, REPS_B, REPS_PHASE_B = 5, 3, 50
JAVA_DEFAULT = {"p_copy": 0.0, "replace_weights": (30, 50, 20), "p_same_position": 28}
RT_PUBLISHED = {"p_copy": 0.0, "replace_weights": RT_WEIGHTS, "p_same_position": 28}

TARGETS = {"share_near": 2.0, "hapax_share": 2.0, "H2": 2.0, "n_types": 1.0, "mean_len": 1.0}  # loss B (001)
LOSS_A_SCALE = {"mean_len": 0.25, "n_types": 500.0, "hapax_share": 0.03}  # loss A (R&T calibrate())
COMPONENT_LIMIT = 0.05
LOSS_LIMIT = 0.15
EXTRA = ["near_ratio"]
NONTARGET = ["edge_mi_excess", "token_mi_excess_share", "bpe_gap64", "pos_jsd_first_mid", "cmi_excess",
             "frac_stems_multi_suffix", "mean_suffixes_per_stem", "frac_multi_suffix_5to9", "H_suffix_given_stem",
             "mi_prefix_suffix", "productivity_slope", "stem_dispersion", "jaccard_stem_sections",
             "jsd_stem_sections", "jsd_suffix_AB"]
CHECKS = ["n_stems_k", "top50_coverage"]
PARADIGM = ["cmi_excess", "frac_stems_multi_suffix", "mean_suffixes_per_stem", "frac_multi_suffix_5to9",
            "H_suffix_given_stem", "mi_prefix_suffix", "productivity_slope", "stem_dispersion",
            "jaccard_stem_sections", "jsd_stem_sections", "n_stems_k"]
MIN_STEM, MIN_INV = 5, 3  # as in 002

# docs/reference_values.md, Rozanova & Temerev (2026) tab. 7 (ZL3b, paragraph text, split-like)
REFERENCE = {
    "ref_RT2026_voynich": {"H2": (2.69, 2.69), "bpe_gap64": (1.05, 1.05), "token_mi_excess_share": (0.0079, 0.0079),
                           "hapax_share": (0.70, 0.70), "edge_mi_excess": (0.197, 0.197), "near_ratio": (1.12, 1.12)},
    "ref_RT2026_selfcitation": {"H2": (2.6, 2.7), "bpe_gap64": (0.65, 0.86), "token_mi_excess_share": (0.005, 0.006),
                                "hapax_share": (0.59, 0.60), "edge_mi_excess": (0.006, 0.009),
                                "near_ratio": (1.06, 1.23)},
}
# spec criterion 1: tolerances of the validation
TOLERANCE = {"H2": 0.1, "hapax_share": 0.05, "bpe_gap64": 0.1, "edge_mi_excess": 0.01, "token_mi_excess_share": 0.003}


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


def config_params(cid: int) -> dict:
    return dict(zip(AXES, CONFIGS[cid]))


def config_id(params: dict) -> int:
    return CONFIGS.index(tuple(params[k] for k in AXES))


def config_columns(cid: int) -> dict:
    p = config_params(cid)
    return {"config_id": cid, "p_copy": p["p_copy"], "replace_weights": "/".join(map(str, p["replace_weights"])),
            "p_same_position": p["p_same_position"], "in_spec_grid": int(p["replace_weights"] != RT_WEIGHTS)}


def reps_of(cid: int) -> int:
    return REPS_A if config_params(cid)["p_same_position"] == 28 else REPS_B


# --------------------------------------------------------------------------- #
# data and metrics
# --------------------------------------------------------------------------- #


class Context:
    """Real ZL paragraph text of one space policy and its learned affix lists."""

    def __init__(self, policy: str):
        self.policy = policy
        self.real = filter_tokens(load_tables(DATA / FILE, space_policy=policy).tokens, generic_type="P", clean=True)
        d = decompose_table(self.real)
        self.prefixes, self.suffixes = d.attrs["prefixes"], d.attrs["suffixes"]

    def generate(self, cid: int, seed: int) -> pd.DataFrame:
        return selfcitation_like(self.real, seed=seed, **config_params(cid))


def lines_of(df: pd.DataFrame) -> list[list[str]]:
    return df.groupby("locus_id", sort=False)["token"].apply(list).tolist()


def _mean_pairs(dfd: pd.DataFrame, by: str, column: str, groups=None) -> tuple[float, float]:
    t = inventory_jaccard(dfd, by, column=column, min_count=MIN_INV, groups=groups)
    return float(t["jaccard"].mean()), float(t["jsd"].mean())


def paradigm_block(df: pd.DataFrame, prefixes, suffixes, seed: int, sections: bool = True) -> dict:
    """002 metric set (identical code and thresholds)."""
    dfd = decompose_table(df, prefixes=prefixes, suffixes=suffixes)
    pm = paradigm_metrics(dfd, min_stem_count=MIN_STEM, seed=seed)
    s = pm.summary
    low = pm.per_stem[pm.per_stem["count"] < 2 * MIN_STEM]
    out = {k: float(s[k]) for k in ("n_stems_k", "mean_suffixes_per_stem", "frac_stems_multi_suffix",
                                     "H_suffix_given_stem", "mi_prefix_suffix", "cmi_excess")}
    out["frac_multi_suffix_5to9"] = float((low["n_suffixes"] > 1).mean()) if len(low) else float("nan")
    out["productivity_slope"] = productivity_slope(dfd, min_stem_count=MIN_STEM)["slope"]
    out["stem_dispersion"] = stem_dispersion(dfd, min_stem_count=MIN_STEM)
    if sections and dfd["illustration"].notna().any():
        out["jaccard_stem_sections"], out["jsd_stem_sections"] = _mean_pairs(dfd, "illustration", "stem")
        _, out["jsd_suffix_AB"] = _mean_pairs(dfd, "language", "suffix", ["A", "B"])
    return out


def stats_for(df: pd.DataFrame, seed: int, full: bool, ctx: Context | None = None) -> dict:
    """Target metrics (001/003 code) and R&T profile metrics; ``full`` adds BPE, position, paradigms."""
    near = neighbor_similarity(df, window=5, shuffles=20, seed=seed).iloc[0]
    vs = vocabulary_summary(df).iloc[0]
    ge = glyph_entropy_report(df, preset="extended", scope="all").iloc[0]
    out = {
        "share_near": float(near["share_near"]),
        "near_ratio": float(near["ratio"]),
        "n_types": float(vs["n_types"]),
        "hapax_share": float(vs["hapax_share"]),
        "mean_len": float(vs["mean_len"]),
        "H2": float(ge["H2"]),
        "top50_coverage": float(vs["top50_coverage"]),
        "edge_mi_excess": float(edge_glyph_mi(df, shuffles=20, seed=seed).iloc[0]["excess_bits"]),
        "token_mi_excess_share": float(adjacent_token_mi(df, cap=2000, shuffles=20, seed=seed).iloc[0]["excess_share"]),
    }
    if full:
        curve = bpe_dependence_curve(lines_of(df), merge_checkpoints=(0, 32, 64, 128)).set_index("merges")["gap"]
        out["bpe_gap64"] = float(curve.get(64, np.nan))
        out["pos_jsd_first_mid"] = float(positional_profile(df, which="first").jsd.loc["line_first", "mid"])
        if ctx is not None:
            out.update(paradigm_block(df, ctx.prefixes, ctx.suffixes, seed))
    return out


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


# --------------------------------------------------------------------------- #
# real text
# --------------------------------------------------------------------------- #


def cmd_real(args) -> None:
    ctx = Context(args.space_policy)
    out = OUT / "real" / f"ZL_{args.space_policy}"
    out.mkdir(parents=True, exist_ok=True)
    point = {"file": FILE, "space_policy": args.space_policy, "source": "real", "n_tokens": len(ctx.real),
             **stats_for(ctx.real, 0, True, ctx), "prefixes": " ".join(ctx.prefixes), "suffixes": " ".join(ctx.suffixes)}
    pd.DataFrame([point]).to_csv(out / "point.csv", index=False)
    p003 = pd.read_csv(REF003 / "real" / f"ZL_{args.space_policy}" / "point.csv").iloc[0]
    common = [m for m in list(TARGETS) + EXTRA + NONTARGET if m in p003.index]
    diffs = {m: float(point[m]) - float(p003[m]) for m in common}
    (out / "check_vs_003.json").write_text(json.dumps(diffs, indent=2))
    print(pd.Series(point).to_string())
    print("max |diff| vs 003:", max(abs(v) for v in diffs.values()))


# --------------------------------------------------------------------------- #
# step 1: validation
# --------------------------------------------------------------------------- #

# The template of Rozanova & Temerev, rebuilt from their bundle copy of ZL3b with their token cleaning
# (voynich_decipherment_repro_bundle/decipherment_attack/plant_crib_attack.py, MIT): needed only to check
# that the port reproduces their published samples exactly.
_RT_LOCUS_RE = re.compile(r"^<([^>.]+)\.(\d+),([@+*=~&])([A-Za-z])[^>]*>\s*(.*)$")
_RT_BAD = set("?*<>{}[]()|@")


def _rt_clean_tokens(s: str) -> list[str]:
    s = re.sub(r"<[^>]*>", "", s)
    s = re.sub(r"\[([^:\]]*):[^\]]*\]", r"\1", s)
    s = re.sub(r"\{[^}]*\}", "", s)
    s = s.replace("'", "")
    out = []
    for w in re.split(r"[.,\s]+", s):
        w = w.strip().lower()
        if w and not any(c in _RT_BAD for c in w) and re.fullmatch(r"[a-z]+", w):
            out.append(w)
    return out


def rt_template(rt_repo: Path, min_tokens: int = 37500) -> list[tuple[str, int, bool, bool]]:
    """``voynich_template`` + ``extend_template(…, 37500)`` of their driver."""
    zl = rt_repo / "voynich_decipherment_repro_bundle" / "voynich_calibration_sources" / "ZL3b.txt"
    rows = []
    for raw in zl.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _RT_LOCUS_RE.match(raw)
        if not m or m.group(4) != "P":
            continue
        toks = _rt_clean_tokens(m.group(5))
        if toks:
            rows.append((m.group(1), len(toks), "<%>" in raw, "<$>" in raw))
    out, total, cycle = list(rows), sum(r[1] for r in rows), 0
    while total < min_tokens:
        cycle += 1
        for page, n, a, b in rows:
            out.append((f"{page}#{cycle}", n, a, b))
            total += n
            if total >= min_tokens:
                break
    return out


def table_from_lines(lines: list[list[str]], template: list[tuple], n_tokens: int) -> pd.DataFrame:
    """Token table of a plain-text sample: first ``n_tokens`` tokens; pages and paragraph flags from the
    generating template (other metadata synthetic, see :func:`voytools.nulls.as_table`)."""
    kept, kept_rows, count = [], [], 0
    for line, row in zip(lines, template):
        if count >= n_tokens:
            break
        take = line[: n_tokens - count]
        if take:
            kept.append(take)
            kept_rows.append(row)
            count += len(take)
    df = as_table(kept)
    lid = np.repeat(np.arange(len(kept)), [len(x) for x in kept])
    pg = np.asarray([r[0] for r in kept_rows], dtype=object)[lid]
    df["folio"] = pg
    df["locus_id"] = [f"{p}.{i}" for p, i in zip(pg, lid)]
    df["is_par_first_line"] = np.asarray([bool(r[2]) for r in kept_rows])[lid]
    df["is_par_last_line"] = np.asarray([bool(r[3]) for r in kept_rows])[lid]
    df.attrs["source"] = "selfcitation"
    return df


def cmd_validate(args) -> None:
    rt_repo = Path(args.rt_repo)
    tmpl = rt_template(rt_repo)
    ctx = Context("split")
    n_real = len(ctx.real)
    rows, ident = [], []
    for seed in (1, 2, 3):
        theirs = [ln.split() for ln in (rt_repo / "data" / "controls" / "selfcitation" / f"selfcite_seed{seed}.txt")
                  .read_text().splitlines()]
        ours = selfcitation_generate(tmpl, seed=seed, **{k: v for k, v in RT_PUBLISHED.items()})
        same = sum(a == b for a, b in zip(ours, theirs))
        ident.append({"seed": seed, "lines_theirs": len(theirs), "lines_port": len(ours), "lines_identical": same,
                      "tokens_theirs": sum(map(len, theirs)), "identical": int(same == len(theirs) == len(ours)),
                      "params": json.dumps({k: list(v) if isinstance(v, tuple) else v for k, v in RT_PUBLISHED.items()})})
        print(ident[-1], flush=True)
        rows.append({"source": "rt_samples", "seed": seed, "n_tokens": n_real,
                     **stats_for(table_from_lines(theirs, tmpl, n_real), seed, True)})
    defaults = config_id(JAVA_DEFAULT)
    for name, cid in (("port_rt_params", config_id(RT_PUBLISHED)), ("port_java_default", defaults)):
        for rep in range(5):
            seed = cid * 100 + rep  # = the grid seeds of this configuration
            rows.append({"source": name, "seed": seed, "n_tokens": n_real,
                         **stats_for(ctx.generate(cid, seed), seed, True)})
            print(name, rep, flush=True)
    rows.append({"source": "real_ZL_split", "seed": 0, "n_tokens": n_real, **stats_for(ctx.real, 0, True)})
    vdir = OUT / "validation"
    vdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ident).to_csv(vdir / "identity.csv", index=False)
    runs = pd.DataFrame(rows)
    runs.to_csv(vdir / "runs.csv", index=False, float_format="%.6g")
    write_validation(runs, pd.DataFrame(ident))


def write_validation(runs: pd.DataFrame, ident: pd.DataFrame) -> None:
    metrics = ["H2", "hapax_share", "bpe_gap64", "edge_mi_excess", "token_mi_excess_share", "near_ratio",
               "share_near", "n_types", "mean_len", "top50_coverage", "pos_jsd_first_mid"]
    ref = REFERENCE["ref_RT2026_selfcitation"]
    s = runs[runs["source"] == "rt_samples"]
    out = []
    for m in metrics:
        lo_s, hi_s, mean_s = s[m].min(), s[m].max(), s[m].mean()
        tab = ref.get(m, (np.nan, np.nan))
        tol = TOLERANCE.get(m, np.nan)
        for src, r in runs.groupby("source", sort=False):
            v = r[m].astype(float)
            row = {"metric": m, "source": src, "n": len(v), "mean": v.mean(), "min": v.min(), "max": v.max(),
                   "sd": v.std(ddof=1) if len(v) > 1 else np.nan, "table7_lo": tab[0], "table7_hi": tab[1],
                   "tolerance": tol, "samples_min": lo_s, "samples_max": hi_s}
            row["diff_vs_samples"] = v.mean() - mean_s
            row["within_samples_tol"] = (int(lo_s - tol <= v.mean() <= hi_s + tol) if not np.isnan(tol) else np.nan)
            row["diff_vs_table7"] = (v.mean() - (tab[0] + tab[1]) / 2) if not np.isnan(tab[0]) else np.nan
            row["within_table7_tol"] = (int(tab[0] - tol <= v.mean() <= tab[1] + tol)
                                        if not (np.isnan(tol) or np.isnan(tab[0])) else np.nan)
            out.append(row)
    val = pd.DataFrame(out)
    note = (f"port reproduces the published R&T samples token for token: "
            f"{int(ident['identical'].sum())}/{len(ident)} seeds identical "
            f"(p_copy 0, replace_weights 0/0/100, their template)")
    val.insert(0, "identity_check", note)
    val.to_csv(OUT / "validation.csv", index=False, float_format="%.6g")
    print(note)
    print(val.drop(columns=["identity_check"]).to_string(float_format=lambda x: f"{x:.4g}"))


# --------------------------------------------------------------------------- #
# step 2: calibration
# --------------------------------------------------------------------------- #


def grid_seeds() -> list[int]:
    return [cid * 100 + r for cid in range(len(CONFIGS)) for r in range(reps_of(cid))]


def cmd_grid(args) -> None:
    ctx = Context("split")
    keep = set(parse_seeds(args.seeds)) if args.seeds else None
    seeds = [s for s in grid_seeds() if keep is None or s in keep]

    def one_run(seed: int) -> dict:
        cid = seed // 100
        return {"space_policy": "split", "rep": seed % 100, **config_columns(cid),
                **stats_for(ctx.generate(cid, seed), seed, False)}

    run_seeds("004_grid", seeds, one_run, args.out or OUT / "grid" / "part0",
              params={"file": FILE, "space_policy": "split", "seeds": args.seeds,
                      "axes": {k: [list(x) if isinstance(x, tuple) else x for x in v] for k, v in AXES.items()},
                      "reps": {"p_same_position == 28": REPS_A, "other": REPS_B},
                      "metrics": "targets + R&T loss inputs + edge/token MI"},
              progress=progress)


def real_point(policy: str) -> pd.Series:
    return pd.read_csv(OUT / "real" / f"ZL_{policy}" / "point.csv").iloc[0]


def add_loss(means: pd.DataFrame, target: pd.Series) -> pd.DataFrame:
    """Loss B (001/003) and loss A (R&T) with their components."""
    out = means.copy()
    out["loss_B"] = 0.0
    out["n_missed_B"] = 0
    for f, w in TARGETS.items():
        comp = ((out[f] - target[f]) / target[f]) ** 2
        out[f"compB_{f}"] = comp
        out["loss_B"] += w * comp
        out["n_missed_B"] += (comp > COMPONENT_LIMIT).astype(int)
    out["loss_A"] = 0.0
    for f, scale in LOSS_A_SCALE.items():
        comp = ((out[f] - target[f]) / scale) ** 2
        out[f"compA_{f}"] = comp
        out["loss_A"] += comp
    return out


ID = ["config_id", "p_copy", "replace_weights", "p_same_position", "in_spec_grid"]


def config_means(runs: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    metrics = [m for m in metrics if m in runs.columns]
    g = runs.groupby(ID, sort=True)
    means = g[metrics].mean()
    means.insert(0, "n_seeds", g.size())
    return means.reset_index()


def cmd_rank(args) -> None:
    parts = sorted((OUT / "grid").glob("part*"))
    runs = merge_runs(parts, OUT / "grid")
    target = real_point("split")
    metrics = list(TARGETS) + EXTRA + ["top50_coverage", "edge_mi_excess", "token_mi_excess_share"]
    ra = runs[(runs["p_same_position"] == 28) & (runs["rep"] < REPS_A)]
    rb = runs[runs["rep"] < REPS_B]
    choice = {}
    for name, r, key in (("A", ra, "loss_A"), ("B", rb, "loss_B")):
        g = add_loss(config_means(r, metrics), target).sort_values(key, kind="stable").reset_index(drop=True)
        g.insert(0, "rank", np.arange(1, len(g) + 1))
        g.to_csv(OUT / f"grid_{name}.csv", index=False, float_format="%.6g")
        spec = g[g["in_spec_grid"] == 1]
        choice[f"cal{name}"] = int(spec.iloc[0]["config_id"])
        print(f"--- grid {name} (top 8)")
        print(g.head(8).to_string(float_format=lambda x: f"{x:.4g}"))
    choice["default"] = config_id(JAVA_DEFAULT)
    choice["rt_published"] = config_id(RT_PUBLISHED)
    choice["n_grid_runs"] = len(runs)
    (OUT / "phaseB_choice.json").write_text(json.dumps(choice, indent=2))
    print(choice)


# --------------------------------------------------------------------------- #
# step 2, phase B
# --------------------------------------------------------------------------- #

SOURCE_NAMES = {"default": "ts2020_default", "calA": "ts2020_calA", "calB": "ts2020_calB", "rt_published": "ts2020_RT"}


def phase_b_configs() -> dict[str, int]:
    c = json.loads((OUT / "phaseB_choice.json").read_text())
    return {SOURCE_NAMES[k]: int(c[k]) for k in SOURCE_NAMES}


def cmd_phase_b(args) -> None:
    ctx = Context(args.space_policy)
    cids = sorted(set(phase_b_configs().values()))
    seeds = [cid * 100 + r for cid in cids for r in range(args.reps)]
    if args.seeds:
        keep = set(parse_seeds(args.seeds))
        seeds = [s for s in seeds if s in keep]

    def one_run(seed: int) -> dict:
        cid = seed // 100
        return {"space_policy": args.space_policy, "rep": seed % 100, **config_columns(cid),
                **stats_for(ctx.generate(cid, seed), seed, True, ctx)}

    run_seeds("004_phaseB", seeds, one_run, args.out or OUT / "phaseB" / f"ZL_{args.space_policy}",
              params={"file": FILE, "space_policy": args.space_policy, "configs": phase_b_configs(),
                      "reps": args.reps, "prefixes": ctx.prefixes, "suffixes": ctx.suffixes},
              progress=progress)


def cmd_trigram86(args) -> None:
    """003 recommended null (trigram config 86) with this experiment's full metric set (same seeds as 003)."""
    ctx = Context(args.space_policy)
    p = pd.read_csv(REF003 / "best_configs.csv").iloc[0]
    assert int(p["config_id"]) == 86
    params = {"edit_model": p["edit_model"], "new_token_source": p["new_token_source"], "p_new": float(p["p_new"]),
              "window": int(p["window"]), "p_copy_exact": float(p["p_copy_exact"]),
              "p_two_edits": float(p["p_two_edits"]), "recency": float(p["recency"]), "edit_type_probs": (0.6, 0.2, 0.2)}
    allowed = ctx.real["token"].tolist()
    grammar = learn_word_grammar(allowed, ctx.real.attrs.get("glyph_preset"))

    def one_run(seed: int) -> dict:
        df = autocopy_like(ctx.real, allowed_from=allowed, word_grammar=grammar, seed=seed, **params)
        return {"space_policy": args.space_policy, "source": "trigram_86", **stats_for(df, seed, True, ctx)}

    seeds = [86 * 100 + r for r in range(args.reps)]  # 003 phase-B seeds of config 86
    run_seeds("004_trigram86", seeds, one_run, OUT / "trigram86" / f"ZL_{args.space_policy}",
              params={"config_003": 86, **{k: str(v) for k, v in params.items()}}, progress=progress)


# --------------------------------------------------------------------------- #
# step 3: vocabulary size
# --------------------------------------------------------------------------- #


def page_subsample(df: pd.DataFrame, frac: float, seed: int) -> pd.DataFrame:
    """Random subset of whole pages (without replacement, text order kept) with about ``frac`` of the tokens."""
    rng = np.random.default_rng(seed)
    folios = df["folio"].drop_duplicates().tolist()
    sizes = df.groupby("folio", sort=False).size()
    order = rng.permutation(len(folios))
    goal, total, keep = frac * len(df), 0, set()
    for i in order:
        if total >= goal:
            break
        keep.add(folios[i])
        total += int(sizes[folios[i]])
    out = df[df["folio"].isin(keep)].reset_index(drop=True)
    out.attrs.update(df.attrs)
    return out


def scaled_template(df: pd.DataFrame, scale: float) -> pd.DataFrame:
    """The real table cut to its first pages (``scale`` < 1) or extended by cycling its pages (> 1)."""
    goal = scale * len(df)
    parts, total, cycle = [], 0, 0
    while total < goal:
        for folio, g in df.groupby("folio", sort=False):
            g = g.copy()
            if cycle:
                g["folio"] = f"{folio}#{cycle}"
                g["locus_id"] = g["locus_id"].astype(str) + f"#{cycle}"
            parts.append(g)
            total += len(g)
            if total >= goal:
                break
        cycle += 1
    out = pd.concat(parts, ignore_index=True)
    out["idx"] = np.arange(len(out))
    out.attrs.update(df.attrs)
    return out


def _bisect(f, target: float, lo: float, hi: float, steps: int = 7) -> float:
    """``x`` in [lo, hi] with f(x) ~ target for an increasing f."""
    for _ in range(steps):
        mid = (lo + hi) / 2
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def cmd_vocab(args) -> None:
    ctx = Context("split")
    cid = phase_b_configs()[args.source]
    runs_b = load_runs(OUT / "phaseB" / "ZL_split")
    gen = runs_b[runs_b["config_id"] == cid]
    real = real_point("split")
    n_gen, n_real = float(gen["n_stems_k"].mean()), float(real["n_stems_k"])
    print(f"n_stems_k: real {n_real:.0f}, generator {n_gen:.1f} (config {cid})", flush=True)
    rows = []

    def add(cond: str, rep: int, scale: float, df: pd.DataFrame, seed: int) -> None:
        rows.append({"condition": cond, "rep": rep, "scale": scale, "n_tokens": len(df),
                     **paradigm_block(df, ctx.prefixes, ctx.suffixes, seed)})

    # (a) real text subsampled by pages to the generator's n_stems_k (only possible if the generator has fewer)
    if n_gen < n_real:
        def real_k(frac):
            return np.mean([paradigm_block(page_subsample(ctx.real, frac, 1000 + i), ctx.prefixes, ctx.suffixes, 0,
                                           sections=False)["n_stems_k"] for i in range(5)])
        frac = _bisect(real_k, n_gen, 0.05, 1.0)
        for i in range(args.reps):
            add("real_subsampled", i, frac, page_subsample(ctx.real, frac, i), i)
    # (b) generator on a template scaled so that its n_stems_k matches the real text
    def gen_k(scale):
        t = scaled_template(ctx.real, scale)
        return np.mean([paradigm_block(selfcitation_like(t, seed=cid * 100 + 90 + i, **config_params(cid)),
                                       ctx.prefixes, ctx.suffixes, 0, sections=False)["n_stems_k"] for i in range(3)])
    lo, hi = (0.1, 1.0) if n_gen > n_real else (1.0, 6.0)
    scale = _bisect(gen_k, n_real, lo, hi)
    t = scaled_template(ctx.real, scale)
    for i in range(args.reps):
        seed = cid * 100 + i
        add("generator_scaled", i, scale, selfcitation_like(t, seed=seed, **config_params(cid)), seed)
    vdir = OUT / "vocab"
    vdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(vdir / f"runs_{args.source}.csv", index=False, float_format="%.6g")
    write_vocab()


def write_vocab() -> None:
    """vocab_matching.csv: real vs generator, before and after matching n_stems_k."""
    real = real_point("split")
    runs_b = load_runs(OUT / "phaseB" / "ZL_split")
    boot = load_runs(REF002 / "parts" / "ZL_split_learned" / "real_boot")
    out = []
    for path in sorted((OUT / "vocab").glob("runs_*.csv")):
        source = path.stem[len("runs_"):]
        cid = phase_b_configs()[source]
        gen = runs_b[runs_b["config_id"] == cid]
        v = pd.read_csv(path)
        sub, scl = v[v["condition"] == "real_subsampled"], v[v["condition"] == "generator_scaled"]
        for m in PARADIGM:
            x = float(real[m])
            g = gen[m].astype(float)
            row = {"generator": source, "config_id": cid, "metric": m, "real_full": x,
                   "real_full_sd_boot": float(boot[m].std(ddof=1)) if m in boot else np.nan,
                   "gen_full": g.mean(), "gen_full_sd": g.std(ddof=1),
                   "diff_full": x - g.mean()}
            if len(sub):
                row.update({"real_sub_scale": float(sub["scale"].iloc[0]), "real_sub": sub[m].mean(),
                            "real_sub_sd": sub[m].std(ddof=1), "diff_matched_a": sub[m].mean() - g.mean()})
            else:
                row.update({"real_sub_scale": np.nan, "real_sub": np.nan, "real_sub_sd": np.nan,
                            "diff_matched_a": np.nan})
            row.update({"gen_scaled_scale": float(scl["scale"].iloc[0]), "gen_scaled": scl[m].mean(),
                        "gen_scaled_sd": scl[m].std(ddof=1), "diff_matched_b": x - scl[m].mean()})
            for k in ("a", "b"):
                d = row[f"diff_matched_{k}"]
                row[f"diff_change_{k}"] = (abs(d) - abs(row["diff_full"])) / abs(row["diff_full"]) \
                    if row["diff_full"] and not np.isnan(d) else np.nan
            out.append(row)
    pd.DataFrame(out).to_csv(OUT / "vocab_matching.csv", index=False, float_format="%.6g")
    print(pd.DataFrame(out).to_string(float_format=lambda x: f"{x:.4g}"))


# --------------------------------------------------------------------------- #
# step 4: summaries
# --------------------------------------------------------------------------- #


def _source_row(row: dict, name: str, r: pd.DataFrame, cid, m: str, x: float, note: str = "") -> dict:
    if m not in r.columns or r[m].isna().all():
        return {**row, "source": name, "config_id": cid, "n": 0, "mean": np.nan, "lo": np.nan, "hi": np.nan,
                "sd": np.nan, "interval": "not computed", "diff_real": np.nan, "z_real": np.nan}
    v = r[m].dropna().astype(float)
    sd = v.std(ddof=1)
    return {**row, "source": name, "config_id": cid, "n": len(v), "mean": v.mean(),
            "lo": v.quantile(0.025), "hi": v.quantile(0.975), "sd": sd,
            "interval": "2.5/97.5 % over seeds" + (f"; {note}" if note else ""),
            "diff_real": x - v.mean(), "z_real": (x - v.mean()) / sd if sd > 1e-9 * max(abs(x), 1) else np.nan}


def _reference_sources(policy: str) -> list[tuple[str, pd.DataFrame, str, str]]:
    """trigram_86 (003, completed here), bigram_573 (001 + 002), within-line shuffle (001 + 002)."""
    tag = f"ZL_{policy}"
    out = []
    t86 = OUT / "trigram86" / tag
    r = load_runs(t86) if (t86 / "runs.csv").exists() else load_runs(REF003 / "phaseB" / tag).query("config_id == 86")
    out.append(("trigram_86", r, "003:86", "seeds 8600..8649 as in 003"))
    b = load_runs(REF001 / "phaseB" / tag)
    b = b[b["config_id"] == 573].reset_index(drop=True)
    p2 = load_runs(REF002 / "parts" / f"{tag}_learned" / "autocopy_573")
    for m in PARADIGM + ["jsd_suffix_AB"]:  # paradigm metrics of other samples of the same generator (002)
        if m in p2.columns and (m not in b.columns or b[m].isna().all()):
            b[m] = pd.Series(p2[m].to_numpy()[: len(b)]) if len(p2) >= len(b) else np.nan
    out.append(("bigram_573", b, "001:573", "metrics of 002 on other seeds of the same config"))
    s = load_runs(REF001 / "shuffle" / tag)
    s2 = load_runs(REF002 / "parts" / f"{tag}_learned" / "shuffle_lines")
    for m in PARADIGM + ["jsd_suffix_AB"]:
        if m in s2.columns and m not in s.columns:
            s[m] = pd.Series(s2[m].to_numpy()[: len(s)])
    out.append(("shuffle_lines", s, "", "001 (order metrics) + 002 (paradigms)"))
    return out


def real_sd(policy: str) -> dict[str, float]:
    tag = f"ZL_{policy}"
    sd = {}
    boot1 = load_runs(REF001 / "real" / tag / "bootstrap")
    for m in boot1.columns:
        if m not in ("seed", "elapsed_s"):
            sd[m] = float(boot1[m].std(ddof=1))
    boot2 = load_runs(REF002 / "parts" / f"{tag}_learned" / "real_boot")
    for m in PARADIGM + ["jsd_suffix_AB"]:
        if m in boot2.columns and m not in sd:
            sd[m] = float(boot2[m].std(ddof=1))
    return sd


def kind_of(m: str) -> str:
    return "target" if m in TARGETS else "extra" if m in EXTRA else "nontarget" if m in NONTARGET else "check"


def cmd_summarize(args) -> None:
    configs = phase_b_configs()
    grid_a, grid_b = pd.read_csv(OUT / "grid_A.csv"), pd.read_csv(OUT / "grid_B.csv")
    metrics = list(TARGETS) + EXTRA + NONTARGET + CHECKS
    info: dict = {"configs": configs}
    profile_rows, summary_rows, cand = [], [], []
    runs_all = []
    for policy in ("split", "join"):
        tag = f"ZL_{policy}"
        runs_b = load_runs(OUT / "phaseB" / tag)
        runs_all.append(runs_b)
        x_real = real_point(policy)
        sds = real_sd(policy)
        info[f"phaseB_{tag}_runs"] = len(runs_b)
        info[f"phaseB_{tag}_suspicious_runs"] = int(((runs_b["share_near"] > 0.9) | (runs_b["H2"] < 2.0)).sum())
        sources = [(name, runs_b[runs_b["config_id"] == cid], cid, "") for name, cid in configs.items()]
        sources += _reference_sources(policy)
        for m in metrics:
            x = float(x_real[m]) if m in x_real and pd.notna(x_real[m]) else np.nan
            row = {"dataset": tag, "metric": m, "kind": kind_of(m)}
            sd = sds.get(m, np.nan)
            profile_rows.append({**row, "source": "real", "config_id": "", "n": 1, "mean": x,
                                 "lo": x - 1.96 * sd if not np.isnan(sd) else np.nan,
                                 "hi": x + 1.96 * sd if not np.isnan(sd) else np.nan, "sd": sd,
                                 "interval": "page bootstrap of 001/002, normal (x +- 1.96 SD)" if not np.isnan(sd) else "",
                                 "z_real": np.nan, "diff_real": np.nan})
            for name, r, cid, note in sources:
                profile_rows.append(_source_row(row, name, r, cid, m, x, note))
            if policy == "split":
                for ref, vals in REFERENCE.items():
                    if m in vals:
                        lo, hi = vals[m]
                        profile_rows.append({**row, "source": ref, "config_id": "", "n": np.nan, "mean": (lo + hi) / 2,
                                             "lo": lo, "hi": hi, "sd": np.nan,
                                             "interval": "docs/reference_values.md (range if given)",
                                             "z_real": np.nan, "diff_real": x - (lo + hi) / 2})
        for name, cid in configs.items():
            r = runs_b[runs_b["config_id"] == cid]
            summ = r[[m for m in metrics if m in r]].describe(percentiles=[0.025, 0.5, 0.975]).T
            summ.insert(0, "dataset", tag)
            summ.insert(1, "source", name)
            summ.insert(2, "config_id", cid)
            summary_rows.append(summ.reset_index().rename(columns={"index": "metric"}))
            if policy == "split":
                means = add_loss(config_means(r, list(TARGETS) + EXTRA + NONTARGET + CHECKS), x_real)
                means.insert(0, "source", name)
                cand.append(means)
    profile = pd.DataFrame(profile_rows)
    profile.to_csv(OUT / "profile.csv", index=False, float_format="%.6g")
    pd.concat(summary_rows, ignore_index=True).to_csv(OUT / "summary.csv", index=False, float_format="%.6g")
    merge = pd.concat(runs_all, ignore_index=True)
    merge.to_csv(OUT / "runs.csv", index=False, float_format="%.6g")
    write_best(pd.concat(cand, ignore_index=True), grid_a, grid_b, info)
    info["lists"] = three_lists(profile)
    (OUT / "summary_info.json").write_text(json.dumps(info, indent=2, default=str))
    plot_profile(profile)
    print(json.dumps(info, indent=2, default=str))


def write_best(cand: pd.DataFrame, grid_a: pd.DataFrame, grid_b: pd.DataFrame, info: dict) -> None:
    """best_configs.csv: calB first if it meets the criterion, else the better of A/B (loss B), provisional."""
    comps = [f"compB_{f}" for f in TARGETS]
    cand = cand.copy()
    cand["meets_criterion"] = ((cand["loss_B"] < LOSS_LIMIT) & (cand[comps] < COMPONENT_LIMIT).all(axis=1)).astype(int)
    cand["missed_targets"] = cand.apply(lambda r: " ".join(f for f in TARGETS if r[f"compB_{f}"] >= COMPONENT_LIMIT), axis=1)
    cal = cand[cand["source"].isin(["ts2020_calA", "ts2020_calB"])].drop_duplicates("source")
    calb = cal[cal["source"] == "ts2020_calB"].iloc[0]
    if bool(calb["meets_criterion"]):
        first, prov = "ts2020_calB", 0
    else:
        first, prov = cal.sort_values("loss_B", kind="stable").iloc[0]["source"], 1
    order = [first] + [s for s in ("ts2020_calA", "ts2020_calB", "ts2020_RT", "ts2020_default") if s != first]
    best = cand.set_index("source").loc[[s for s in order if s in set(cand["source"])]].reset_index()
    best.insert(1, "rank", np.arange(1, len(best) + 1))
    best["recommended"] = (best["source"] == first).astype(int)
    best["provisional"] = np.where(best["source"] == first, prov, 1)
    best["note"] = [
        ("recommended null; " if s == first else "")
        + ("meets criterion (loss_B < 0.15, all components < 0.05)" if ok else f"misses: {miss}")
        + ("; R&T published parameters, outside the spec grid" if s == "ts2020_RT" else "")
        + ("; Java/API defaults" if s == "ts2020_default" else "")
        for s, ok, miss in zip(best["source"], best["meets_criterion"], best["missed_targets"])]
    for name, g, key in (("A", grid_a, "loss_A"), ("B", grid_b, "loss_B")):
        best[f"grid{name}_rank"] = best["config_id"].map(g.set_index("config_id")["rank"])
        best[f"grid{name}_{key}"] = best["config_id"].map(g.set_index("config_id")[key])
    lead = ["source", "rank", "config_id", "p_copy", "replace_weights", "p_same_position", "in_spec_grid", "recommended",
            "provisional", "note", "loss_B", "loss_A", "n_missed_B", "n_seeds"]
    best = best[lead + [c for c in best.columns if c not in lead]]
    best.to_csv(OUT / "best_configs.csv", index=False, float_format="%.6g")
    info["recommended"] = first
    info["recommended_provisional"] = prov
    print(best.head(4).to_string(float_format=lambda x: f"{x:.4g}"))


def three_lists(profile: pd.DataFrame, z_limit: float = 3.0) -> dict:
    """(a) reproduced in split and join, (b) not reproduced in both, (c) inconsistent - per generator source."""
    out = {}
    nt = profile[(profile["kind"] == "nontarget")]
    for src in ("ts2020_calB", "ts2020_calA", "ts2020_default", "ts2020_RT", "trigram_86"):
        sub = nt[nt["source"] == src]
        lists = {"reproduced": [], "not_reproduced": [], "inconsistent": []}
        for m, g in sub.groupby("metric", sort=False):
            z = g.set_index("dataset")["z_real"].abs()
            ok = [bool(v < z_limit) for v in z.reindex(["ZL_split", "ZL_join"])]
            if any(np.isnan(v) for v in z.reindex(["ZL_split", "ZL_join"])):
                continue
            key = "reproduced" if all(ok) else "not_reproduced" if not any(ok) else "inconsistent"
            lists[key].append(m)
        out[src] = lists
    return out


def plot_profile(profile: pd.DataFrame) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; profile.png skipped")
        return
    profile = profile[profile["kind"] != "check"]
    metrics = list(dict.fromkeys(profile["metric"]))
    datasets = list(dict.fromkeys(profile["dataset"]))
    colors = {"real": "#222222", "ts2020_calB": "#c0392b", "ts2020_calA": "#e67e22", "ts2020_default": "#8e44ad",
              "ts2020_RT": "#d35400", "trigram_86": "#2c8d5a", "bigram_573": "#8a8a8a", "shuffle_lines": "#3a78b5",
              "ref_RT2026_voynich": "#000000", "ref_RT2026_selfcitation": "#c9302c"}
    markers = {"ref_RT2026_voynich": "*", "ref_RT2026_selfcitation": "s"}
    ncol = 4
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.6 * nrow))
    for ax, m in zip(axes.flat, metrics):
        sub = profile[profile["metric"] == m]
        kind = sub["kind"].iloc[0]
        for di, ds in enumerate(datasets):
            d = sub[(sub["dataset"] == ds) & sub["mean"].notna()]
            for si, (_, r) in enumerate(d.iterrows()):
                y = di + 0.08 * si - 0.4
                lo, hi = r["lo"], r["hi"]
                xerr = None if pd.isna(lo) else [[max(r["mean"] - lo, 0)], [max(hi - r["mean"], 0)]]
                ax.errorbar(r["mean"], y, xerr=xerr, fmt=markers.get(r["source"], "o"), ms=3.5,
                            color=colors.get(r["source"], "grey"), label=r["source"], capsize=2)
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets, fontsize=8)
        ax.set_title(f"{m} ({kind})", fontsize=9, color="#000000" if kind == "target" else "#7a2a00")
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.3)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    handles, labels = [], []
    for ax in axes.flat:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    fig.legend(handles, labels, loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("004: real ZL text vs ported Timm & Schinner generator (default, calA, calB, R&T parameters), "
                 "003 trigram 86, 001 bigram 573, within-line shuffle, literature (95 % intervals)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "profile.png", dpi=100)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["real", "validate", "grid", "rank", "phaseB", "trigram86", "vocab", "summarize"])
    ap.add_argument("--space-policy", default="split", choices=["split", "join"])
    ap.add_argument("--rt-repo", default="/tmp/voynich-units", help="clone of github.com/lrozanova/voynich-units")
    ap.add_argument("--seeds", default=None, help="encoded seeds config_id*100+rep: 'a:b', 'a:b/chunks:k' or list")
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--source", default="ts2020_calB", help="vocab: which phase-B generator to match")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.command == "phaseB":
        args.reps = args.reps or REPS_PHASE_B
        cmd_phase_b(args)
    elif args.command == "trigram86":
        args.reps = args.reps or REPS_PHASE_B
        cmd_trigram86(args)
    elif args.command == "vocab":
        args.reps = args.reps or 20
        cmd_vocab(args)
    else:
        {"real": cmd_real, "validate": cmd_validate, "grid": cmd_grid, "rank": cmd_rank,
         "summarize": cmd_summarize}[args.command](args)


if __name__ == "__main__":
    main()
