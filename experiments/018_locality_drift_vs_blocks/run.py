"""Experiment 018: locality of the page vocabulary - smooth drift versus blocks, and physical order.

Usage (from the repository root; see ``spec.md``)::

    python experiments/018_locality_drift_vs_blocks/run.py run --combo ZL_split                   # all sources
    python experiments/018_locality_drift_vs_blocks/run.py run --combo GC_join --sources real,perm_global
    python experiments/018_locality_drift_vs_blocks/run.py run --combo ctl:pliny_prose             # control (009)
    python experiments/018_locality_drift_vs_blocks/run.py memory --phase grid                     # ts2020-mem grid
    python experiments/018_locality_drift_vs_blocks/run.py memory --phase winner                   # 30 seeds
    python experiments/018_locality_drift_vs_blocks/run.py extras --combo ZL_split                 # matrices, lists
    python experiments/018_locality_drift_vs_blocks/run.py summarize

Combinations: ``{ZL, GC}_{split, join}`` (real paragraph text, clean) and ``ctl:<name>`` for the 009
controls (``pliny_prose``, ``celsus_prose``, ``linnaeus_items``; ZL split layout, sections = books /
orders) and ``ctl:gibberish_pooled`` (013; 60 pages, points 1-2 only).

Sources (``seed`` = repetition):

* ``real`` - the real text; every seed is one random split of the stems into halves (cross-validated
  block gain), 50 seeds;
* ``perm_global`` / ``perm_section`` - whole pages permuted among page positions globally / inside the
  section (page content moves, position metadata stay), 200 seeds each; 5 splits averaged per seed;
* ``selfcitation_21`` - ts2020 config 21 (004), seeds 2100 + rep, 30 reps, Eva only;
* ``selfcitation_mem`` - ts2020 config 21 with memory across pages (``p_prev_page`` / ``page_memory``
  from ``memory_calibration.csv``, phase ``winner``), seeds 2300 + rep, 30 reps, Eva only;
* ``template_005`` - the 005 template with topic (``best`` configuration), seeds 100000 + rep, 30 reps.

Every run reports, per measure (``rare`` = rare stem sharing 2-5, ``cos``, ``jac``) and cut (``all`` =
whole manuscript, ``within`` = pairs inside sections only, ``herbalA`` = herbal pages of language A),
the pooled similarity profile per distance (``prof_d<d>``), the cross-validated block gain curve
(``cv_K<K>``), the in-sample gain (``fit_K<K>``), the change-point precision/recall against section /
language / hand / quire boundaries at ``K = K_ref, 2 K_ref, 3 K_ref`` with random baselines
(``cp<m>x_<ref>_prec`` ...), and (cut ``all``) the physical-order metrics: distance-matched excess of
same-bifolio / same-quire pairs, nearest-neighbour distance shares and quire-order agreement.
Column names: ``<measure>__<cut>__<metric>``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_runs, load_tables, run_seeds  # noqa: E402
from voytools import locality as L  # noqa: E402
from voytools import recurrence as R  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.selfcitation import selfcitation_like  # noqa: E402
from voytools.templates import ParadigmTemplate  # noqa: E402


def _load(name: str, rel: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M6 = _load("exp006", "006_topic_at_distance")

DATA = ROOT / "data"
OUT = ROOT / "results" / "018_locality_drift_vs_blocks"
FILES = {"ZL": "ZL3b-n.txt", "GC": "GC2a-n.txt"}
COMBOS = ["ZL_split", "ZL_join", "GC_split", "GC_join"]
CONTROLS_009 = ("pliny_prose", "celsus_prose", "linnaeus_items")
CONTROLS_013 = ("gibberish_pooled",)
MEASURES = L.MEASURES
CUTS = ("all", "within", "herbalA")
K_MAX = 30
DMAX = L.DMAX
REFS = ("section", "language", "hand", "quire")
CP_MULT = (1, 2, 3)
N_SPLITS_NULL = 5
MIN_CUT_PAGES = 12
SOURCES = {"real": range(50), "perm_global": range(200), "perm_section": range(200),
           "selfcitation_21": range(30), "selfcitation_mem": range(30), "template_005": range(30)}
CONTROL_SOURCES = ("real", "perm_global", "perm_section")
GIBBERISH_SOURCES = ("real", "perm_global")
TS_SEED0, TSM_SEED0, TPL_SEED0 = 2100, 2300, 100_000
MEM_GRID = {"p_prev_page": (0.1, 0.2, 0.3, 0.5), "page_memory": (2, 4, 8, 16)}
MEM_GRID_SEEDS, MEM_WINNER_SEEDS = 5, 30
MEM_GRID_SEED0 = 2400
CAL_BINS = ("1", "2-3", "4-7", "8-15")


# --------------------------------------------------------------------------- #
# corpora
# --------------------------------------------------------------------------- #


class Context:
    """One corpus (real combination or control): token table, page table, affix lists."""

    def __init__(self, combo: str):
        self.combo = combo
        self.is_control = combo.startswith("ctl:")
        if self.is_control:
            name = combo[4:]
            self.key, self.policy = "ZL", "split"
            if name in CONTROLS_013:
                M13 = _load("exp013", "013_gibberish_and_cipher_controls")
                self.real = M13.build(name).reset_index(drop=True)
            else:
                M9 = _load("exp009", "009_register_profile")
                self.real = M9.control_table(name).reset_index(drop=True)
            self.file = name
            self.rank = R.page_order(self.real)
        else:
            self.key, self.policy = combo.split("_")
            self.file = FILES[self.key]
            tables = load_tables(DATA / self.file, space_policy=self.policy)
            self.rank = R.page_order(tables.tokens)
            self.real = filter_tokens(tables.tokens, generic_type="P", clean=True).reset_index(drop=True)
        self.real_d = decompose_table(self.real)
        self.prefixes, self.suffixes = self.real_d.attrs["prefixes"], self.real_d.attrs["suffixes"]
        self.pages = R.page_table(self.real, self.rank)
        self.folios = self.pages.index.to_numpy()
        self.dist = L.distance_matrix(self.pages["rank"].to_numpy())
        self.cuts = self.make_cuts()
        self._tpl = None

    def columns(self) -> dict:
        return {"combo": self.combo, "file": self.file, "space_policy": self.policy, "n_pages": len(self.pages)}

    def decompose(self, df: pd.DataFrame) -> pd.DataFrame:
        return decompose_table(df, prefixes=self.prefixes, suffixes=self.suffixes)

    def counts(self, df: pd.DataFrame) -> np.ndarray:
        return L.stem_page_counts(self.decompose(df), self.folios)[0]

    def make_cuts(self) -> dict:
        """``{cut: (page index array, same-section-mask-or-None, reference label for K_ref)}``."""
        p = self.pages
        n = len(p)
        sec = p["section"].astype(str).to_numpy()
        cuts = {"all": (np.arange(n), None, "section"),
                "within": (np.arange(n), sec[:, None] == sec[None, :], "section")}
        if not self.is_control:
            idx = np.flatnonzero((sec == "H") & (p["language"].astype(str).to_numpy() == "A"))
            if len(idx) >= MIN_CUT_PAGES:
                cuts["herbalA"] = (idx, None, "quire")
        return cuts

    def template(self) -> ParadigmTemplate:
        if self._tpl is None:
            self._tpl = ParadigmTemplate.learn(self.real_d, **M6.template_config(self.key))
        return self._tpl


# --------------------------------------------------------------------------- #
# metrics of one page x stem count matrix
# --------------------------------------------------------------------------- #


def _r(v) -> float:
    return round(float(v), 6) if np.isfinite(v) else float("nan")


def cut_metrics(measure: str, counts: np.ndarray, weight: np.ndarray, dist: np.ndarray, mask, labels: pd.DataFrame,
                ref: str, split_seeds: list[int], physical: bool) -> dict:
    """All metrics of one measure on one cut (pages already restricted, ``weight`` renormalised)."""
    out = {}
    n_mat, w_mat = L.similarity_matrices(measure, counts, weight)
    prof = L.pooled_profile(n_mat, w_mat, dist, mask, DMAX)
    for d, ratio in zip(prof["d"], prof["ratio"]):
        out[f"prof_d{int(d)}"] = _r(ratio)
    r_full, w_full, f = L.residuals(n_mat, w_mat, dist, mask, DMAX)
    out["n_pairs_w"] = _r(w_full[np.triu_indices(len(w_full), 1)].sum())
    # cross-validated and in-sample gain curves
    cv = np.nanmean([L.crossval_gain(measure, counts, weight, dist, mask, s, K_MAX, DMAX, f) for s in split_seeds], axis=0)
    w_model = f[dist] * w_full if measure == "rare" else w_full
    num, den, kind = L.segmentation_inputs(r_full, w_model, measure)
    segs = L.segment(num, den, K_MAX, kind)
    for k in range(1, K_MAX + 1):
        out[f"cv_K{k}"] = _r(cv[k - 1])
        out[f"fit_K{k}"] = _r(L.block_gain(r_full, w_full, segs[k])) if k in segs else float("nan")
    # change points
    n = len(labels)
    refs = {c: L.boundaries_of(labels[c].astype(str).to_numpy()) for c in REFS}
    k_ref = len(refs[ref]) + 1
    out["K_ref"] = k_ref
    for m in CP_MULT:
        k = min(max(m * k_ref, 2), K_MAX)
        found = segs[k]
        for c in REFS:
            prec, rec = L.match_boundaries(found, refs[c], 1)
            prec0, rec0 = L.random_boundary_baseline(n, len(found), refs[c], 1, 200, seed=k)
            out[f"cp{m}x_{c}_prec"], out[f"cp{m}x_{c}_rec"] = _r(prec), _r(rec)
            out[f"cp{m}x_{c}_prec0"], out[f"cp{m}x_{c}_rec0"] = _r(prec0), _r(rec0)
            out[f"cp{m}x_{c}_n_ref"] = int(len(refs[c]))
        out[f"cp{m}x_K"] = int(k)
    if not physical:
        return out
    # physical order (distance-matched excess, residual in units of the pooled similarity)
    for c, dmax_phys in (("bifolio", 15), ("quire", 15)):
        v = labels[c].astype(str).to_numpy()
        same = v[:, None] == v[None, :]
        res = L.matched_excess(r_full, w_model, dist, same, dmax_phys)
        for k, val in res.items():
            out[f"phys_{c}_{k}"] = val if k == "n_same" else _r(val)
    # neighbours (same language), score = standardised residual for rare, residual for cos / jac
    lang = labels["language"].astype(str).to_numpy()
    allowed = (lang[:, None] == lang[None, :]) & (w_full > 0)
    if measure == "rare":
        with np.errstate(invalid="ignore", divide="ignore"):
            score = np.where(w_model > 0, r_full / np.sqrt(np.where(w_model > 0, w_model, 1.0)), -np.inf)
    else:
        score = np.where(w_full > 0, r_full, -np.inf)
    nb = L.neighbours(score, dist, allowed, 3, L.NEAR_EDGES, 8)
    for k in ("n_pages", "share_le2", "share_le5", "share_gt16", "n_anomalous"):
        out[f"nb_{k}"] = nb[k] if k.startswith("n_") else _r(nb[k])
    q = L.quire_matrix(r_full, w_model, labels["quire"].astype(str).to_numpy())
    for k, val in L.quire_order_metrics(q).items():
        out[f"qo_{k}"] = _r(val)
    return out


def all_metrics(ctx: Context, counts: np.ndarray, weight: np.ndarray, split_seeds: list[int],
                cuts: tuple = CUTS) -> dict:
    out = {}
    for cut in cuts:
        if cut not in ctx.cuts:
            continue
        idx, mask, ref = ctx.cuts[cut]
        sub = counts[idx]
        w = weight[idx] / weight[idx].sum()
        dist = ctx.dist[np.ix_(idx, idx)]
        labels = ctx.pages.iloc[idx]
        for measure in MEASURES:
            res = cut_metrics(measure, sub, w, dist, mask, labels, ref, split_seeds, physical=(cut == "all"))
            for k, v in res.items():
                out[f"{measure}__{cut}__{k}"] = v
    return out


def page_permutation(ctx: Context, within: str | None, seed: int) -> np.ndarray:
    """``perm[p]`` = index of the page whose content goes to position ``p`` (inside groups ``within``)."""
    rng = np.random.default_rng(seed)
    n = len(ctx.pages)
    if within is None:
        return rng.permutation(n)
    groups = ctx.pages[within].astype(str).to_numpy()
    perm = np.arange(n)
    for g in np.unique(groups):
        i = np.flatnonzero(groups == g)
        perm[i] = i[rng.permutation(len(i))]
    return perm


def memory_params() -> dict:
    """Winner of the ts2020-mem calibration (``memory_calibration.csv``, column ``winner``)."""
    cal = pd.read_csv(OUT / "memory_calibration.csv")
    row = cal[cal["winner"] == 1].iloc[0]
    return {"p_prev_page": float(row["p_prev_page"]), "page_memory": float(row["page_memory"])}


def make_runner(ctx: Context, source: str, cuts: tuple):
    weight = ctx.pages["weight"].to_numpy()
    real_counts = ctx.counts(ctx.real)

    def one(seed: int) -> dict:
        base = {**ctx.columns(), "source": source}
        null_splits = [seed * 10 + s for s in range(N_SPLITS_NULL)]
        if source == "real":
            return {**base, **all_metrics(ctx, real_counts, weight, [seed], cuts)}
        if source.startswith("perm_"):
            perm = page_permutation(ctx, {"perm_global": None, "perm_section": "section"}[source], seed)
            return {**base, **all_metrics(ctx, real_counts[perm], weight[perm], null_splits, cuts)}
        if source in ("selfcitation_21", "selfcitation_mem"):
            if ctx.key != "ZL" or ctx.is_control:
                raise ValueError("selfcitation_like generates basic Eva; no null for a non-Eva file")
            params = M6.selfcitation_params()
            if source == "selfcitation_mem":
                sample = selfcitation_like(ctx.real, seed=TSM_SEED0 + seed, **params, **memory_params())
            else:
                sample = selfcitation_like(ctx.real, seed=TS_SEED0 + seed, **params)
        elif source == "template_005":
            sample = ctx.template().sample(ctx.real, seed=TPL_SEED0 + seed, pool_seed=seed)
        else:
            raise ValueError(source)
        return {**base, **all_metrics(ctx, ctx.counts(sample), weight, null_splits, cuts)}

    return one


def progress(i: int, n: int, s: int) -> None:
    if i % 10 == 0 or i == n:
        print(f"  {i}/{n} seed={s}", flush=True)


def cmd_run(args) -> None:
    if args.combo.startswith("ctl:"):
        default = GIBBERISH_SOURCES if args.combo[4:] in CONTROLS_013 else CONTROL_SOURCES
    else:
        default = list(SOURCES)
    sources = list(default) if args.sources == "all" else args.sources.split(",")
    ctx = Context(args.combo)
    cuts = ("all", "within") if args.combo[4:] in CONTROLS_013 else CUTS
    parts = OUT / "runs" / args.combo.replace(":", "_")
    parts.mkdir(parents=True, exist_ok=True)
    ctx.pages.reset_index().to_csv(parts / "pages.csv", index=False)
    for source in sources:
        seeds = M6.parse_seeds(args.seeds) if args.seeds else list(SOURCES[source])
        print(f"{args.combo} {source}: {len(seeds)} seeds", flush=True)
        params = {"combo": args.combo, "source": source, "measures": list(MEASURES), "cuts": list(cuts),
                  "k_max": K_MAX, "dmax": DMAX, "n_splits_null": N_SPLITS_NULL, "n_pages": len(ctx.pages),
                  "cut_pages": {c: int(len(ctx.cuts[c][0])) for c in cuts if c in ctx.cuts}}
        if source == "selfcitation_21":
            params["selfcitation_004"] = M6.selfcitation_params()
            params["seed_offset"] = TS_SEED0
        if source == "selfcitation_mem":
            params["selfcitation_004"] = M6.selfcitation_params()
            params["memory"] = memory_params()
            params["seed_offset"] = TSM_SEED0
        if source == "template_005":
            params["template_005"] = M6.template_config(ctx.key)
            params["seed_offset"], params["pool_seed"] = TPL_SEED0, "repetition"
        if ctx.is_control:
            params["control"] = {k: v for k, v in ctx.real.attrs.get("control", {}).items()}
        run_seeds(f"018_{args.combo}_{source}", seeds, make_runner(ctx, source, cuts), parts / source,
                  params=params, progress=progress)


# --------------------------------------------------------------------------- #
# ts2020 with memory across pages: calibration on the 015 bins
# --------------------------------------------------------------------------- #


def recurrence_ratios(ctx: Context, df: pd.DataFrame) -> dict:
    """015 profile of a sample on the pages of ``ctx``: observed/expected per bin, far-section excess."""
    occ = R.stem_occurrences(ctx.decompose(df), ctx.pages, 2, 5)
    m = R.recurrence_metrics(occ, ctx.pages)
    out = {"n_stems": m["n_stems"], "n_pairs": m["n_pairs"]}
    for b in R.DISTANCE_LABELS:
        out[f"ratio_{b}"] = _r(m[f"p_{b}"] / m[f"e_{b}"]) if m[f"e_{b}"] > 0 else float("nan")
    out["x_far"], out["x_far_nl"] = _r(m["x_far"]), _r(m["x_far_nl"])
    for c in ("section", "language", "hand"):
        out[f"far_{c}_excess"] = _r(m[f"far_{c}"] - m[f"far_{c}_exp"])
    for c in R.FORM_LABELS:
        out[f"form_{c}"] = _r(m[f"form_same_{c}"] - m[f"form_exp_{c}"])
    return out


def loss_004(df: pd.DataFrame, seed: int) -> dict:
    """004 loss B inputs and value of a sample against the real ZL split point of 004."""
    M4 = _load("exp004", "004_selfcitation_port")
    st = M4.stats_for(df, seed, full=False)
    target = M4.real_point("split")
    out = {}
    loss = 0.0
    for f, w in M4.TARGETS.items():
        comp = ((st[f] - float(target[f])) / float(target[f])) ** 2
        out[f"t_{f}"] = _r(st[f])
        out[f"compB_{f}"] = _r(comp)
        loss += w * comp
    out["loss_B"] = _r(loss)
    return out


def mem_grid() -> list[dict]:
    return [{"p_prev_page": p, "page_memory": m} for p in MEM_GRID["p_prev_page"] for m in MEM_GRID["page_memory"]]


def cmd_memory(args) -> None:
    ctx = Context("ZL_split")
    parts = OUT / "runs" / "memory"
    parts.mkdir(parents=True, exist_ok=True)
    params21 = M6.selfcitation_params()
    real = recurrence_ratios(ctx, ctx.real)
    if args.phase == "table":
        write_calibration(parts, real)
        return
    pd.DataFrame([{"source": "real", **real}]).to_csv(parts / "real_profile.csv", index=False)
    if args.phase == "grid":
        grid = mem_grid()

        def one(seed: int) -> dict:
            cfg, rep = divmod(seed, 100)
            p = grid[cfg]
            sample = selfcitation_like(ctx.real, seed=MEM_GRID_SEED0 + seed, **params21, **p)
            return {"config": cfg, "rep": rep, **p, **recurrence_ratios(ctx, sample)}

        seeds = [c * 100 + r for c in range(len(grid)) for r in range(MEM_GRID_SEEDS)]
        run_seeds("018_memory_grid", seeds, one, parts / "grid",
                  params={"grid": MEM_GRID, "config21": params21, "seeds_per_config": MEM_GRID_SEEDS,
                          "seed_offset": MEM_GRID_SEED0, "loss": "sum of squared deviations of the bin ratios "
                          + "/".join(CAL_BINS) + " from the real ZL split"}, progress=progress)
        write_calibration(parts, real)
        return
    if args.phase == "winner":
        cal = pd.read_csv(OUT / "memory_calibration.csv")
        win = cal[cal["winner"] == 1].iloc[0]
        p = {"p_prev_page": float(win["p_prev_page"]), "page_memory": float(win["page_memory"])}

        def one(seed: int) -> dict:
            sample = selfcitation_like(ctx.real, seed=TSM_SEED0 + seed, **params21, **p)
            return {**p, **recurrence_ratios(ctx, sample), **loss_004(sample, seed)}

        def base(seed: int) -> dict:
            sample = selfcitation_like(ctx.real, seed=TS_SEED0 + seed, **params21)
            return {"p_prev_page": 0.0, "page_memory": 1.0, **recurrence_ratios(ctx, sample), **loss_004(sample, seed)}

        run_seeds("018_memory_winner", list(range(MEM_WINNER_SEEDS)), one, parts / "winner",
                  params={"memory": p, "config21": params21, "seed_offset": TSM_SEED0}, progress=progress)
        run_seeds("018_memory_base21", list(range(MEM_WINNER_SEEDS)), base, parts / "base21",
                  params={"config21": params21, "seed_offset": TS_SEED0}, progress=progress)
        write_calibration(parts, real)
        return
    raise ValueError(args.phase)


def write_calibration(parts: Path, real: dict) -> None:
    """``memory_calibration.csv``: grid means, loss on the bins 1 / 2-3 / 4-7 / 8-15, winner, 30-seed rows."""
    g = load_runs(parts / "grid")
    if g.empty:
        return
    cols = [c for c in g.columns if c.startswith("ratio_") or c.startswith("x_far") or c.startswith("far_")
            or c.startswith("form_") or c in ("n_stems", "n_pairs")]
    means = g.groupby(["config", "p_prev_page", "page_memory"])[cols].mean().reset_index()
    means["n_seeds"] = g.groupby("config").size().to_numpy()
    means["loss_bins"] = sum((means[f"ratio_{b}"] - real[f"ratio_{b}"]) ** 2 for b in CAL_BINS)
    means["phase"] = "grid"
    means = means.sort_values("loss_bins", kind="stable").reset_index(drop=True)
    means["rank"] = np.arange(1, len(means) + 1)
    means["winner"] = (means["rank"] == 1).astype(int)
    rows = [means]
    for name in ("winner", "base21"):
        w = load_runs(parts / name)
        if w.empty:
            continue
        wc = [c for c in cols if c in w.columns] + [c for c in w.columns if c.startswith("t_") or c.startswith("compB_") or c == "loss_B"]
        m = w[wc].mean().to_dict()
        lo = w[wc].quantile(0.025).to_dict()
        hi = w[wc].quantile(0.975).to_dict()
        row = {"config": -1, "p_prev_page": float(w["p_prev_page"].iloc[0]), "page_memory": float(w["page_memory"].iloc[0]),
               **m, "n_seeds": int(len(w)), "phase": name, "winner": 0, "rank": 0}
        row["loss_bins"] = sum((m[f"ratio_{b}"] - real[f"ratio_{b}"]) ** 2 for b in CAL_BINS)
        for b in R.DISTANCE_LABELS:
            row[f"ratio_{b}_lo"], row[f"ratio_{b}_hi"] = lo[f"ratio_{b}"], hi[f"ratio_{b}"]
        if "loss_B" in m:
            row["loss_B_lo"], row["loss_B_hi"] = lo["loss_B"], hi["loss_B"]
            # 004 computes the loss on the seed means of the targets (not the mean of per-seed losses)
            M4 = _load("exp004", "004_selfcitation_port")
            target = M4.real_point("split")
            row["loss_B_of_means"] = float(sum(
                wgt * ((m[f"t_{f}"] - float(target[f])) / float(target[f])) ** 2 for f, wgt in M4.TARGETS.items()))
            row["far_section_excess_lo"], row["far_section_excess_hi"] = lo["far_section_excess"], hi["far_section_excess"]
        rows.append(pd.DataFrame([row]))
    rows.append(pd.DataFrame([{"config": -2, "p_prev_page": np.nan, "page_memory": np.nan, **real, "n_seeds": 1,
                               "phase": "real", "winner": 0, "rank": 0, "loss_bins": 0.0}]))
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "memory_calibration.csv", index=False)


# --------------------------------------------------------------------------- #
# extras of the real text: residual matrices, change points, neighbours, quire matrix
# --------------------------------------------------------------------------- #


def cmd_extras(args) -> None:
    ctx = Context(args.combo)
    tag = args.combo.replace(":", "_")
    parts = OUT / "runs" / tag
    parts.mkdir(parents=True, exist_ok=True)
    counts = ctx.counts(ctx.real)
    weight = ctx.pages["weight"].to_numpy()
    cps, nbs, qms, mats = [], [], [], {}
    for cut, (idx, mask, ref) in ctx.cuts.items():
        sub, w = counts[idx], weight[idx] / weight[idx].sum()
        dist = ctx.dist[np.ix_(idx, idx)]
        labels = ctx.pages.iloc[idx]
        folios = labels.index.to_numpy()
        refs = {c: L.boundaries_of(labels[c].astype(str).to_numpy()) for c in REFS}
        k_ref = len(refs[ref]) + 1
        for measure in MEASURES:
            n_mat, w_mat = L.similarity_matrices(measure, sub, w)
            r_full, w_full, f = L.residuals(n_mat, w_mat, dist, mask, DMAX)
            w_model = f[dist] * w_full if measure == "rare" else w_full
            num, den, kind = L.segmentation_inputs(r_full, w_model, measure)
            segs = L.segment(num, den, K_MAX, kind)
            if cut == "all":
                with np.errstate(invalid="ignore", divide="ignore"):
                    z = np.where(w_model > 0, r_full / np.where(w_model > 0, w_model, 1.0), np.nan)
                mats[measure] = pd.DataFrame(z, index=folios, columns=folios)
                if measure == "rare":   # 7 x 7 page window pooled residual (sum R / sum fW): readable in a figure
                    mats["rare_smooth"] = pd.DataFrame(L.window_pooled(r_full, w_model, 3), index=folios, columns=folios)
            for m in CP_MULT:
                k = min(max(m * k_ref, 2), K_MAX)
                for b in segs[k]:
                    row = {"combo": args.combo, "cut": cut, "measure": measure, "K_mult": m, "K": k, "position": int(b),
                           "folio_before": folios[b - 1], "folio_after": folios[b]}
                    for c in REFS:
                        row[f"hits_{c}"] = int((np.abs(refs[c] - b) <= 1).any()) if len(refs[c]) else 0
                        row[f"{c}_before"], row[f"{c}_after"] = labels[c].iloc[b - 1], labels[c].iloc[b]
                    cps.append(row)
            if cut != "all":
                continue
            lang = labels["language"].astype(str).to_numpy()
            allowed = (lang[:, None] == lang[None, :]) & (w_full > 0)
            if measure == "rare":
                score = np.where(w_model > 0, r_full / np.sqrt(np.where(w_model > 0, w_model, 1.0)), -np.inf)
            else:
                score = np.where(w_full > 0, r_full, -np.inf)
            nb = L.neighbours(score, dist, allowed, 3, L.NEAR_EDGES, 8)
            nd = nb["neighbour_dist"]
            for i in range(len(folios)):
                if nd[i, 0] < 0:
                    continue
                cand = np.flatnonzero(allowed[i])
                cand = cand[cand != i]
                top = cand[np.argsort(-score[i, cand], kind="stable")[:3]]
                nbs.append({"combo": args.combo, "measure": measure, "folio": folios[i], "rank": int(labels["rank"].iloc[i]),
                            "section": labels["section"].iloc[i], "language": labels["language"].iloc[i],
                            "hand": labels["hand"].iloc[i], "quire": labels["quire"].iloc[i],
                            "neighbours": " ".join(folios[top]), "distances": " ".join(str(int(x)) for x in nd[i]),
                            "min_dist": int(nd[i].min()), "anomalous": int(nd[i].min() > 8)})
            q = L.quire_matrix(r_full, w_model, labels["quire"].astype(str).to_numpy())
            qm = q.reset_index().rename(columns={"index": "quire"})
            qm.insert(0, "measure", measure)
            qm.insert(0, "combo", args.combo)
            qms.append(qm)
    pd.DataFrame(cps).to_csv(parts / "changepoints.csv", index=False)
    pd.DataFrame(nbs).to_csv(parts / "neighbours.csv", index=False)
    pd.concat(qms).to_csv(parts / "quire_matrix.csv", index=False)
    for measure, m in mats.items():
        m.round(4).to_csv(parts / f"residual_{measure}.csv")
    print(f"extras written to {parts}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--combo", required=True)
    p.add_argument("--sources", default="all")
    p.add_argument("--seeds", default=None, help="'0:30', '0:200/4:1' or '3,5'")
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("memory")
    p.add_argument("--phase", default="grid", choices=["grid", "winner", "table"])
    p.set_defaults(func=cmd_memory)
    p = sub.add_parser("extras")
    p.add_argument("--combo", required=True)
    p.set_defaults(func=cmd_extras)
    p = sub.add_parser("summarize")
    p.set_defaults(func=lambda a: __import__("summarize").summarize())
    args = ap.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    args.func(args)


if __name__ == "__main__":
    main()
