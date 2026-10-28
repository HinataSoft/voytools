"""Experiment 017: does a local rewrite give position, seam and word cohesion at once?

ts2020 (004 config 21) post-processed with the conditional place tables of ``voytools.placerules``
(:func:`learn_conditional_rules` / :func:`apply_conditional_rules`; R1-R3 of 014 unchanged), measured with
the full profile of 009 (``M9.profile``), 010 (``M10.line_profile``), 011 (``M11.measure``) and the
calibration targets of 004 (``M14.calib_metrics``) - the same code path as 014 (``measure_all``).

Sets (``SETS``; 30 seeds each, ZL split and join): ``base`` plain ts2020; ``R12`` +R1+R2 and ``R123``
+R1+R2+R3 (014, recomputed); ``R12pA`` +R1+R2'+R3'a; ``R12pB`` +R1+R2'+R3'b; ``R12pC`` +R1+R2'+R3'c;
``R12pA3`` +R1+R2'+R3'a+R3''; ``R12pA3L`` the same with the type licence.  Tables are learned on the real
text with the affix lists (b) of 011 (as in 014).  Further commands: ``template`` (the same on the 005
template with topic, sets ``base`` / ``R12pA`` / ``R12pA3``, ZL split only - the template of 005 is
calibrated for ZL_split and GC_split), ``selfcheck`` (set ``R12pA`` with tables learned on another ts2020
sample, seed + 1000, instead of the real text), ``bootstrap`` (page bootstrap of the real text, 100 draws,
same metrics), ``rules`` (table support and the learned tables), ``summarize`` (tables, lists, plots).

Usage (from the repository root)::

    python experiments/017_conditional_seam_rewrite/run.py rules
    python experiments/017_conditional_seam_rewrite/run.py augment --set R12pA --tag ZL_split [--seeds 0:30]
    python experiments/017_conditional_seam_rewrite/run.py template --set R12pA --tag ZL_split
    python experiments/017_conditional_seam_rewrite/run.py selfcheck --tag ZL_split
    python experiments/017_conditional_seam_rewrite/run.py bootstrap --tag ZL_split [--seeds -1:100]
    python experiments/017_conditional_seam_rewrite/run.py summarize
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import placerules as PR  # noqa: E402
from voytools import run_seeds  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M14 = _load("exp014", "014_place_rules_inventory")
M9 = M14.M9

OUT = ROOT / "results" / "017_conditional_seam_rewrite"
TAGS = ("ZL_split", "ZL_join")          # ts2020 is Eva-only
TEMPLATE_TAGS = ("ZL_split",)           # the 005 template has a calibrated configuration for ZL_split only
N_NULL = 30
N_BOOT = 100
MIN_CELL = 5
MIN_CONTEXT = M14.MIN_COUNT             # R3 of 014
END_RATIO = 5.0
LICENCE_MIN = 2
MAX_TRIES = 5
SELFCHECK_OFFSET = 1000
# set name -> (rules, licence); 'base' = the plain generator through the same code path
SETS = {
    "base": ((), False),
    "R12": (("R1", "R2"), False),
    "R123": (("R1", "R2", "R3"), False),
    "R12pA": (("R1", "R2p", "R3a"), False),
    "R12pB": (("R1", "R2p", "R3b"), False),
    "R12pC": (("R1", "R2p", "R3c"), False),
    "R12pA3": (("R1", "R2p", "R3a", "R3pp"), False),
    "R12pA3L": (("R1", "R2p", "R3a", "R3pp"), True),
}
OLD_SETS = ("R12", "R123")              # applied with apply_place_rules (014 code path)
TEMPLATE_SETS = ("base", "R12pA", "R12pA3")
RULE_NAMES = ("R1", "R2", "R2p", "R3", "R3a", "R3b", "R3c", "R3pp")
LEVELS = ("R3", "R3a", "R3b", "R3c", "R3pp", "uncond:R3", "uncond:R3pp")


# --------------------------------------------------------------------------- #
# rules and sets
# --------------------------------------------------------------------------- #

_RULES: dict = {}


def rules_for(ctx) -> PR.ConditionalRules:
    """Conditional tables learned on the real text of ``ctx`` with the affix lists (b)."""
    if ctx.source not in _RULES:
        pre, suf = M14.lists_of(ctx, "b")
        _RULES[ctx.source] = PR.learn_conditional_rules(ctx.real, pre, suf, preset="extended", min_cell=MIN_CELL,
                                                        min_context=MIN_CONTEXT, end_ratio=END_RATIO,
                                                        licence_min_count=LICENCE_MIN)
    return _RULES[ctx.source]


def apply_set(gen: pd.DataFrame, rules: PR.ConditionalRules, set_name: str, seed: int) -> pd.DataFrame:
    which, licence = SETS[set_name]
    if set_name == "base":
        out = gen.copy()
        out.attrs.update(gen.attrs)
        out.attrs["place_rules_changed"] = {}
        out.attrs["rewrite_report"] = {}
        out.attrs["n_tokens_changed"] = 0
        return out
    if set_name in OLD_SETS:
        out = PR.apply_place_rules(gen, rules.base, which, seed=seed)
        out.attrs["rewrite_report"] = {k: {"n_changed": v} for k, v in out.attrs["place_rules_changed"].items()}
        out.attrs["n_tokens_changed"] = int((out["token"].to_numpy() != gen["token"].to_numpy()).sum())
        return out
    return PR.apply_conditional_rules(gen, rules, which, seed=seed, licence=licence, max_tries=MAX_TRIES)


def cost_columns(gen: pd.DataFrame, aug: pd.DataFrame) -> dict:
    """Rewrite cost of one corpus: changed tokens, backoff support, licence rejections, new types."""
    rep = aug.attrs.get("rewrite_report", {})
    row = {"n_tokens": int(len(aug)), "n_tokens_changed": int(aug.attrs.get("n_tokens_changed", 0))}
    row["share_tokens_changed"] = row["n_tokens_changed"] / len(aug) if len(aug) else np.nan
    for r in RULE_NAMES:
        d = rep.get(r, {})
        row[f"n_cand_{r}"] = int(d.get("n_candidates", 0))
        row[f"n_changed_{r}"] = int(d.get("n_changed", 0))
        row[f"n_rejected_{r}"] = int(d.get("n_rejected", 0))
    sup = {lv: 0 for lv in LEVELS}
    n_cond = 0
    for r in ("R3", "R3a", "R3b", "R3c", "R3pp"):
        for lv, n in rep.get(r, {}).get("support", {}).items():
            sup[lv] = sup.get(lv, 0) + int(n)
            n_cond += int(n)
    for lv, n in sup.items():
        row[f"support_{lv.replace(':', '_')}"] = n
    full = sum(int(rep.get(r, {}).get("support", {}).get(r, 0)) for r in ("R3", "R3a", "R3b", "R3c", "R3pp"))
    row["n_conditional_candidates"] = n_cond
    row["share_backoff"] = (n_cond - full) / n_cond if n_cond else np.nan
    row["n_rejected_total"] = sum(row[f"n_rejected_{r}"] for r in RULE_NAMES)
    row["n_cand_total"] = sum(row[f"n_cand_{r}"] for r in RULE_NAMES)
    row["share_rejected"] = row["n_rejected_total"] / row["n_cand_total"] if row["n_cand_total"] else np.nan
    v0 = set(gen["token"].tolist())
    v1 = pd.Series(aug["token"].tolist()).value_counts()
    row["n_types_base"] = int(len(v0))
    row["n_types"] = int(len(v1))
    row["n_types_new"] = int((~v1.index.isin(list(v0))).sum())
    row["n_tokens_new_types"] = int(v1[~v1.index.isin(list(v0))].sum())
    row["end_glyphs"] = " ".join(rep.get("R2p", {}).get("end_glyphs", []))
    return row


def measure_row(source: str, kind: str, gen: pd.DataFrame, aug: pd.DataFrame, ctx, seed: int, set_name: str,
                tag: str, generator: str) -> dict:
    row = {"source": source, "kind": kind, "set": set_name, "tag": tag, "generator": generator,
           "licence": bool(SETS[set_name][1]), "rules": "+".join(SETS[set_name][0])}
    row.update(cost_columns(gen, aug))
    row.update(M14.measure_all(aug, ctx, seed))
    return row


def params_for(ctx, rules: PR.ConditionalRules, set_name: str, extra: dict | None = None) -> dict:
    which, licence = SETS[set_name]
    p = {"set": set_name, "rules": list(which), "licence": licence, "min_cell": MIN_CELL, "min_context": MIN_CONTEXT,
         "end_ratio": END_RATIO, "end_glyphs": list(rules.end_glyphs), "licence_min_count": LICENCE_MIN,
         "max_tries": MAX_TRIES, "ts2020_config": M9.TS2020_CONFIG, "n_null": N_NULL,
         "rule_lists": {"prefixes": rules.prefixes, "suffixes": rules.suffixes},
         "profile_prefixes": ctx.prefixes, "profile_suffixes": ctx.suffixes}
    p.update(extra or {})
    return p


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def cmd_augment(args) -> None:
    tag, set_name = args.tag, args.set
    ctx = M14.ctx_of(tag)
    rules = rules_for(ctx)
    source = f"ts2020+{set_name}:{tag}" if set_name != "base" else f"ts2020:{tag}"

    def one(seed: int) -> dict:
        gen = ctx.corpus("ts2020", seed)
        aug = apply_set(gen, rules, set_name, seed)
        return measure_row(source, "null", gen, aug, ctx, seed, set_name, tag, "ts2020")

    seeds = M9.parse_seeds(args.seeds, list(range(N_NULL)))
    run_seeds("017_augment", seeds, one, OUT / "runs" / f"{set_name}__{tag}",
              params=params_for(ctx, rules, set_name, {"generator": "ts2020"}), progress=M9.progress)


def cmd_template(args) -> None:
    tag, set_name = args.tag, args.set
    ctx = M14.ctx_of(tag)
    rules = rules_for(ctx)
    source = f"tmpl+{set_name}:{tag}" if set_name != "base" else f"tmpl:{tag}"

    def one(seed: int) -> dict:
        gen = ctx.corpus("tmpl", seed)
        aug = apply_set(gen, rules, set_name, seed)
        return measure_row(source, "null", gen, aug, ctx, seed, set_name, tag, "tmpl")

    seeds = M9.parse_seeds(args.seeds, list(range(N_NULL)))
    run_seeds("017_template", seeds, one, OUT / "runs" / f"tmpl_{set_name}__{tag}",
              params=params_for(ctx, rules, set_name, {"generator": "tmpl", "tmpl_seed0": M9.TMPL_SEED0}),
              progress=M9.progress)


def cmd_selfcheck(args) -> None:
    """Set R12pA with tables learned on ts2020 (seed + 1000) instead of the real text, applied to ts2020 (seed)."""
    tag, set_name = args.tag, "R12pA"
    ctx = M14.ctx_of(tag)
    pre, suf = M14.lists_of(ctx, "b")
    source = f"selfcheck+{set_name}:{tag}"

    def one(seed: int) -> dict:
        train = ctx.corpus("ts2020", SELFCHECK_OFFSET + seed)
        rules = PR.learn_conditional_rules(train, pre, suf, preset="extended", min_cell=MIN_CELL, min_context=MIN_CONTEXT,
                                           end_ratio=END_RATIO, licence_min_count=LICENCE_MIN)
        gen = ctx.corpus("ts2020", seed)
        aug = apply_set(gen, rules, set_name, seed)
        row = measure_row(source, "null", gen, aug, ctx, seed, set_name, tag, "ts2020")
        row["train_seed"] = SELFCHECK_OFFSET + seed
        return row

    seeds = M9.parse_seeds(args.seeds, list(range(N_NULL)))
    rules0 = rules_for(ctx)
    run_seeds("017_selfcheck", seeds, one, OUT / "runs" / f"selfcheck_{set_name}__{tag}",
              params=params_for(ctx, rules0, set_name, {"generator": "ts2020", "tables_from": "ts2020 seed+1000"}),
              progress=M9.progress)


def page_bootstrap(df: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, dict]:
    """Pages drawn with replacement (as many as there are); the k-th copy of a page gets folio ``<folio>~k``
    so that its lines stay separate units.  Returns the table in page order and the extra folio -> original."""
    pages = list(pd.unique(df["folio"].to_numpy()))
    rng = np.random.default_rng(seed)
    draw = rng.choice(len(pages), size=len(pages), replace=True)
    counts = np.bincount(draw, minlength=len(pages))
    parts, renamed = [], {}
    for i, p in enumerate(pages):
        block = df[df["folio"] == p]
        for k in range(counts[i]):
            b = block.copy()
            if k > 0:
                new = f"{p}~{k}"
                b["folio"] = new
                renamed[new] = p
            parts.append(b)
    out = pd.concat(parts, ignore_index=True)
    out.attrs.update(df.attrs)
    out.attrs["source"] = f"{df.attrs.get('source', 'real')}:bootstrap{seed}"
    return out, renamed


def cmd_bootstrap(args) -> None:
    tag = args.tag
    ctx = M14.ctx_of(tag)

    def one(seed: int) -> dict:
        if seed < 0:
            df, c = ctx.real, ctx
            kind = "point"
        else:
            df, renamed = page_bootstrap(ctx.real, seed)
            c = copy.copy(ctx)
            # page order of the bootstrap table (copies next to their original), re-ranked consecutively so
            # that the locality metrics see distinct positions (a copy sits at distance 1 from its original)
            c.rank = {p: i for i, p in enumerate(pd.unique(df["folio"].to_numpy()))}
            kind = "boot"
        row = {"source": tag, "kind": kind, "set": "real", "tag": tag, "generator": "real", "licence": False, "rules": "",
               "n_tokens": int(len(df)), "n_pages_distinct": int(df["folio"].nunique())}
        row.update(M14.measure_all(df, c, seed))
        return row

    seeds = M9.parse_seeds(args.seeds, [-1] + list(range(N_BOOT)))
    run_seeds("017_bootstrap", seeds, one, OUT / "runs" / f"real__{tag}",
              params={"tag": tag, "n_boot": N_BOOT, "unit": "page (with replacement, copies renamed)",
                      "profile_prefixes": ctx.prefixes, "profile_suffixes": ctx.suffixes}, progress=M9.progress)


def cmd_rules(args) -> None:
    """``table_support.csv`` (real tables: contexts, supported cells, share of real tokens in supported cells)
    and ``rules_tables.csv`` (the learned conditional tables, ZL split and join)."""
    OUT.mkdir(parents=True, exist_ok=True)
    sup, tabs = [], []
    for tag in TAGS:
        ctx = M14.ctx_of(tag)
        r = rules_for(ctx)
        s = r.support()
        s.insert(0, "source", tag)
        s["end_glyphs"] = " ".join(r.end_glyphs)
        s["n_vocab_types_ge2"] = int(sum(1 for v in r.vocab.values() if v >= LICENCE_MIN))
        sup.append(s)
        for name in ("R3", "R3a", "R3b", "R3pp"):
            t = r.table(name)
            t.insert(0, "source", tag)
            tabs.append(t)
    out = pd.concat(sup, ignore_index=True)
    out.to_csv(OUT / "table_support.csv", index=False, float_format="%.5g")
    t = pd.concat(tabs, ignore_index=True)
    t["item"] = t["item"].where(t["item"] != "", M14.EMPTY)
    t["context"] = t["context"].str.replace("||", "|<none>|", regex=False).str.replace(r"\|$", "|<none>", regex=True)
    t.to_csv(OUT / "rules_tables.csv", index=False, float_format="%.5g")
    print(out.to_string())


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("rules")
    for name, sets, tags in (("augment", list(SETS), TAGS), ("template", list(TEMPLATE_SETS), TEMPLATE_TAGS)):
        p = sub.add_parser(name)
        p.add_argument("--set", required=True, choices=sets)
        p.add_argument("--tag", required=True, choices=tags)
        p.add_argument("--seeds", default=None)
    for name in ("selfcheck", "bootstrap"):
        p = sub.add_parser(name)
        p.add_argument("--tag", required=True, choices=TAGS)
        p.add_argument("--seeds", default=None)
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from summarize import summarize  # noqa: E402  (experiments/017_conditional_seam_rewrite/summarize.py)

        summarize()
        return
    {"rules": cmd_rules, "augment": cmd_augment, "template": cmd_template, "selfcheck": cmd_selfcheck,
     "bootstrap": cmd_bootstrap}[args.cmd](args)


if __name__ == "__main__":
    main()
