"""Experiment 014: inventory of place rules - how many items carry position and the word seam.

Steps (``spec.md``), all with ``voytools.placerules``:

1. ``inventory``: pointwise decomposition of I(prefix | suffix | first glyph | last glyph; position class
   of 007) into affixes / glyphs, position classes and pairs, with the number of items covering
   50 / 80 / 95 % of the information; chain decomposition edge of line -> paragraph-first line -> drift
   inside the line; the same for the word seam (last glyph -> first glyph of the next token, with and
   without y->q, suffix -> next prefix), and the test "is the seam only a glyph pair"
   (I(suffix; next prefix | edge glyphs)).  Sources: ZL/GC x split/join, decompositions (a) learned
   affixes and (b) learned + positional glyphs (011), the ts2020 null (3 seeds) and - step 5 - the
   catalogue (Linnaeus per item) and the gibberish in its own lines (controls of 009/013; skipped
   when their repositories are not available).
2. ``stability``: the same tables per Currier language, scribe (>= 3 000 tokens) and section; JSD
   between groups versus between random halves (lines) of one group, 200 draws, and the share of
   rules (pointwise terms covering 80 %) with the same direction in both groups.
3. ``augment``: ts2020 (004 config 21) post-processed with the learned tables R1 / R1+R2 / R1+R2+R3,
   measured with the full profile of 009 (``M9.profile``), 010 (``M10.line_profile``) and 011
   (``M11.measure``), 30 seeds, ZL split and join; the plain ts2020 of the same seeds is taken over
   from ``results/009..011``.
4. ``summarize``: tables, lists (what the augmented generator catches up with, what it breaks), plots.

Usage (from the repository root)::

    python experiments/014_place_rules_inventory/run.py inventory
    python experiments/014_place_rules_inventory/run.py stability
    python experiments/014_place_rules_inventory/run.py augment --variant R1 --tag ZL_split [--seeds 0:30]
    python experiments/014_place_rules_inventory/run.py summarize
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import placerules as PR  # noqa: E402
from voytools import run_seeds  # noqa: E402
from voytools._common import glyph_lists, line_ids  # noqa: E402
from voytools.experiment import write_meta  # noqa: E402
from voytools.morphology import decompose_table  # noqa: E402
from voytools.stats import conditional_entropy  # noqa: E402
from voytools.variation import line_place, position_class  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / rel / "run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M11 = _load("exp011", "011_suffix_variation")
M9 = M11.M9
M10 = _load("exp010", "010_line_as_item")

OUT = ROOT / "results" / "014_place_rules_inventory"
REAL_TAGS = ("ZL_split", "ZL_join", "GC_split", "GC_join")
AUG_TAGS = ("ZL_split", "ZL_join")          # ts2020 is Eva-only
VARIANTS = {"R1": ("R1",), "R12": ("R1", "R2"), "R123": ("R1", "R2", "R3")}
TARGETS = ("prefix", "suffix", "first_glyph", "last_glyph")
N_TS_INVENTORY = 3
N_NULL = 30
N_REP = 200
MIN_HAND = 3000
MIN_SECTION = 2000
MIN_COUNT = 20
SHUFFLES = {"chain": 50, "seam": 50}
CONTROLS = {"ctl:linnaeus_items": "catalogue per item", "ctl:gibberish_native": "gibberish in own lines"}
POS_AFFIXES = M11.POS_AFFIXES
EMPTY = "<none>"            # label of the empty affix in the CSV outputs


# --------------------------------------------------------------------------- #
# sources
# --------------------------------------------------------------------------- #

_CTX: dict = {}


def ctx_of(tag: str):
    if tag not in _CTX:
        _CTX[tag] = M9.Ctx(tag)
    return _CTX[tag]


def lists_of(ctx, variant: str) -> tuple[list, list]:
    """Affix lists of decomposition (a) or (b) (011 ``variants``; controls have only (a))."""
    v = M11.variants(ctx)
    return v[variant if variant in v else "a"]


def control_table(name: str) -> pd.DataFrame | None:
    try:
        if name == "ctl:gibberish_native":
            M13 = _load("exp013", "013_gibberish_and_cipher_controls")
            return M13.Ctx("gibberish_native").real
        return ctx_of(name).real
    except Exception as exc:  # repository of the control not available
        print(f"  {name}: not available ({type(exc).__name__}: {exc})", flush=True)
        return None


def source_tables() -> list[tuple[str, pd.DataFrame, object]]:
    """(source name, token table, context with affix lists) of every inventory source."""
    out = []
    for tag in REAL_TAGS:
        c = ctx_of(tag)
        out.append((tag, c.real, c))
    c = ctx_of("ZL_split")
    for s in range(N_TS_INVENTORY):
        out.append((f"ts2020:ZL_split:s{s}", c.corpus("ts2020", s), c))
    for name in CONTROLS:
        t = control_table(name)
        if t is not None:
            out.append((name, t.reset_index(drop=True), ctx_of(name) if name != "ctl:gibberish_native" else None))
    return out


# --------------------------------------------------------------------------- #
# step 1: inventory
# --------------------------------------------------------------------------- #


def _preset(df: pd.DataFrame):
    return "extended" if str(df.attrs.get("alphabet", "")).lower().startswith("eva") else None


def labels_for(df: pd.DataFrame, ctx, variant: str) -> dict:
    if ctx is None:
        d = decompose_table(df)
        pre, suf = d.attrs["prefixes"], d.attrs["suffixes"]
    else:
        pre, suf = lists_of(ctx, variant)
    return PR.token_labels(df, pre, suf, _preset(df))


def _cov_row(base: dict, items: str, cov: dict) -> dict:
    return {**base, "items": items, **cov}


def inventory_one(source: str, df: pd.DataFrame, ctx, variant: str, seed: int = 0) -> dict:
    """Position and seam decompositions of one corpus; returns frames keyed by output file."""
    lab = labels_for(df, ctx, variant)
    pos = position_class(df)
    base = {"source": source, "variant": variant}
    contrib, cover, chain, seam, seam_test = [], [], [], [], []
    for target in TARGETS:
        if variant == "b" and target.endswith("glyph"):
            continue  # glyph targets do not depend on the affix lists
        v = "-" if target.endswith("glyph") else variant
        b = {**base, "variant": v, "decomposition": "position", "target": target}
        d = PR.item_decomposition(lab[target], pos, "item", "pos")
        ch = PR.position_chain(df, lab[target], shuffles=SHUFFLES["chain"], seed=seed)
        chain.append({**b, **ch})
        for name, items in (("by_item", d["by_x"]), ("by_pos", d["by_y"])):
            t = items.rename(columns={"item": "x", "pos": "x"}).assign(**b, items=name)
            contrib.append(t)
        contrib.append(d["pairs"].rename(columns={"item": "x", "pos": "y"}).assign(**b, items="pair"))
        cov = d["coverage"].replace({"items": {"item": "by_item", "pos": "by_pos"}})
        for _, r in cov.iterrows():
            cover.append({**b, **r.to_dict(), "mi_corrected": ch["I_pos"]})
    # seam: edge glyphs
    if variant == "a":
        L, R, _ = PR.seam_pairs(df, lab["last_glyph"], lab["first_glyph"])
        d = PR.item_decomposition(L, R, "last", "first")
        b = {**base, "variant": "-", "decomposition": "seam_glyph", "target": "last_glyph->first_glyph"}
        seam.append(d["pairs"].assign(**b, items="pair"))
        seam.append(d["by_x"].rename(columns={"last": "x"}).assign(**b, items="by_last"))
        seam.append(d["by_y"].rename(columns={"first": "x"}).assign(**b, items="by_first"))
        st = PR.seam_glyph_vs_affix(df, lab, shuffles=SHUFFLES["seam"], seed=seed)
        for _, r in d["coverage"].replace({"items": {"last": "by_last", "first": "by_first"}}).iterrows():
            cover.append({**b, **r.to_dict(), "mi_corrected": st["I_glyph"]})
        # without the strongest pair (y -> q in Eva; its v101 counterpart in GC): the pair's term is removed
        # from the numerator and from the total
        top = d["pairs"].iloc[:1]
        c_yq = float(top["contribution"].sum())
        rest = d["pairs"].drop(index=top.index)
        b2 = {**b, "decomposition": "seam_glyph_no_yq", "removed_pair": f"{top['last'].iloc[0]}>{top['first'].iloc[0]}"}
        cover.append(_cov_row({**b2, "mi_corrected": st["I_glyph"] - c_yq, "share_yq": c_yq / d["mi"] if d["mi"] else np.nan},
                              "pair", PR.coverage(rest["contribution"], total=d["mi"] - c_yq)))
        for by, name in (("last", "by_last"), ("first", "by_first")):
            m = PR.marginal_contributions(rest, by)
            cover.append(_cov_row({**b2, "mi_corrected": st["I_glyph"] - c_yq}, name,
                                  PR.coverage(m["contribution"], total=d["mi"] - c_yq)))
    # seam: suffix -> next prefix (depends on the lists)
    L, R, _ = PR.seam_pairs(df, lab["suffix"], lab["prefix"])
    d = PR.item_decomposition(L, R, "suffix", "prefix")
    st = PR.seam_glyph_vs_affix(df, lab, shuffles=SHUFFLES["seam"], seed=seed)
    b = {**base, "decomposition": "seam_affix", "target": "suffix->prefix"}
    seam.append(d["pairs"].assign(**b, items="pair"))
    seam.append(d["by_x"].rename(columns={"suffix": "x"}).assign(**b, items="by_suffix"))
    seam.append(d["by_y"].rename(columns={"prefix": "x"}).assign(**b, items="by_prefix"))
    for _, r in d["coverage"].replace({"items": {"suffix": "by_suffix", "prefix": "by_prefix"}}).iterrows():
        cover.append({**b, **r.to_dict(), "mi_corrected": st["I_affix"]})
    seam_test.append({**base, "n_tokens": len(df), **st})
    return {"position_contributions": pd.concat(contrib, ignore_index=True),
            "seam_contributions": pd.concat(seam, ignore_index=True),
            "coverage": pd.DataFrame(cover), "position_chain": pd.DataFrame(chain), "seam_test": pd.DataFrame(seam_test)}


def cmd_inventory(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    acc: dict[str, list] = {}
    t0 = time.time()
    for source, df, ctx in source_tables():
        for variant in ("a", "b") if ctx is not None and not source.startswith("ctl:") else ("a",):
            print(f"{source} ({variant}) {len(df)} tokens", flush=True)
            res = inventory_one(source, df, ctx, variant)
            for k, v in res.items():
                acc.setdefault(k, []).append(v)
    for k, frames in acc.items():
        out = pd.concat(frames, ignore_index=True)
        lead = [c for c in ("source", "variant", "decomposition", "target", "items") if c in out.columns]
        out = out[lead + [c for c in out.columns if c not in lead]]
        for c in ("x", "y"):
            if c in out.columns:
                out[c] = out[c].where(out[c] != "", EMPTY)
        out.to_csv(OUT / f"{k}.csv", index=False, float_format="%.6g")
        print(f"{k}.csv: {len(out)} rows", flush=True)
    print(f"inventory done in {time.time() - t0:.0f} s")


# --------------------------------------------------------------------------- #
# step 2: stability
# --------------------------------------------------------------------------- #


def r1_context_all(df: pd.DataFrame) -> np.ndarray:
    c = PR._r1_context(df)
    return np.where(c == "", "other", c).astype(object)


def stability_tables(df: pd.DataFrame, lab: dict) -> dict:
    """name -> (x, y, row mask or None, contexts, unit ids) of the rule tables compared across groups."""
    ctx1 = r1_context_all(df)
    place = line_place(df)
    lid = line_ids(df)
    rows = np.arange(len(df))
    Lg, Rg, idx = PR.seam_pairs(df, lab["last_glyph"], lab["first_glyph"])
    La, Ra, _ = PR.seam_pairs(df, lab["suffix"], lab["prefix"])
    return {
        "R1 prefix | first token, par-first line": (ctx1, lab["prefix"], rows, ["first|pf"], lid),
        "R1 prefix | first token, other line": (ctx1, lab["prefix"], rows, ["first|other"], lid),
        "R1 prefix | inner token, par-first line": (ctx1, lab["prefix"], rows, ["inner|pf"], lid),
        "first glyph | first token": (place, lab["first_glyph"], rows, ["first"], lid),
        "R2 suffix | last token": (place, lab["suffix"], rows, ["last"], lid),
        "last glyph | last token": (place, lab["last_glyph"], rows, ["last"], lid),
        "R3 first glyph | last glyph before": (Lg, Rg, idx, None, lid[idx]),
        "R3 prefix | last glyph before": (Lg, lab["prefix"][idx + 1], idx, None, lid[idx]),
        "seam prefix | suffix before": (La, Ra, idx, None, lid[idx]),
    }


def groupings(df: pd.DataFrame) -> dict:
    out = {"language": (df["language"].fillna("?").astype(str).to_numpy(), ["A", "B"])}
    hands = df["hand"].fillna("?").astype(str)
    out["hand"] = (hands.to_numpy(), sorted(h for h, n in hands.value_counts().items() if n >= MIN_HAND and h != "?"))
    sec = df["illustration"].fillna("?").astype(str)
    out["section"] = (sec.to_numpy(), sorted(s for s, n in sec.value_counts().items() if n >= MIN_SECTION and s != "?"))
    return out


def cmd_stability(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, cov_rows = [], []
    for tag in REAL_TAGS:
        ctx = ctx_of(tag)
        df = ctx.real
        for variant in ("a", "b"):
            lab = labels_for(df, ctx, variant)
            tables = stability_tables(df, lab)
            pos = position_class(df)
            for gname, (g, keep) in groupings(df).items():
                if len(keep) < 2:
                    continue
                pairs = [(a, b) for i, a in enumerate(keep) for b in keep[i + 1:]]
                for tname, (x, y, idx, contexts, units) in tables.items():
                    if variant == "b" and "prefix" not in tname and "suffix" not in tname:
                        continue
                    gg = g[idx]
                    st = PR.stability(x, y, gg, units, pairs=pairs, n_rep=N_REP, min_count=MIN_COUNT,
                                      contexts=contexts, seed=14)
                    st.insert(0, "table", tname)
                    st.insert(0, "grouping", gname)
                    st.insert(0, "variant", "-" if ("prefix" not in tname and "suffix" not in tname) else variant)
                    st.insert(0, "source", tag)
                    rows.append(st)
                    print(f"{tag} ({variant}) {gname} {tname}: {len(st)} pairs", flush=True)
                # per-group coverage (position of the four targets, seam glyph pairs)
                for grp in keep:
                    m = g == grp
                    sub = df[m]
                    for target in TARGETS:
                        if variant == "b" and target.endswith("glyph"):
                            continue
                        d = PR.item_decomposition(lab[target][m], pos[m], "item", "pos")
                        c = d["coverage"].set_index("items")
                        cov_rows.append({"source": tag, "variant": "-" if target.endswith("glyph") else variant,
                                         "grouping": gname, "group": grp, "n_tokens": int(m.sum()),
                                         "decomposition": "position", "target": target, "mi": d["mi"],
                                         "n_80_items": c.loc["item", "n_80"], "n_80_pos": c.loc["pos", "n_80"],
                                         "n_80_pairs": c.loc["pair", "n_80"],
                                         "top_items": " ".join(f"{r['item'] or '∅'}:{r['share']:.2f}" for _, r in d["by_x"].head(5).iterrows())})
                    if variant == "a":
                        L, R, _ = PR.seam_pairs(sub, lab["last_glyph"][m], lab["first_glyph"][m])
                        d = PR.item_decomposition(L, R, "last", "first")
                        c = d["coverage"].set_index("items")
                        cov_rows.append({"source": tag, "variant": "-", "grouping": gname, "group": grp,
                                         "n_tokens": int(m.sum()), "decomposition": "seam_glyph",
                                         "target": "last_glyph->first_glyph", "mi": d["mi"],
                                         "n_80_items": c.loc["last", "n_80"], "n_80_pos": c.loc["first", "n_80"],
                                         "n_80_pairs": c.loc["pair", "n_80"],
                                         "top_items": " ".join(f"{r['last']}>{r['first']}:{r['contribution'] / d['mi']:.2f}"
                                                              for _, r in d["pairs"].head(5).iterrows())})
    pd.concat(rows, ignore_index=True).to_csv(OUT / "stability.csv", index=False, float_format="%.6g")
    pd.DataFrame(cov_rows).to_csv(OUT / "coverage_by_group.csv", index=False, float_format="%.6g")
    print("stability.csv, coverage_by_group.csv written")


# --------------------------------------------------------------------------- #
# step 3: augmented generator
# --------------------------------------------------------------------------- #


def rules_for(ctx) -> PR.PlaceRules:
    pre, suf = lists_of(ctx, "b")
    return PR.learn_place_rules(ctx.real, pre, suf, preset="extended", min_context=MIN_COUNT)


def calib_metrics(df: pd.DataFrame) -> dict:
    """Calibration targets of 004 that are not in the profiles: H2 of glyphs, mean length, vocabulary."""
    gl = glyph_lists(df, "extended")
    vc = df["token"].value_counts()
    return {"calib__H2_glyph": conditional_entropy(gl, order=1), "calib__mean_len": float(np.mean([len(g) for g in gl])),
            "calib__n_types": int(len(vc)), "calib__hapax_share": float((vc == 1).mean())}


def measure_all(df: pd.DataFrame, ctx, seed: int) -> dict:
    out = {}
    out.update({f"p009__{k}": v for k, v in M9.profile(df, ctx, seed, "sub").items()})
    out.update({f"p010__{k}": v for k, v in M10.line_profile(df, ctx, seed, "sub").items()})
    out.update({f"p011__{k}": v for k, v in M11.measure(df, ctx, seed, "null").items()})
    out.update(calib_metrics(df))
    return out


def cmd_augment(args) -> None:
    tag, variant = args.tag, args.variant
    which = VARIANTS[variant]
    ctx = ctx_of(tag)
    rules = rules_for(ctx)
    source = f"ts2020+{variant}:{tag}"

    def one(seed: int) -> dict:
        gen = ctx.corpus("ts2020", seed)
        aug = PR.apply_place_rules(gen, rules, which, seed=seed)
        row = {"source": source, "kind": "null", "n_changed": sum(aug.attrs["place_rules_changed"].values()),
               **{f"n_changed_{k}": v for k, v in aug.attrs["place_rules_changed"].items()}}
        row.update(measure_all(aug, ctx, seed))
        return row

    seeds = M9.parse_seeds(args.seeds, list(range(N_NULL)))
    params = {"source": source, "rules": list(which), "rule_lists": {"prefixes": rules.prefixes, "suffixes": rules.suffixes},
              "min_context": MIN_COUNT, "ts2020_config": M9.TS2020_CONFIG, "n_null": N_NULL,
              "profile_prefixes": ctx.prefixes, "profile_suffixes": ctx.suffixes}
    run_seeds("014_augment", seeds, one, OUT / "runs" / f"{variant}__{tag}", params=params, progress=M9.progress)


def cmd_rules(args) -> None:
    """Write the learned tables R1-R3 of ZL split / join (``rules_tables.csv``)."""
    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for tag in AUG_TAGS:
        r = rules_for(ctx_of(tag))
        for rule in ("R1", "R2", "R3"):
            frames.append(r.table(rule).assign(source=tag))
    out = pd.concat(frames, ignore_index=True)
    out["item"] = out["item"].where(out["item"] != "", EMPTY)
    out.to_csv(OUT / "rules_tables.csv", index=False, float_format="%.5g")
    print(f"rules_tables.csv: {len(out)} rows")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("inventory")
    sub.add_parser("stability")
    sub.add_parser("rules")
    p = sub.add_parser("augment")
    p.add_argument("--variant", required=True, choices=list(VARIANTS))
    p.add_argument("--tag", required=True, choices=AUG_TAGS)
    p.add_argument("--seeds", default=None)
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "summarize":
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from summarize import summarize  # noqa: E402  (experiments/014_place_rules_inventory/summarize.py)

        summarize()
        return
    {"inventory": cmd_inventory, "stability": cmd_stability, "augment": cmd_augment, "rules": cmd_rules}[args.cmd](args)


if __name__ == "__main__":
    main()
