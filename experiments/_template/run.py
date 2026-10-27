"""Template experiment: one checkpointed run per seed, real text and null models side by side.

Usage (from the repository root)::

    python experiments/_template/run.py --seeds 0:20 --out results/_template
    python experiments/_template/run.py --seeds 20:40 --out results/_template/part1   # another session

Each seed produces one row with the same statistics for the real text (``source=real``
is seed-independent and therefore computed once, seed 0) and for the null models.
Copy this file into ``experiments/NNN_name/run.py`` and replace :func:`one_run`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))  # make voytools importable without pip install

from voytools import filter_tokens, load_tables, run_seeds, seed_range, summarize_runs  # noqa: E402
from voytools.nulls import autocopy_like, shuffle_within_lines  # noqa: E402
from voytools.stats import edge_glyph_mi, neighbor_similarity  # noqa: E402

DATA = ROOT / "data"


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


def load_text(name: str, space_policy: str):
    tables = load_tables(DATA / name, space_policy=space_policy)
    return filter_tokens(tables.tokens, generic_type="P", clean=True)


def stats_for(df, seed: int) -> dict:
    """The statistics compared across sources; returns a flat dict."""
    edge = edge_glyph_mi(df, shuffles=10, seed=seed).iloc[0]
    near = neighbor_similarity(df, window=5, shuffles=5, seed=seed).iloc[0]
    return {
        "edge_mi_excess": float(edge["excess_bits"]),
        "near_share": float(near["share_near"]),
        "near_ratio": float(near["ratio"]),
    }


def make_runner(real, policy: str):
    def one_run(seed: int) -> dict:
        out = {"space_policy": policy}
        for source, df in (
            ("real", real),
            ("shuffle_lines", shuffle_within_lines(real, seed=seed)),
            ("autocopy", autocopy_like(real, seed=seed, window=8, p_copy_exact=0.35)),
        ):
            for k, v in stats_for(df, seed).items():
                out[f"{source}__{k}"] = v
        return out

    return one_run


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", default="ZL3b-n.txt")
    ap.add_argument("--space-policy", default="split", choices=["split", "join", "strict"])
    ap.add_argument("--seeds", default="0:10", help="'a:b', 'a:b/chunks:k' or comma list")
    ap.add_argument("--out", default=str(ROOT / "results" / "_template"))
    args = ap.parse_args()

    real = load_text(args.file, args.space_policy)
    seeds = parse_seeds(args.seeds)
    runs = run_seeds(
        "template",
        seeds,
        make_runner(real, args.space_policy),
        args.out,
        params={"file": args.file, "space_policy": args.space_policy, "seeds": args.seeds},
        progress=lambda i, n, s: print(f"  {i}/{n} seed={s}", flush=True),
    )
    summary = summarize_runs(runs)
    summary.to_csv(Path(args.out) / "summary.csv")
    print(summary.to_string(float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
