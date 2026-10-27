"""Checkpointed, resumable runners for stochastic experiments.

Cloud and local sessions can be interrupted (command time limits, paused or
reclaimed VMs).  The helpers here make a many-run experiment safe to interrupt
and restart:

* :func:`run_seeds` calls ``fn(seed)`` for every seed in a list, appends each
  result as one row to ``<out_dir>/runs.csv`` immediately, and on restart skips
  the seeds that are already in the file.  A ``meta.json`` records the
  experiment name, parameters and timing.
* :func:`load_runs` reads the rows back; :func:`summarize_runs` aggregates the
  numeric columns (mean, sd, quantiles, n).
* :func:`seed_range` splits a global seed range into chunks so that several
  sessions can work on disjoint seed sets and the files can be merged later with
  :func:`merge_runs`.

Rows are written with :mod:`csv` and flushed on every write, so a run that is
killed mid-way loses at most the row being computed.  Results should be
committed to git regularly (see ``experiments/README.md``).
"""

from __future__ import annotations

import csv
import json
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

__all__ = [
    "RUNS_FILE",
    "META_FILE",
    "run_seeds",
    "load_runs",
    "summarize_runs",
    "merge_runs",
    "seed_range",
    "write_meta",
]

RUNS_FILE = "runs.csv"
META_FILE = "meta.json"


def seed_range(start: int, stop: int, chunks: int) -> list[list[int]]:
    """Split ``range(start, stop)`` into ``chunks`` nearly equal lists of seeds.

    Use it to give parallel sessions disjoint seed sets, e.g.
    ``seed_range(0, 1000, 5)[k]`` for session ``k``.
    """
    if chunks < 1:
        raise ValueError("chunks must be >= 1")
    seeds = list(range(start, stop))
    size, rem = divmod(len(seeds), chunks)
    out: list[list[int]] = []
    pos = 0
    for k in range(chunks):
        n = size + (1 if k < rem else 0)
        out.append(seeds[pos : pos + n])
        pos += n
    return out


def write_meta(out_dir: str | Path, name: str, params: Mapping[str, Any] | None = None, **extra: Any) -> Path:
    """Write (or update) ``meta.json`` in ``out_dir``; returns its path."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / META_FILE
    meta: dict[str, Any] = {}
    if path.exists():
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}
    meta["name"] = name
    if params is not None:  # a resumed run without params keeps the stored ones
        meta["params"] = dict(params)
    meta.setdefault("params", {})
    meta.update(extra)
    meta.setdefault("created", time.strftime("%Y-%m-%d %H:%M:%S"))
    meta["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
    path.write_text(json.dumps(meta, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    return path


def _done_seeds(path: Path) -> set[int]:
    if not path.exists():
        return set()
    done: set[int] = set()
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                done.add(int(row["seed"]))
            except (KeyError, ValueError):
                continue
    return done


def run_seeds(
    name: str,
    seeds: Iterable[int],
    fn: Callable[[int], Mapping[str, Any]],
    out_dir: str | Path,
    *,
    params: Mapping[str, Any] | None = None,
    resume: bool = True,
    progress: Callable[[int, int, int], None] | None = None,
) -> pd.DataFrame:
    """Run ``fn(seed)`` for every seed, appending one CSV row per run.

    Parameters
    ----------
    name:
        Experiment name (stored in ``meta.json``).
    seeds:
        Seeds to run.  With ``resume=True`` the seeds already present in
        ``runs.csv`` are skipped.
    fn:
        Callable returning a flat mapping of column -> scalar for one run.  It
        must be deterministic given the seed.  The ``seed`` and ``elapsed_s``
        columns are added automatically.
    out_dir:
        Directory for ``runs.csv`` and ``meta.json`` (created if needed).
    params:
        Parameters to record in ``meta.json``.
    progress:
        Optional callback ``progress(done, total, seed)`` called after each run.

    Returns
    -------
    pandas.DataFrame
        All rows now in ``runs.csv`` (including earlier ones).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    runs_path = out / RUNS_FILE
    seeds = list(seeds)
    done = _done_seeds(runs_path) if resume else set()
    todo = [s for s in seeds if s not in done]
    write_meta(out, name, params, n_seeds_requested=len(seeds), n_seeds_done_before=len(done))

    if not resume and runs_path.exists():
        runs_path.unlink()

    header_written = runs_path.exists() and runs_path.stat().st_size > 0
    fieldnames: list[str] | None = None
    if header_written:
        with runs_path.open("r", encoding="utf-8", newline="") as fh:
            fieldnames = next(csv.reader(fh), None)

    total = len(todo)
    for i, seed in enumerate(todo, 1):
        t0 = time.perf_counter()
        result = dict(fn(seed))
        row = {"seed": seed, **result, "elapsed_s": round(time.perf_counter() - t0, 4)}
        if fieldnames is None:
            fieldnames = list(row.keys())
        missing = [k for k in row if k not in fieldnames]
        if missing:
            raise ValueError(
                f"run for seed {seed} produced new columns {missing}; "
                f"existing runs.csv has columns {fieldnames}"
            )
        with runs_path.open("a", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            if not header_written:
                writer.writeheader()
                header_written = True
            writer.writerow({k: row.get(k, "") for k in fieldnames})
            fh.flush()
        if progress is not None:
            progress(i, total, seed)

    write_meta(out, name, params, n_seeds_requested=len(seeds), n_seeds_done=len(_done_seeds(runs_path)))
    return load_runs(out)


def load_runs(out_dir: str | Path) -> pd.DataFrame:
    """Read ``runs.csv`` from ``out_dir`` (empty DataFrame if absent)."""
    path = Path(out_dir) / RUNS_FILE
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path)


def summarize_runs(runs: pd.DataFrame, quantiles: Sequence[float] = (0.025, 0.5, 0.975)) -> pd.DataFrame:
    """Aggregate the numeric columns of a runs table: n, mean, sd, min, quantiles, max."""
    if runs.empty:
        return pd.DataFrame()
    num = runs.select_dtypes("number").drop(columns=[c for c in ("seed",) if c in runs.columns])
    rows = {}
    for col in num.columns:
        s = num[col].dropna()
        row = {"n": int(s.shape[0]), "mean": s.mean(), "sd": s.std(ddof=1) if s.shape[0] > 1 else float("nan")}
        row["min"] = s.min()
        for q in quantiles:
            row[f"q{q:g}"] = s.quantile(q)
        row["max"] = s.max()
        rows[col] = row
    return pd.DataFrame.from_dict(rows, orient="index")


def merge_runs(dirs: Iterable[str | Path], out_dir: str | Path | None = None) -> pd.DataFrame:
    """Concatenate ``runs.csv`` files from several directories (e.g. parallel sessions).

    Duplicate seeds keep the first occurrence.  When ``out_dir`` is given the
    merged table is written there as ``runs.csv``.
    """
    frames = [load_runs(d) for d in dirs]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame()
    merged = pd.concat(frames, ignore_index=True)
    if "seed" in merged.columns:
        merged = merged.drop_duplicates(subset="seed", keep="first").sort_values("seed").reset_index(drop=True)
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out / RUNS_FILE, index=False)
    return merged
