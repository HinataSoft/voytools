"""Tests for the checkpointed experiment runner."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from voytools.experiment import load_runs, merge_runs, run_seeds, seed_range, summarize_runs


def _fn(seed: int) -> dict:
    return {"value": seed * 2.0, "flag": seed % 2}


def test_run_and_resume(tmp_path):
    out = tmp_path / "exp"
    df = run_seeds("demo", range(5), _fn, out, params={"a": 1})
    assert list(df["seed"]) == [0, 1, 2, 3, 4]
    assert "elapsed_s" in df.columns and "value" in df.columns

    calls: list[int] = []

    def fn2(seed: int) -> dict:
        calls.append(seed)
        return _fn(seed)

    df2 = run_seeds("demo", range(8), fn2, out)
    assert calls == [5, 6, 7]  # earlier seeds were skipped
    assert len(df2) == 8
    meta = json.loads((out / "meta.json").read_text(encoding="utf-8"))
    assert meta["name"] == "demo" and meta["params"] == {"a": 1} and meta["n_seeds_done"] == 8


def test_no_resume_restarts(tmp_path):
    out = tmp_path / "exp"
    run_seeds("demo", range(3), _fn, out)
    df = run_seeds("demo", range(2), _fn, out, resume=False)
    assert list(df["seed"]) == [0, 1]


def test_new_columns_rejected(tmp_path):
    out = tmp_path / "exp"
    run_seeds("demo", range(2), _fn, out)
    with pytest.raises(ValueError):
        run_seeds("demo", range(2, 3), lambda s: {"other": 1}, out)


def test_summary_and_merge(tmp_path):
    a = run_seeds("demo", seed_range(0, 10, 2)[0], _fn, tmp_path / "a")
    b = run_seeds("demo", seed_range(0, 10, 2)[1], _fn, tmp_path / "b")
    assert len(a) == 5 and len(b) == 5
    merged = merge_runs([tmp_path / "a", tmp_path / "b", tmp_path / "a"], tmp_path / "m")
    assert list(merged["seed"]) == list(range(10))
    assert load_runs(tmp_path / "m").shape[0] == 10
    s = summarize_runs(merged)
    assert s.loc["value", "n"] == 10
    assert s.loc["value", "mean"] == pytest.approx(9.0)
    assert {"q0.025", "q0.5", "q0.975", "sd", "min", "max"} <= set(s.columns)


def test_seed_range_partition():
    parts = seed_range(0, 7, 3)
    assert parts == [[0, 1, 2], [3, 4], [5, 6]]
    assert sum(len(p) for p in seed_range(0, 1000, 7)) == 1000


def test_empty_helpers(tmp_path):
    assert load_runs(tmp_path).empty
    assert summarize_runs(pd.DataFrame()).empty
    assert merge_runs([tmp_path]).empty
