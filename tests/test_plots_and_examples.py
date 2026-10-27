"""Optional plots and the example scripts (the latter run on the real data and are skipped without it)."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from voytools import stats as S
from voytools.morphology import bpe_dependence_curve
from tests.conftest import make_tokens

ROOT = Path(__file__).resolve().parents[1]


def _env() -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("VOYTOOLS_DATA", None)
    return env


def _run(script: str, *args: str, cwd: Path = ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(ROOT / "examples" / script), *args],
        capture_output=True, text=True, cwd=str(cwd), timeout=300, env=_env(),
    )


# ------------------------------------------------------------------ plots


def test_plots_return_figures_and_save_files(tmp_path):
    pytest.importorskip("matplotlib")
    from voytools.plots import have_matplotlib, plot_bpe_curve, plot_positional, plot_zipf

    assert have_matplotlib()
    lines = [["qokeedy", "chedy", "shedy", "daiin", "ol", "or", "chol"] * 2 for _ in range(30)]
    df = make_tokens(lines)
    f1 = plot_zipf({"a": S.token_frequency(df)}, path=tmp_path / "zipf.png")
    f2 = plot_positional(S.positional_profile(df, which="first"), path=tmp_path / "pos.png")
    f3 = plot_bpe_curve(bpe_dependence_curve(lines, (0, 2, 4), preset="basic"), path=tmp_path / "bpe.png")
    for f, name in ((f1, "zipf.png"), (f2, "pos.png"), (f3, "bpe.png")):
        assert hasattr(f, "savefig") and (tmp_path / name).stat().st_size > 1000
    plot_zipf(S.token_frequency(df))  # a plain frame and no path are fine


# ------------------------------------------------------------------ examples


def test_example_reports_searched_paths_when_data_is_missing(tmp_path):
    r = _run("01_summary.py", cwd=tmp_path)
    if r.returncode == 0:
        pytest.skip("ZL3b-n.txt is found by the default search on this machine")
    assert "Tried" in r.stderr + r.stdout and "ZL3b-n.txt" in r.stderr + r.stdout


def test_example_rejects_a_missing_explicit_path(tmp_path):
    r = _run("01_summary.py", str(tmp_path / "nope.txt"))
    assert r.returncode != 0 and "not found" in (r.stderr + r.stdout).lower()


@pytest.mark.parametrize("script, expect", [
    ("01_summary.py", "Glyph entropy of the space-less text"),
    ("02_compare_sections.py", "Labels vs paragraph text"),
    ("03_morphology_and_nulls.py", "BPE dependence gap"),
])
def test_examples_run_end_to_end(zl_path, script, expect):
    args = [str(zl_path)] + (["--shuffles", "3"] if script.startswith("03") else [])
    r = _run(script, *args)
    assert r.returncode == 0, r.stderr[-2000:]
    assert expect in r.stdout and "Traceback" not in r.stderr
    r.stdout.encode("ascii")  # console friendly
