"""Command line interface, exercised through ``main([...])`` on the synthetic files."""

from __future__ import annotations

import io

import pandas as pd
import pytest

from voytools.cli import build_parser, main
from tests.conftest import SYNTH_INTERLINEAR


@pytest.fixture()
def interlinear_path(tmp_path):
    p = tmp_path / "interlinear.txt"
    p.write_text(SYNTH_INTERLINEAR, encoding="utf-8")
    return p


def run(capsys, *argv):
    code = main([str(a) for a in argv])
    out = capsys.readouterr()
    return code, out.out, out.err


def test_summary(capsys, synth_path):
    code, out, err = run(capsys, "summary", synth_path)
    assert code == 0 and err == ""
    assert "Pages: 2" in out and "Loci: 12" in out and "Tokens in loci: 34" in out
    assert "Alphabet: Eva-" in out and "Vocabulary" in out and "Glyph entropy" in out
    assert "Loci by generic type" in out
    out.encode("ascii")  # console friendly: no non-ASCII characters


def test_summary_renders_high_ascii_as_codes(capsys, synth_path, tmp_path):
    p = tmp_path / "hi.txt"
    p.write_text("#=IVTFF Eva- 2.0 M\n<f1r>\n<f1r.1,@P0>   fa@150;ar.sory\n", encoding="utf-8")
    code, out, _ = run(capsys, "freq", p, "--no-clean")
    assert code == 0 and "@150;" in out
    out.encode("ascii")


def test_freq_and_csv_output(capsys, synth_path):
    code, out, _ = run(capsys, "freq", synth_path, "--top", "3", "--csv")
    assert code == 0
    body = [l for l in out.splitlines() if l and not l.startswith("#")]
    frame = pd.read_csv(io.StringIO("\n".join(body)))
    assert list(frame.columns) == ["token", "count", "freq", "rank"] and len(frame) == 3
    assert frame.loc[0, "token"] == "kor" and frame.loc[0, "count"] == 2


def test_freq_by_language(capsys, synth_path):
    code, out, _ = run(capsys, "freq", synth_path, "--by", "language", "--top", "2", "--generic-type", "all")
    assert code == 0 and "language" in out


def test_entropy_both_presets_and_options(capsys, synth_path):
    code, out, _ = run(capsys, "entropy", synth_path)
    assert code == 0 and "basic" in out and "extended" in out and "H2" in out
    code, out, _ = run(capsys, "entropy", synth_path, "--preset", "none", "--scope", "all")
    assert code == 0 and "none" in out
    code, out, _ = run(capsys, "entropy", synth_path, "--boundary", "_", "--preset", "basic")
    assert code == 0 and "basic" in out


def test_positions_affixes_and_mi(capsys, synth_path):
    code, out, _ = run(capsys, "positions", synth_path, "--which", "both", "--top", "2", "--permutations", "5")
    assert code == 0 and "Position classes (first" in out and "Position classes (last" in out
    assert "Permutation p-values" in out
    code, out, _ = run(capsys, "affixes", synth_path, "--min-count", "1", "--top", "3")
    assert code == 0 and "Prefixes" in out and "Suffixes" in out
    code, out, _ = run(capsys, "affixes", synth_path, "--min-count", "1", "--side", "suffix")
    assert code == 0 and "Prefixes" not in out and "Suffixes" in out
    code, out, _ = run(capsys, "mi", synth_path, "--shuffles", "3")
    assert code == 0 and "Adjacent-token MI" in out and "Edge-glyph MI" in out and "Neighbour similarity" in out
    code, out, _ = run(capsys, "mi", synth_path, "--shuffles", "2", "--no-neighbor")
    assert code == 0 and "Neighbour similarity" not in out


def test_export_csv_with_tables(capsys, synth_path, tmp_path):
    out_path = tmp_path / "tokens.csv"
    code, out, _ = run(capsys, "export", synth_path, out_path, "--tables")
    assert code == 0 and "wrote 34 tokens" in out
    df = pd.read_csv(out_path, keep_default_na=False)
    assert len(df) == 34 and "token" in df.columns and "folio" in df.columns
    assert (tmp_path / "tokens_loci.csv").exists() and (tmp_path / "tokens_pages.csv").exists()
    # high-ASCII characters are written as codes, not as private-use characters
    assert "@150;" in "".join(df["token_raw"].astype(str))


def test_token_selection_options(capsys, synth_path):
    code, out, _ = run(capsys, "summary", synth_path, "--language", "B", "--generic-type", "P")
    assert code == 0 and "language=B" in out
    code, out, err = run(capsys, "summary", synth_path, "--language", "Q")
    assert code == 0 and "Selected tokens" in out and "no tokens left" in err


def test_interlinear_picks_most_complete_transcriber(capsys, interlinear_path):
    code, out, err = run(capsys, "freq", interlinear_path, "--top", "2")
    assert code == 0 and "using 'H'" in err
    code, out, err = run(capsys, "freq", interlinear_path, "--transcriber", "C", "--top", "2")
    assert code == 0 and err == ""


def test_interlinear_export_keeps_all_transcribers(capsys, interlinear_path, tmp_path):
    p = tmp_path / "all.csv"
    code, out, _ = run(capsys, "export", interlinear_path, p)
    assert code == 0 and "wrote 20 tokens" in out
    assert set(pd.read_csv(p, keep_default_na=False)["transcriber"]) == {"H", "C", "Z"}  # U has only filler characters, hence no tokens


def test_errors_give_exit_code_2(capsys, tmp_path):
    code, _, err = run(capsys, "summary", tmp_path / "does_not_exist.txt")
    assert code == 2 and "error" in err
    bad = tmp_path / "bad.txt"
    bad.write_text("#=IVTFF Eva- 2.0 M\n<f1r>\n<f1r.1 @P0> x\n", encoding="utf-8")
    code, _, err = run(capsys, "summary", bad)
    assert code == 2 and "error" in err


def test_parser_and_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0 and "voytools 0.1.0" in capsys.readouterr().out
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 2
    ns = build_parser().parse_args(["freq", "x.txt", "--no-clean"])
    assert ns.clean is False and ns.space_policy == "split" and ns.alternatives == "first"


def test_python_dash_m_entry_point(synth_path):
    import os
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    r = subprocess.run(
        [sys.executable, "-m", "voytools", "freq", str(synth_path), "--top", "2"],
        capture_output=True, text=True, cwd=str(root), timeout=120,
        env={**os.environ, "PYTHONPATH": str(root) + os.pathsep + os.environ.get("PYTHONPATH", "")},
    )
    assert r.returncode == 0 and "kor" in r.stdout
