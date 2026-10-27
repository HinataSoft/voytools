"""Helpers shared by the examples: locate the transliteration file and make ``voytools`` importable.

Search order for the data file (first existing path wins):

1. the path given on the command line;
2. ``<repo>/data/<name>`` (the copy committed to the git repository);
3. ``<examples>/../../transliterations/<name>``  (layout ``D:\\Work\\Voynich\\analysis\\examples`` +
   ``D:\\Work\\Voynich\\transliterations``);
4. ``<examples>/../transliterations/<name>``;
5. ``./transliterations/<name>`` and ``../transliterations/<name>`` relative to the working directory;
6. the folder named by the environment variable ``VOYTOOLS_DATA``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_NAME = "ZL3b-n.txt"

try:  # the examples also work from a source checkout that has not been pip-installed
    import voytools  # noqa: F401
except ImportError:
    sys.path.insert(0, str(ROOT))


def candidates(name: str = DEFAULT_NAME) -> list[Path]:
    """Default locations tried when no path is given."""
    out = [
        HERE.parent / "data" / name,  # copy inside the git repository (analysis/data)
        HERE.parent.parent / "transliterations" / name,
        HERE.parent / "transliterations" / name,
        Path.cwd() / "transliterations" / name,
        Path.cwd().parent / "transliterations" / name,
    ]
    env = os.environ.get("VOYTOOLS_DATA")
    if env:
        out.append(Path(env) / name)
    return out


def resolve_data(arg: str | None = None, name: str = DEFAULT_NAME) -> Path:
    """Return the transliteration file to use or exit with a message listing what was tried."""
    if arg:
        p = Path(arg)
        if p.is_dir():
            p = p / name
        if p.exists():
            return p
        sys.exit(f"File not found: {p}")
    tried = candidates(name)
    for p in tried:
        if p.exists():
            return p
    listing = "\n  ".join(str(p) for p in tried)
    sys.exit(
        f"Could not find {name}. Pass the path as the first argument, or set VOYTOOLS_DATA.\nTried:\n  {listing}"
    )


def show(title: str, obj, digits: int = 4) -> None:
    """Print a titled DataFrame/Series with a fixed number of decimals."""
    import pandas as pd

    print(f"\n== {title} ==")
    if isinstance(obj, (pd.DataFrame, pd.Series)):
        print(obj.to_string(float_format=lambda x: f"{x:.{digits}f}"))
    else:
        print(obj)
