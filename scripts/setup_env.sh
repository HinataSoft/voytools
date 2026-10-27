#!/bin/bash
# Setup script for a Claude Code cloud environment (paste into the environment's
# "Setup script" field, or run once in a fresh session).  Finishes in well under
# five minutes so the environment cache can be built.
set -e
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
pip install -q -e ".[dev,plot]"
python -c "import voytools, numpy, pandas; print('voytools', voytools.__version__, '| numpy', numpy.__version__, '| pandas', pandas.__version__)"
python -m pytest -q -x --no-header 2>&1 | tail -3
