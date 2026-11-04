#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$REPO_ROOT/.venv"
REQ="$REPO_ROOT/requirements.txt"
PYTHON_BIN="${PYTHON:-python3}"

unset PYTHONPATH

if [[ ! -x "$VENV/bin/python" ]]; then
    "$PYTHON_BIN" -m venv "$VENV"
fi

"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -r "$REQ"
"$VENV/bin/python" -c "import numpy, nibabel, matplotlib, h5py"

echo "neuro-tools setup complete."
