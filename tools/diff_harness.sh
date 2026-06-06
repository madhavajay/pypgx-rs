#!/usr/bin/env bash
# Differential test harness: regenerate the ground-truth fixtures from the
# Python reference, then run the Rust parity suite against them. A green run
# means the Rust port reproduces PyPGx's computed outputs exactly.
#
# Requires the reference venv at ../.refenv (see TODO.md Phase 0):
#   uv venv --python 3.10 .refenv && source .refenv/bin/activate && uv pip install ./repos/pypgx
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [[ ! -d .refenv ]]; then
  echo "Reference env .refenv not found. Create it with:"
  echo "  uv venv --python 3.10 .refenv && source .refenv/bin/activate && uv pip install ./repos/pypgx"
  exit 1
fi

source .refenv/bin/activate

echo "==> Regenerating ground-truth fixtures from Python reference"
python tools/gen_truth.py  > tests/fixtures/truth.json
python tools/gen_truth2.py > tests/fixtures/truth2.json

echo "==> Running Rust parity suite against the fixtures"
cargo test --quiet

echo "==> Differential parity OK"
