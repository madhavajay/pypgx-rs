#!/usr/bin/env bash
# Characterize the upstream PyPGx (Python) test suite as a CI sanity check.
#
# On the vendored v0.26.0 data, 3 of the 6 tests pass and 3 are
# data-consistency assertions that the shipped tables themselves violate
# (duplicate allele 19-39738787-C-T; the ACYP2 variant-table diffs;
# the MT-RNR1 priority mismatch). The Rust parity suite asserts the reference's
# *computed* values, reproducing those exact discrepancies — so `cargo test`
# is fully green. This script passes iff the reference is in that documented
# state: the 3 self-consistent tests pass and the 3 data tests fail.
set -uo pipefail
cd "$(dirname "$0")/../repos/pypgx"

echo "== Self-consistent PyPGx reference tests (must pass) =="
python -m unittest -v \
  test.TestPypgx.test_diplotype_table \
  test.TestPypgx.test_equation_table \
  test.TestPypgx.test_predict_alleles
pass_rc=$?

echo ""
echo "== Data-consistency assertions the shipped v0.26.0 tables violate (expected to FAIL) =="
if python -m unittest \
  test.TestPypgx.test_allele_table \
  test.TestPypgx.test_definition_table \
  test.TestPypgx.test_priority_table >/dev/null 2>&1; then
  echo "UNEXPECTED: the data-consistency tests passed — upstream data changed; revisit parity notes."
  exit 1
fi
echo "Confirmed: the 3 documented data-consistency failures are present"
echo "(the Rust suite reproduces these exact discrepancies)."

exit "$pass_rc"
