# pypgx-rs

A 1-for-1 Rust port of [**PyPGx**](https://github.com/sbslee/pypgx), a package
for pharmacogenomics (PGx) research. The upstream Python package is vendored as a
git submodule at [`repos/pypgx/`](./repos/pypgx) and used as the reference for byte-for-byte
parity.

See [`TODO.md`](./TODO.md) for the full port plan and status.

## Status

The **pure analytical port is complete and verified**: reference tables, the
`fuc` primitives, `Archive` read/write, the full allele/phenotype/score logic,
and the whole `predict_alleles → call_genotypes → call_phenotypes →
combine_results → count_alleles`/`compare_genotypes` pipeline — all reproducing
the Python reference's computed output exactly. Everything that needs an external
program (Beagle, samtools, bcftools), scikit-learn (the CNV model), or matplotlib
(plots) is deferred by design with a documented stub in `src/external.rs`; see
TODO.md.

## Layout

| Module | Mirrors | Purpose |
|---|---|---|
| `src/fuc.rs` | the used slice of `fuc` | `parse_variant`, `sort_variants`, `VcfFrame` |
| `src/core.rs` | `pypgx/api/core.py` | reference tables + allele/phenotype/score logic |
| `src/sdk.rs` | `pypgx/sdk` | `Archive`, semantic types, `PgxError` |
| `src/api.rs` | `pypgx/api/utils.py` | `predict_alleles`, `call_phenotypes`, `combine_results`, `compare_genotypes`, `count_alleles` |
| `src/genotype.rs` | `pypgx/api/genotype.py` | `call_genotypes` (+ all SV genotypers) |
| `src/external.rs` | external/sklearn/plot fns | deferred stubs (`PgxError::NotPorted`) |
| `src/table.rs` | — | pandas-like frame with exact NaN/bool semantics |
| `src/bin/pypgx.rs` | `pypgx/cli` + `__main__.py` | clap CLI (11 subcommands) |

## Build & test

```sh
cargo test          # 24 tests: parity + breadth + round-trip + pipeline
cargo run -- predict-alleles tests/fixtures/CYP4F2-GRCh37.zip
```

## Differential parity vs Python

Parity is enforced against ground truth dumped from the Python reference. To
re-verify against upstream:

```sh
# one-time: reference env (Python 3.10 + pypgx 0.26.0 + fuc + scikit-learn)
uv venv --python 3.10 .refenv
source .refenv/bin/activate && uv pip install ./pypgx

# regenerate ground truth from Python, then run the Rust suite against it
./tools/diff_harness.sh
```

> Note: on the vendored v0.26.0 data, 3 of upstream's 6 `test.py` checks are
> data-consistency assertions the shipped tables violate (a duplicate allele,
> two variant-table diffs, and a phenotype-table mismatch). The Rust tests assert
> the reference's *computed* values, so they reproduce these discrepancies
> exactly rather than masking them.
