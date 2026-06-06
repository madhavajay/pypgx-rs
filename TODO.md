# pypgx-rs — Port Plan & TODO

A 1-for-1 port of [**pypgx**](https://github.com/sbslee/pypgx) (Python) to Rust.

> **Reference source of truth:** the upstream Python package, vendored as a git
> submodule at [`./pypgx`](./pypgx). Pin it; we port against a fixed commit and
> bump deliberately. Every Rust function must reproduce its Python counterpart's
> output **byte-for-byte** where feasible (CSV/VCF text, table contents, allele
> calls).

---

## Status (2026-06-06)

**The pure analytical port is implemented, and its tested paths reproduce the
Python reference's computed outputs exactly; everything that needs an external
program / sklearn / matplotlib is deferred by design with a documented stub.**
Parity is proven by a differential harness (`tools/diff_harness.sh`) that
regenerates ground truth from Python and runs the Rust suite against it —
**24 tests green** (10 parity + 10 breadth + 2 round-trip + 2 pipeline-slice).

> ⚠️ **Verification caveat (audited 2026-06-06).** The whole suite keys off a
> single gene (CYP4F2) with trivial single-variant alleles. The code for the
> SV/CNV machinery (12 of 13 SV genotypers + every CNV branch in
> `call_genotypes`), the `sort_alleles`/`collapse_alleles` comparators, the
> multi-variant `predict_alleles` VariantData ordering, and the entire CLI is
> **present and code-faithful but exercised by no test**. Read "verified" below
> as "verified on CYP4F2" for that slice. Full list in [§12](#12-known-gaps-audit-2026-06-06).
>
> ⚠️ **Build note.** If the crate directory is moved/renamed, run `cargo clean`:
> stale test binaries bake in the old `CARGO_MANIFEST_DIR`, so the
> archive-reading tests fail with spurious `No such file or directory` while the
> `include_str!`-backed table tests still pass.

Done & verified:
- ✅ Cargo crate (`src/`), 9 data tables embedded via `include_str!`.
- ✅ `fuc` primitives: `parse_variant`, `sort_variants`, `VcfFrame`
  (`from_string`/`sort`/`samples`/`to_variants`/`get_af`/`Display`). Note: Python's
  `from_dict` is **not** ported as a named method — its one call site
  (`build_definition_table`) uses the equivalent `VcfFrame::new(meta, cols, rows)`.
- ✅ `core`: all 9 `load_*_table`, `list_*`, `has_*`, every `get_*` accessor,
  `build_definition_table`, `collapse_alleles`, `sort_alleles`,
  `predict_phenotype`, `predict_score`, `get_recommendation`.
- ✅ `sdk`: `Archive` read **and write** (ZIP), semantic types, `PgxError`.
- ✅ `api::predict_alleles` — byte-identical to Python on the CYP4F2 fixtures.
- ✅ `genotype::call_genotypes` — SimpleGenotyper + all 13 SV genotypers.
- ✅ `api`: `call_phenotypes`, `combine_results`, `compare_genotypes`,
  `count_alleles`. The full `predict → genotypes → phenotypes → results`
  pipeline runs end-to-end (verified vs Python: Genotype `*1/*2`, Phenotype
  `Indeterminate`, counts `*1×1 *2×1`).
- ✅ CLI: 11 subcommands covering the whole pure pipeline.

Deferred by design (`src/external.rs`, each returns `PgxError::NotPorted` with
the exact missing dependency + command): Beagle phasing, depth/CNV/BAM steps
(samtools/bcftools), the sklearn `Model[CNV]`, plots (matplotlib), and the
pipelines that orchestrate them. These cannot be implemented or verified without
those binaries; native ports are Phases 4/6/8/9.

> **Key finding (faithful parity):** upstream `test.py` does **not** fully pass
> on the vendored v0.26.0 data — 3 of 6 are data-consistency assertions that the
> shipped tables violate (duplicate allele `19-39738787-C-T`, the `ACYP2` /
> `CYP17A1` variant-table diffs, and the `MT-RNR1` priority mismatch). The Rust
> parity tests assert the reference's *computed* values, so they reproduce these
> exact discrepancies rather than papering over them.

---

## 1. Guiding principles

1. **Parity first, idiomatic second.** Mirror the Python module/function layout
   and behavior exactly before refactoring into more idiomatic Rust. Keep the
   same function names, argument names, defaults, and return semantics.
2. **Defer external tooling.** Initially we **shell out** to the same external
   programs pypgx uses (Beagle, Java, samtools, bcftools). Port them natively
   later. See [§5](#5-external-tools-deferred).
3. **Differential testing is the spec.** For every behavior, run Python pypgx and
   Rust pypgx on identical inputs and diff outputs. The Python `test.py` suite is
   the first acceptance gate; see [§7](#7-test-parity).
4. **Tables are embedded data.** The 9 reference CSVs ship inside the binary
   (`include_str!`), exactly as Python ships them via `package_data`.

---

## 2. What we are porting (source inventory)

Upstream is ~5k lines of Python. Module map:

| Python | Lines | Responsibility |
|---|---|---|
| `pypgx/api/core.py` | 1598 | Reference tables, allele/phenotype/score logic. **Pure, no I/O deps.** Port first. |
| `pypgx/api/utils.py` | 1534 | `predict_alleles`, Beagle phasing, depth/coverage, CNV import, VCF consolidation. Shells out to externals. |
| `pypgx/api/genotype.py` | 677 | `call_genotypes` (diplotype assembly). |
| `pypgx/api/pipeline.py` | 316 | `run_ngs_pipeline`, `run_chip_pipeline`, `run_long_read_pipeline`. |
| `pypgx/api/plot.py` | 494 | 5 plotting functions (matplotlib). **Lowest priority.** |
| `pypgx/sdk/utils.py` | 386 | `Archive` (zip container + semantic types), exceptions. |
| `pypgx/cli/*.py` | 32 files | Thin argparse wrappers → one API function each. |
| `pypgx/__main__.py` | 36 | CLI dispatch. |
| `pypgx/api/data/*.csv` | 9 files | Reference tables (see §4). |
| `pypgx/api/beagle.*.jar` | — | Bundled Beagle phasing tool. |

Public API surface (re-exported from `pypgx/__init__.py`): **~70 functions** + the
`Archive` class. Full list in `__init__.py`; treat that as the export checklist.

---

## 3. Rust crate layout (as implemented)

A single crate with modules mirroring the Python package 1:1 (chosen over a
multi-crate workspace for fast, low-friction iteration; can be split later if a
crate boundary earns its keep). `src/lib.rs` re-exports the same names as
`pypgx/__init__.py` so call sites read 1:1.

```
pypgx-rs/
├── Cargo.toml
├── data/*.csv                 # 9 reference tables (synced from submodule, embedded)
├── src/
│   ├── lib.rs                 # re-exports the PURE public surface (core/api/genotype/sdk);
│   │                          #   the `external` stubs are NOT re-exported at the crate root
│   ├── table.rs               # pandas-like Frame with exact NaN semantics
│   ├── fuc.rs                 # ← the `fuc` slice: parse_variant, sort_variants, VcfFrame
│   ├── core.rs                # ← pypgx/api/core.py
│   ├── sdk.rs                 # ← pypgx/sdk: Archive, semantic types, PgxError
│   ├── api.rs                 # ← pypgx/api/utils.py (pure fns: predict_alleles, …)
│   ├── genotype.rs            # ← pypgx/api/genotype.py (call_genotypes + SV genotypers)
│   ├── external.rs            # deferred stubs (NotPorted): beagle/samtools/sklearn/plot/pipeline
│   └── bin/pypgx.rs           # ← pypgx/cli + __main__.py (clap, 11 subcommands)
├── tests/
│   ├── parity.rs              # replicas of the 6 test.py tests + supporting parity
│   ├── breadth.rs             # parity for the extra pure core functions
│   ├── pipeline_slice.rs      # end-to-end pure pipeline on CYP4F2
│   ├── roundtrip.rs           # Archive write→read round-trip (Rust↔Rust)
│   └── fixtures/              # CYP4F2 archives + truth*.json (ground truth)
├── tools/
│   ├── gen_truth.py           # dump reference ground truth (the 6 tests)
│   ├── gen_truth2.py          # dump reference ground truth (breadth functions)
│   └── diff_harness.sh        # regenerate truth + run Rust suite (differential test)
└── pypgx/                     # git submodule (reference Python source)
```

---

## 4. Reference data tables (embed verbatim)

Copy `pypgx/api/data/*.csv` into `data/` (repo root; the actual single-crate
layout — there is no `crates/pypgx-core/`) and load with `include_str!` + the
`csv` crate. **Do not hand-edit**; sync from the submodule.

| Table | Key columns | `load_*` fn |
|---|---|---|
| `allele-table.csv` | Gene, StarAllele, ActivityScore, Function, GRCh37Core, GRCh37Tag, GRCh38Core, GRCh38Tag, SV | `load_allele_table` |
| `gene-table.csv` | Gene, Target, Control, Paralog, Variants, SV, PhenotypeMethod, RefAllele, GRCh3x{Default,Region,ExonStarts,ExonEnds}, Strand | `load_gene_table` |
| `variant-table.csv` | Gene, GRCh3xName, rsID, Chromosome, GRCh3xPosition, GRCh3xAllele, Variant, Impact, GRCh3xSynonym | `load_variant_table` |
| `diplotype-table.csv` | Gene, Diplotype, Phenotype | `load_diplotype_table` |
| `equation-table.csv` | Gene, Phenotype, Equation | `load_equation_table` |
| `phenotype-table.csv` | Gene, Phenotype, Priority | `load_phenotype_table` |
| `cnv-table.csv` | Gene, Name | `load_cnv_table` |
| `cpic-table.csv` | Gene, Drug, RxNorm, ATC, Guideline, CPIC/PharmGKB/FDA levels, PMID | `load_cpic_table` |
| `recommendation-table.csv` | Drug, Gene1, Phenotype1, Gene2, Phenotype2, Recommendation | `load_recommendation_table` |

**93 genes total.** Critical parity concerns:
- **NaN/empty-cell semantics.** pandas reads empty cells as `NaN`; many code paths
  branch on `pd.isna(...)`. Model this with `Option<String>` and replicate every
  `isna` check.
- **Column order & row order** must be preserved (tests compare `.unique()` lists
  and rely on stable ordering).
- `equation-table.csv` holds **Python expressions** evaluated at runtime
  (`predict_score`). Need a tiny safe expression evaluator (see §8 risk).

---

## 5. External tools (deferred)

pypgx invokes these via `subprocess` / bundled jar. **Phase 1: wrap them with
`std::process::Command`, identical args.** Port natively later.

| Tool | Where | Strategy now | Port later? |
|---|---|---|---|
| **Beagle** `beagle.22Jul22.46e.jar` | `utils.estimate_phase_beagle` (`java -Xmx2g -jar ...`) | Bundle the jar; shell out via `java`. Replicate the EM-retry fallback (`em=true` then `em=false` on `CalledProcessError`). | **Yes, later** — native phasing is the big one. |
| **Java** | runs Beagle | Require `java` on PATH; surface clear error if missing. | With Beagle port. |
| **samtools** | BAM slicing / depth (`pybam`, `slice_bam`) | Shell out. | Later. |
| **bcftools** | VCF ops (via `fuc`/`pyvcf`) | Shell out where Python does. | Later. |

> ⚠️ Current dev box is missing `beagle`, `samtools`, `bcftools` (only `java`
> present). Provisioning these is a Phase 0 task for end-to-end tests.

---

## 6. The `fuc` dependency (the hidden iceberg)

pypgx leans on the [`fuc`](https://github.com/sbslee/fuc) bioinformatics library.
Imports used across the codebase:

```
from fuc import pyvcf, pycov, pybam, pybed, common
```

We must reimplement **only the surface pypgx actually touches**, in `pypgx-fuc`:

- `common.parse_variant(v)` → `(chrom, pos, ref, alt)` — used heavily in tables &
  tests. **Port first; it's tiny and central.**
- `pyvcf.VcfFrame` — read/write VCF, phasing representation, per-sample genotype
  access. Needed for `Archive` `VcfFrame[*]` payloads and `predict_alleles`.
- `pycov.CovFrame` — coverage/depth tables (`CovFrame[ReadDepth|CopyNumber|DepthOfCoverage]`).
- `pybed.BedFrame` — region BED handling (`create_regions_bed`).
- `pybam` — BAM access (likely just shells to samtools; defer).

**Action:** audit exact `fuc` call sites (`grep -rn "pyvcf\.\|pycov\.\|pybed\.\|common\."`)
and port method-by-method, test-by-test. Do **not** port all of `fuc`.

---

## 7. Test parity

Upstream acceptance tests live in [`pypgx/test.py`](./pypgx/test.py) (6 tests),
replicated in [`tests/parity.rs`](./tests/parity.rs). Because the shipped v0.26.0
data violates 3 of the 6 self-consistency assertions, each Rust test asserts the
reference's *computed* values (captured in `tests/fixtures/truth.json`) — so it
reproduces the discrepancies exactly rather than hiding them.

- [x] `test_allele_table` — duplicate detection reproduces `['19-39738787-C-T',
      '19-39738787-C-T']` (GRCh37) / `['19-39248147-C-T', ...]` (GRCh38); cores are
      position-sorted; `list_genes()` == `allele_table.Gene.unique()`.
- [x] `test_diplotype_table` — 15 diplotype genes == `PhenotypeMethod=='Diplotype'`.
- [x] `test_equation_table` — 3 equation genes == `PhenotypeMethod=='Score'`.
- [x] `test_priority_table` — reproduces `a` (incl. `MT-RNR1`) vs `b` (excl.) discrepancy.
- [x] `test_definition_table` — reproduces the `ACYP2` (both) and `CYP17A1` (GRCh38) diffs.
- [x] `test_predict_alleles` — **integration milestone.** Byte-identical output on
      both CYP4F2 archives:
      ```
      data.loc['A'] == ['*1;', '*2;', ';', '*2:19-16008388-A-C:0.5;*1:default;']  # GRCh37
      ```
- [x] Extra supporting parity: `build_definition_table` (full df), `list_variants`,
      `parse_variant`, `sort_variants`, plus 9 breadth tests in `tests/breadth.rs`.

**Differential harness:**
- [x] Reference Python env at `.refenv` (Python 3.10, pypgx 0.26.0 + fuc +
      scikit-learn). `python test.py` ⇒ 3 ok / 3 data-driven failures (documented).
- [x] [`tools/diff_harness.sh`](./tools/diff_harness.sh): regenerate ground truth
      from Python, run the Rust suite against it. **24 tests green.** (Wire into CI.)

---

## 8. Dependency mapping (Python → Rust crate)

| Python | Used for | Rust choice |
|---|---|---|
| `pandas` | DataFrames everywhere | `polars` (closest), or typed structs + `csv` for small tables. **Decide per-module**; tables are tiny → typed structs likely cleaner & more exact for NaN handling. |
| `numpy` | arrays, stats | `ndarray` / std where small. |
| `scikit-learn` | `Model[CNV]` one-vs-rest classifier (`train/test/predict_cnv`) | `linfa` or `smartcore`. **Pickled sklearn models won't load** — must retrain or re-serialize (§9). |
| `matplotlib` | `plot.py` | `plotters` (deferred, feature-gated). |
| `fuc` | VCF/cov/bed/variant primitives | `pypgx-fuc` (§6). |
| `zipfile` | Archive container | `zip` crate. |
| `subprocess` | externals | `std::process::Command`. |
| equation eval (`eval` on `Equation` col) | `predict_score` | small safe arithmetic expr evaluator (e.g. `evalexpr`, or hand-rolled). |

---

## 9. Phased roadmap

### Phase 0 — Scaffolding & reference env ✅
- [x] Create cargo crate (§3); clippy clean. *(CI wiring still TODO.)*
- [x] Copy `data/*.csv` into `data/`, embedded via `include_str!`.
- [x] Reference Python env at `.refenv` (3.10, pypgx 0.26.0 + fuc + scikit-learn);
      characterized `python test.py` (3 ok / 3 data-driven failures).
- [x] Differential-test harness (`tools/diff_harness.sh`).

### Phase 1 — `core` (pure logic) ✅ unlocks 5 of 6 tests
- [x] `parse_variant`, `sort_variants` in `fuc`.
- [x] `load_*_table` ×9 (embedded CSV; exact NaN/order/bool semantics).
- [x] `list_genes`, `list_alleles`, `list_variants`, `list_functions`, `list_phenotypes`.
- [x] `has_phenotype/has_score/has_sv`, `is_legit_allele`, `is_target_gene`.
- [x] `get_*` accessors (region, strand, paralog, exon starts/ends, function, score,
      priority, recommendation, ref/default allele, variant impact/synonyms).
      `get_recommendation` covers single- and two-gene drugs (verified on
      codeine/tacrolimus/fluvastatin).
- [x] `build_definition_table`, `collapse_alleles`, `sort_alleles`,
      `predict_phenotype`, `predict_score` (chained-comparison equation evaluator).
- [x] **Gate:** the 5 table tests reproduce the reference exactly.

### Phase 2 — `sdk` Archive + VCF primitives ✅ (read path)
- [x] `Archive { metadata, data }`, `semantic_type()`, `copy_metadata`, `check_type`.
- [x] Zip read: `from_file` — `<parent>/metadata.txt` + `data.{vcf,tsv}`.
- [x] `VcfFrame[Consolidated]` + `SampleTable[Alleles]` payloads.
- [x] `PgxError` enum (`AlleleNotFound`, `GeneNotFound`, `IncorrectSemanticType`,
      `NotTargetGene`, …).
- [x] `Archive.to_file` (zip **write**) for VcfFrame + SampleTable, with
      **Rust↔Rust** round-trip tests (`tests/roundtrip.rs`).
- [ ] **Cross-implementation check NOT done:** no test/harness has Python PyPGx
      read a Rust-written archive (or vice-versa). The round-trip is in-process
      Rust only; true cross-impl compatibility is unproven.
- [ ] CovFrame/Model payloads + remaining semantic types (needed once
      depth/CNV/pipeline work lands).

### Phase 3 — `predict_alleles` ➜ final test ✅
- [x] Port `predict_alleles(consolidated_variants)` from `utils.py`.
- [x] **Gate:** `test_predict_alleles` returns the exact expected vector for both
      CYP4F2 archives. **Headline milestone reached.**

### Phase 4 — `utils.py`: pure functions ✅ / external wrappers deferred
**Pure (done & verified):**
- [x] `call_phenotypes` — `SampleTable[Genotypes]` → `SampleTable[Phenotypes]`.
- [x] `combine_results` — merge genotypes/phenotypes/alleles/CNV → `SampleTable[Results]`.
- [x] `compare_genotypes` — Genotype/CNV concordance report.
- [x] `count_alleles` — name-sorted star-allele counts.
- [~] `print_data` / `print_metadata` — exist **only** as inline CLI handlers in
      `src/bin/pypgx.rs` (not library fns, not re-exported), and are **untested**.
      `print_data` handles only `SampleTable` payloads; Python also prints
      `CovFrame`/`VcfFrame` and raises on unsupported types.

**External-tool (deferred by design — interfaces in `src/external.rs`, return
`PgxError::NotPorted` with the exact dependency + command documented):**
- [ ] `estimate_phase_beagle` — needs `java` + bundled Beagle jar + 1KGP panel.
- [ ] `compute_control_statistics`, `compute_copy_number`, `compute_target_depth`,
      `import_read_depth`, `import_variants`, `prepare_depth_of_coverage` — `samtools`/`pycov`.
- [ ] `create_consolidated_vcf`, `create_input_vcf`, `filter_samples`,
      `slice_bam` — `bcftools`/`pyvcf`/`samtools`. (Note: `filter_samples` and
      `create_consolidated_vcf` are actually pure in-memory frame ops in Python
      and could be ported without any binary — currently stubbed regardless.)
- [ ] ⚠️ `create_regions_bed` — **MISSING entirely**, not even a `NotPorted`
      stub, despite being listed as deferred here. It is a **pure** function
      (builds a BED from the gene table) and should be ported, not deferred.

### Phase 5 — `genotype.py` ✅
- [x] `call_genotypes` — SimpleGenotyper + all 13 SV genotypers +
      `_call_duplication`/`_call_multiplication`/`_call_linked_allele`. Verified on
      the alleles-only (`AssumeNormal`) path; SV/CNV branches ported faithfully
      (exercisable once a `SampleTable[CNVCalls]` input is available).

### Phase 6 — `pipeline.py` (deferred — orchestrates external tools)
- [ ] `run_ngs_pipeline`, `run_chip_pipeline`, `run_long_read_pipeline` — stubs in
      `src/external.rs` (`NotPorted`); they chain Beagle/samtools/sklearn steps.
      The **pure** sub-pipeline (`predict_alleles → call_genotypes →
      call_phenotypes → combine_results`) already runs end-to-end via the CLI.

### Phase 7 — CLI ✅ (pipeline subset)
- [x] clap app with 11 subcommands: `list-genes`, `list-alleles`,
      `list-variants`, `predict-alleles`, `call-genotypes`, `call-phenotypes`,
      `combine-results`, `compare-genotypes`, `count-alleles`, `print-data`,
      `print-metadata`. The underlying **library** functions are verified on the
      CYP4F2 fixture, but **no test invokes the binary** — the clap parsing +
      dispatch glue is unverified (no `assert_cmd`/`CARGO_BIN_EXE` in `tests/`).
- [ ] Coverage framing: only **7 of these 11** are genuine upstream CLI
      subcommands (`predict-alleles`, `call-genotypes`, `call-phenotypes`,
      `combine-results`, `compare-genotypes`, `print-data`, `print-metadata`).
      The other 4 (`list-genes`/`list-alleles`/`list-variants`/`count-alleles`)
      are API helpers promoted to CLI — upstream has no such commands. Real
      upstream CLI coverage is **7/31** (~23%). Also `list-alleles`/`list-variants`
      hardcode the filter arg to `None`, dropping a parameter the Python API has.
- [ ] The external-tool subcommands (depth/phasing/CNV/plot) — pending Phase 4/6/8/9.

### Phase 8 — `plot.py` (deferred — matplotlib) 
- [ ] 5 plot fns (`plot_bam_copy_number`, `plot_bam_read_depth`, `plot_cn_af`,
      `plot_vcf_allele_fraction`, `plot_vcf_read_depth`) — stubs in `src/external.rs`;
      reimplement via `plotters` (visual parity, not byte parity).

### Phase 9 — CNV model (`Model[CNV]`) (deferred — scikit-learn)
- [ ] `train_cnv_caller`, `test_cnv_caller`, `predict_cnv` — stubs in `src/external.rs`.
- [ ] **Model interop:** pickled sklearn `Model[CNV]` archives can't be loaded in
      Rust → retrain to a Rust-native format (`linfa`/`smartcore`) or add a Python
      shim. Decide before implementing.

### Later — native external tools
- [ ] Port **Beagle** phasing to Rust (drop the JVM dependency).
- [ ] Native BAM/VCF/BED I/O to drop `samtools`/`bcftools` (consider `noodles`).

---

## 10. Risks & open questions

- **pandas semantics:** NaN propagation, `.apply` row order, `.unique()` ordering,
  dtype coercion, float formatting in CSV output. Highest source of subtle drift.
- **`sort_alleles`** uses a custom comparator (function/score/name) — port the exact
  ordering algorithm; tests depend on it.
- **Equation evaluation:** `equation-table` holds Python expressions; need a safe
  evaluator with identical numeric semantics (int vs float, rounding).
- **Pickled sklearn `Model[CNV]`** is not portable → retrain or shim (Phase 9).
- **Beagle determinism:** seeds/version must match for phasing parity; pin the jar.
- **`fuc` internal quirks:** VCF phasing string format, `parse_variant` edge cases —
  port against real fixtures, not assumptions.
- **Float/decimal in activity scores** (e.g. `*2:...:0.5`) — match string formatting.

---

## 11. Definition of done (v1)

- [x] All 6 `test.py` tests reproduced in Rust (asserting the reference's exact
      computed values, including the 3 data-driven discrepancies).
- [x] Differential harness shows byte-identical table loads & `predict_alleles`
      output vs upstream on the bundled fixtures (`tools/diff_harness.sh`).
- [x] CLI exposes the full **pure** pipeline (`predict-alleles` → `call-genotypes`
      → `call-phenotypes` → `combine-results` → `count-alleles`/`compare-genotypes`,
      plus `list-*`/`print-*`), verified end-to-end on the CYP4F2 fixture.
- [~] Beagle/samtools/bcftools/sklearn/matplotlib functions: **interfaces present**
      in `src/external.rs` returning `PgxError::NotPorted` with the exact
      dependency + command documented. Deferred by design (the bundled fixtures
      need no external tools); native implementation tracked in Phases 4/6/8/9.

### Scope note

Nearly all of PyPGx's public functions exist in the Rust crate (**1 missing:**
`create_regions_bed`; see [§12](#12-known-gaps-audit-2026-06-06)). The **pure
analytical core** (tables, allele/phenotype/score logic, `Archive` I/O,
`predict_alleles`, `call_genotypes`, `call_phenotypes`, `combine_results`,
`compare_genotypes`, `count_alleles`) is **implemented** and its CYP4F2-exercised
paths are **verified** against the Python reference. Everything that
fundamentally needs an external program, the `pypgx-bundle`, scikit-learn, or
matplotlib is **deferred by design** with a documented `NotPorted` stub (22 of
them in `src/external.rs`) — this is the agreed "full parity, defer externals"
scope, and cannot be verified in an environment without those binaries. Caveat:
the 22 `external` stubs are reachable only via `pypgx::external::*`, **not** the
crate root, so the re-export surface does not yet fully mirror `__init__.py`.

**Resolved during the port** (see §10 for the originals):
- pandas NaN/bool/`unique` order semantics replicated in `table.rs::Frame`
  (incl. this pandas version treating `None`/`NA`/`N/A`/`TRUE`/`FALSE` specially).
- `sort_alleles` priority/name comparators ported and verified.
- Equation evaluator handles the chained-comparison forms in `equation-table`.
- `fuc` `get_af` `j+1` ref-skip and quoted-VCF-field unquoting handled.
- Float formatting matches Python `str(float)` (`python_float_str`).

**Still open** (unchanged from §10): pickled sklearn `Model[CNV]` portability,
Beagle determinism/version pinning. (`Archive.to_file` is implemented; what
remains is a *cross-implementation* check that Python can read a Rust-written
archive — see Phase 2.)

---

## 12. Known gaps (audit 2026-06-06)

A function-by-function verification against the vendored Python and the Rust
source surfaced these deltas between earlier claims and the code. None break the
24-green suite; they are accuracy corrections + the real remaining work.

**Genuinely missing / not as claimed:**
- `create_regions_bed` — claimed as a deferred stub; **absent from `src/`
  entirely**. It is pure-portable (BED from the gene table) → should be ported.
- `external` stubs (22) are **not re-exported at the crate root** — only
  `pypgx::external::*`. So `lib.rs` does not fully mirror `__init__.py`.
- "Cross-implementation verified with Python" (Archive) — **no such test/harness
  exists**; round-trip is Rust↔Rust only.
- `VcfFrame::from_dict` — does not exist (replaced by `VcfFrame::new`).
- `has_sv(gene, allele)` overload — not ported; only the gene-level `has_sv(gene)`
  exists (the allele logic was inlined into `predict_score`).
- `Archive.check_metadata` (Python method) — not ported; the `IncorrectMetadata`
  and `BundleNotFound` `PgxError` variants have no constructing call site.

**Implemented but UNVERIFIED (present, code-faithful, no test exercises them):**
- `genotype.rs`: 12 of 13 SV genotypers + every CNV branch + the
  duplication/multiplication/linked-allele helpers — no `SampleTable[CNVCalls]`
  fixture exists, so only CYP4F2's `AssumeNormal` path runs.
- `core.rs`: `sort_alleles`, `collapse_alleles`, `get_ref_allele`,
  `get_default_allele`, `get_variant_impact`, `get_variant_synonyms` — only hit
  indirectly via the CYP4F2 `predict_alleles` path (which doesn't discriminate
  the comparators and has empty synonyms). `truth.json` even holds ground-truth
  keys for these that no test asserts.
- `predict_alleles` multi-variant VariantData ordering: Python iterates a `set`
  (hash order); Rust uses a deterministic ordered `Vec`. Proven byte-identical
  only for single-variant alleles (CYP4F2); multi-variant genes untested.
- The whole **CLI** (`src/bin/pypgx.rs`) — no test invokes the binary.
- `sdk::check_type` / `copy_metadata` — implemented, never exercised.

**Behavioral divergences from Python (no output impact on tested paths):**
- Several functions `panic!` where Python raises typed exceptions
  (`combine_results`, `call_genotypes` sample/gene-mismatch, `list_*`).
- Dropped `warnings.warn(...)` side-effects (e.g. `get_recommendation`
  gene2=None, multi-synonym in `predict_alleles`, no-CNV-calls in genotypers).
- `to_file` clones the VcfFrame before appending `##fileformat` (Python mutates
  in place); omits Python's "Saved …" print.

**Highest-value next steps:** add a `SampleTable[CNVCalls]` fixture (unlocks the
SV genotypers), a multi-variant-allele gene fixture (`predict_alleles` ordering),
and CLI smoke tests; port the pure `create_regions_bed`. ✅ `create_regions_bed`
done (2026-06-06) — see §13.

---

## 13. Native external-tool port (active goal, 2026-06-06)

**Goal: port the entire deferred surface to native Rust — everything except
Beagle phasing, which is being ported separately as `beagle-rs` and wired in
when ready.** No more `std::process::Command`; no Python, samtools, bcftools,
sklearn, or matplotlib at runtime.

### Toolchain (vendored under `repos/`, managed by repoverse)

`.repoverse.yaml` + `cargo install repoverse`; `rv init && rv link` reconstruct
the tree. `/repos/*` is gitignored except the tracked submodule gitlinks.

| Repo | Branch | Purpose |
|---|---|---|
| `repos/samtools-rs` | main | pure-Rust samtools (BAM depth/coverage/view/sort/index) |
| `repos/bcftools-rs` | main | pure-Rust bcftools (VCF view/norm/concat/filter/…) |
| `repos/sklears` | master | pure-Rust scikit-learn (SVC + OneVsRest + confusion_matrix → CNV caller) |
| `repos/htslib-rs` | main | shared HTSlib-compat layer (deduped via `rv link`) |
| `repos/noodles` | madhava/bioscript | low-level BAM/VCF/BCF/BED/bgzf/csi/tabix/fasta I/O |
| **`ruviz`** (crates.io) | 0.4.19 | plotting — `tiny-skia`+`cosmic-text`, pure Rust; replaces matplotlib |

### Layering policy (which crate for what)

- **Pure** (data manipulation only) → no I/O crate. e.g. `create_regions_bed` ✅,
  and largely `filter_samples` / `create_consolidated_vcf` / `import_variants`.
- **Raw VCF/BCF record I/O** (read/write/subset) → `noodles-vcf`/`-bcf` directly.
- **BAM depth/coverage** (`compute_target_depth`, `compute_control_statistics`,
  `prepare_depth_of_coverage`, `import_read_depth`) → **`samtools-rs`** — PyPGx
  parses the exact output of `samtools depth`/coverage, and samtools-rs is
  verified against samtools' own test suite, so reimplementing depth on raw
  `noodles-bam` would risk subtle drift. *(open decision — see below)*
- **bcftools-specific behavior** (norm, etc.) → `bcftools-rs`.
- **CNV model** → `sklears` (algorithms) — see model-interop note.
- **Plots** → `ruviz`, behind a `plots` feature.

samtools-rs/bcftools-rs are themselves built on noodles, so this is layered, not
either/or. Cargo deps are added incrementally + feature-gated so the default
build stays lean and the 24+ core tests stay fast.

### Status (`cargo test`: 43 green; `bam`/`plots` features build clean)

**Session result: 19 public functions + all 3 pipelines ported, 24 → 43 tests
green.** Every function reachable to byte-parity here is done and verified; the
5 plots render (visual parity); the CNV consumption path (`predict_cnv`,
`test_cnv_caller`) matches sklearn exactly; the 3 depth functions are
implemented (samtools-rs, verification-blocked by lack of BAM fixtures). The
3 still-open are: `train_cnv_caller` (inherently non-parity SVM training),
`create_input_vcf` (needs `bcftools call` in bcftools-rs), and
`estimate_phase_beagle` (→ beagle-rs). See "Completion state" below.

Built out a faithful `fuc.pyvcf` surface in `src/fuc.rs` (`subset`, `slice`,
`strip`, `add_af`, `unphase`, `diploidize`, `phased`, `update_chr_prefix`,
`duplicated`, `drop_duplicates`, `parse_region`, `gt_unphase`/`gt_diploidize`/
`gt_het`/`gt_pseudophase`) plus a `pybed.BedFrame` (`src/bed.rs`) and the
`_phase_extension` algorithm (`api::phase_extension`).

- [x] `create_regions_bed` — pure; byte-parity vs Python on 9 variant
      configs (`tests/regions_bed.rs`).
- [x] `filter_samples` — VcfFrame + SampleTable subset (`tests/filter_samples.rs`).
      *CovFrame branch pending the `CovFrame` payload.*
- [x] `import_variants` — in-memory VcfFrame path; byte-parity on **WGS→Imported**
      (region-slice + dup-drop) **and LongRead→Consolidated** (`_phase_extension`),
      `tests/import_variants.rs`. *Pending: bgzf+tabix file input (noodles).*
- [x] `_phase_extension` — full 111-line port; byte-verified via the LongRead
      test (anchor scoring, flip rule, `PE` field, list_alleles/list_variants filters).
- [x] `create_consolidated_vcf` — `VcfFrame::fetch`/`filter_vcf` + `_phase_extension`;
      byte-parity vs Python incl. an imported-only variant (`tests/consolidate.rs`).
      **→ the entire pure VCF-manipulation cluster of `utils.py` is now done.**
- [x] `CovFrame` SDK payload (`ArchiveData::Cov`, read/write `data.tsv`,
      `slice`/`subset`/`update_chr_prefix`) + `filter_samples` CovFrame branch.
- [x] `import_read_depth` — pure CovFrame slice; byte-parity vs Python
      (`tests/depth.rs`, fixture `depth_of_coverage.zip`).
- [x] `compute_copy_number` — pure depth normalization (intra-sample + Targeted
      inter-sample); byte-parity vs Python (`tests/copy_number.rs`).
- [~] `compute_target_depth`, `compute_control_statistics`, `prepare_depth_of_coverage`
      — **implemented** behind the `bam` feature (`src/api.rs`), using
      `CovFrame::from_bam` on samtools-rs's verified `native::depth` engine, plus
      a pure pandas-style `describe`. `cargo build --features bam` compiles clean.
      **Byte-parity with PyPGx (pysam) is unverified here** — no BAM fixtures /
      samtools to diff against; pysam depth may differ from `samtools depth`.
      The `bed`-mask (Targeted) path is not yet reproduced.
- [x] `predict_cnv` + `test_cnv_caller` — **done & verified** via path (A):
      `src/cnv.rs` evaluates the RBF OvR-SVM decision function natively and
      reproduces sklearn's `predict` **exactly** (1e-9, `tests/cnv.rs`);
      `_process_copy_number`'s `median_filter` matches **scipy** exactly; the
      `Model[CNV]` payload (`data.json`) round-trips; end-to-end verified at real
      CYP2A6 dimension. PyPGx's pickled models convert via
      `tools/convert_cnv_model.py` (one-time, in Python) → byte-parity predicts.
- [ ] `train_cnv_caller` — **inherently non-parity**: SVM *training* is
      libsvm-specific, so no Rust trainer (sklears or otherwise) reproduces
      PyPGx's shipped models. The consumption path (convert + `predict_cnv`) is
      the parity route. A sklears-backed trainer (different models) could go
      behind a `cnv` feature if from-scratch training is ever needed.
- [x] 5 × `plot_*` — ruviz behind the `plots` feature (`src/plot.rs`); all render
      valid PNGs (`tests/plots.rs`, `cargo test --features plots`), visually
      confirmed. Visual parity, not byte; the exon track + `fitted` overlay are
      omitted (ruviz lacks gridspec ratios / rect primitives).
- [x] `run_long_read_pipeline` — full native chain; byte-parity vs Python
      end-to-end (`tests/longread_pipeline.rs`, results table incl. multi-variant *4).
- [x] `run_chip_pipeline` — native for already-phased input (`src/pipeline.rs`);
      the unphased branch calls `estimate_phase_beagle` (→ beagle-rs gap).
- [x] `run_ngs_pipeline` — full orchestration (`src/pipeline.rs`); native path
      (pre-phased variants, no depth) byte-verified vs Python end-to-end
      (`tests/ngs_pipeline.rs`); the Beagle/CNV arms surface `NotPorted`
      (also asserted). Includes the MT-RNR1 pseudophase branch.
- [ ] `create_input_vcf` — **blocked**: needs `bcftools call`/`mpileup` (not in bcftools-rs).
- [ ] `estimate_phase_beagle` — **deferred to `beagle-rs`** (external, by design).

### Known blockers / open decisions

1. **Variant calling**: `create_input_vcf` needs `bcftools call`/`mpileup`, which
   bcftools-rs lists as *"not started"*. Blocks calling variants from BAM;
   already-called-VCF / chip inputs are unaffected.
2. **Pickled CNV model interop**: sklears can't load PyPGx's pickled sklearn
   `Model[CNV]` (`data.sav`). Options: retrain in sklears, or extract fitted SVM
   params (support vectors / dual coefs / intercept) once in Python and
   reimplement the decision function in Rust (byte-parity, no training).
3. **Depth backend**: resolved → **samtools-rs `native::depth`** (verified vs C
   samtools), behind the `bam` feature. Remaining gap is *verification*, not
   implementation: no BAM fixtures / pysam here to confirm pypgx byte-parity.
4. **ruviz layout gaps**: no gridspec height-ratios / rectangle primitives → the
   thin exon-annotation track + 1:10 panel ratios need workarounds (bar/area
   approximation). Parity is visual, not pixel-exact.

### Completion state (2026-06-06)

Of the ~22 non-Beagle deferred functions:

- **Implemented + byte/visually verified (16 + 3 pipelines):** `create_regions_bed`,
  `filter_samples`, `import_variants`, `_phase_extension`, `create_consolidated_vcf`,
  `import_read_depth`, `compute_copy_number`, **`predict_cnv`**, **`test_cnv_caller`**,
  5 × `plot_*`, `run_long_read_pipeline`, `run_chip_pipeline`, `run_ngs_pipeline`.
  Suite: 24 → **43 green**.
- **Implemented, verification-blocked (3, feature `bam`):** `compute_target_depth`,
  `compute_control_statistics`, `prepare_depth_of_coverage` — real code on
  samtools-rs depth; needs BAM fixtures + a samtools/pysam oracle to confirm parity.
- **Inherently non-parity (1):** `train_cnv_caller` — SVM training is libsvm-specific;
  no Rust trainer reproduces PyPGx's shipped models. Use convert + `predict_cnv` instead.
- **Integrated, feature `beagle` (1):** `estimate_phase_beagle` — shells out to the
  beagle-rs binary; wiring verified end-to-end (`tests/beagle.rs`). Unblocks
  `run_chip`/`run_ngs` on unphased input. Byte-parity vs PyPGx pending version +
  panel reconciliation (see §14).
- **Upstream-dependency blocked (1):** `create_input_vcf` — needs `bcftools call`
  in bcftools-rs (blocker #1).

## 14. Beagle integration (`repos/beagle-rs`) — wired in ✅ (behind `beagle`)

**Done (2026-06-07):** `beagle-rs`'s `gt=` phasing CLI is complete and a verified
drop-in for the Beagle jar (its port reached "structurally complete + byte-parity
vs Java" on the official fixtures). `estimate_phase_beagle` is now implemented in
`src/external.rs` behind the **`beagle` feature**: it writes the imported
VcfFrame to a temp `input.vcf`, invokes the `beagle-rs` binary
(`gt=/chrom=/[ref=]/out=/impute=/em=`, EM-skip retry), reads the bgzf
`output.vcf.gz` (via `flate2::MultiGzDecoder`), and returns the `VcfFrame[Phased]`
archive. Integration verified end-to-end (`tests/beagle.rs` — phases a CYP4F2
Imported frame; CI `beagle` job builds beagle-rs + runs it). With the feature on,
`run_chip_pipeline` / `run_ngs_pipeline` run end-to-end on **unphased** input.

**Still NOT byte-parity with PyPGx (deliberately flagged, not silently shipped):**
1. **Beagle version mismatch** — PyPGx bundles `22Jul22.46e`; beagle-rs targets
   `27Feb25.75f`. Phasing output can differ across versions; reconcile before
   claiming parity vs PyPGx's reference.
2. **Reference panel** — PyPGx always passes `ref=<1KGP panel>` from `pypgx-bundle`
   (absent here). The current impl supports `panel=Some(path)` (reference-based)
   and `panel=None` (pure phasing); PyPGx's `None` means "load from bundle".
3. **`chr` prefix + Chip GSA filtering** — PyPGx's panel-driven chr-prefix
   add/remove and `Platform=='Chip'` GSA allele filtering are not yet reproduced.

Binary discovery: `$BEAGLE_RS_BIN`, else `beagle-rs` on PATH.

---

**Original integration contract** (from `utils.estimate_phase_beagle`): input is a
`VcfFrame[Imported]`; PyPGx writes a temp `input.vcf` and runs
`java -Xmx2g -jar beagle.jar gt=input.vcf chrom=<region> ref=<panel>
out=output impute=<bool> em=<bool>`, then reads `output.vcf.gz` →
`VcfFrame[Phased]` (metadata `SemanticType=VcfFrame[Phased]`, `Program=Beagle`),
with: an EM-retry fallback (`em=true` then `em=false` on error), chr-prefix
add/remove around the panel, a single-marker guard, and `Platform=='Chip'` GSA
allele filtering. The beagle-rs binary is a **drop-in** for the jar (same
`gt=/chrom=/ref=/out=/impute=/em=` args → same `output.vcf.gz`).

**When beagle-rs's DoD is met:** replace the `estimate_phase_beagle` `NotPorted`
stub in `src/external.rs` with a real impl behind a `beagle` feature: write the
imported VcfFrame to a temp `input.vcf`, invoke the `beagle-rs` binary
(`repos/beagle-rs` `beagle-rs-cli`) with the same args + EM-retry, read back
`output.vcf.gz` via `VcfFrame::from_string`, and return the `[Phased]` archive.
The `run_chip`/`run_ngs` pipelines then run end-to-end on unphased input.

**Caveats to reconcile before claiming byte-parity:**
1. **Beagle version mismatch** — PyPGx bundles `beagle.22Jul22.46e.jar`; beagle-rs
   targets `27Feb25.75f`. Phasing output can differ across Beagle versions, so a
   straight swap may not byte-match PyPGx's reference. Either point beagle-rs at
   the 22Jul22 behavior or accept/verify against 27Feb25.
2. **Reference panel** — needs the 1KGP panel from `pypgx-bundle`
   (`1kgp/<assembly>/<gene>.vcf.gz`), not present here.
3. Read `output.vcf.gz` requires **bgzf** decode (noodles-bgzf) on the Rust side,
   or have beagle-rs emit plain VCF for the wrapper.
