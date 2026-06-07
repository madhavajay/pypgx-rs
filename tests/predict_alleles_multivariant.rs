//! Multi-variant `predict_alleles` parity. The core suite only covers CYP4F2
//! (single-variant alleles), so the `VariantData` for an allele defined by
//! *multiple* variants was unverified.
//!
//! Investigated against PyPGx 0.26.0: each allele's defining variants are stored
//! in a Python **set** (`utils.py:1131`), so the order they appear in
//! `VariantData` is `PYTHONHASHSEED`-dependent — it is NOT stable across PyPGx
//! runs (seeds 0/1/2 → `45411941,45412079`; seed 3 → the reverse). Exact
//! byte-parity on that order is therefore neither achievable nor meaningful.
//!
//! What IS invariant: the *set* of variants and their fractions. Rust emits them
//! in a deterministic position-sorted order (a strict improvement over upstream's
//! nondeterminism), which equals PyPGx under the majority of hash seeds. This
//! test asserts the haplotypes match exactly and the VariantData variant set +
//! fractions match, order-insensitively.

use std::collections::HashSet;

use pypgx::fuc::VcfFrame;
use pypgx::sdk::{Archive, ArchiveData};

#[test]
fn apoe_multivariant_predict_alleles_matches_pypgx() {
    let cols: Vec<String> = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let row = |pos: &str, r: &str, a: &str| -> Vec<String> {
        vec!["19", pos, ".", r, a, ".", ".", ".", "GT", "1|1"]
            .into_iter()
            .map(|s| s.to_string())
            .collect()
    };
    let rows = vec![row("45411941", "T", "C"), row("45412079", "C", "T")];

    let consolidated = Archive::new(
        vec![
            ("Platform".into(), "WGS".into()),
            ("Gene".into(), "APOE".into()),
            ("Assembly".into(), "GRCh37".into()),
            ("SemanticType".into(), "VcfFrame[Consolidated]".into()),
        ],
        ArchiveData::Vcf(VcfFrame::new(Vec::new(), cols, rows)),
    );

    let result = pypgx::predict_alleles(&consolidated).expect("predict_alleles");
    let t = result.as_sample_table();
    let r = t.loc("S1");
    let col = |name: &str| t.columns.iter().position(|c| c == name).unwrap();

    // Haplotypes are order-stable in both impls.
    assert_eq!(r[col("Haplotype1")], "E3;");
    assert_eq!(r[col("Haplotype2")], "E3;");
    assert_eq!(r[col("AlternativePhase")], ";");

    // VariantData: "E3:<v1,v2,...>:<f1,f2,...>;" — compare order-insensitively.
    let vd = r[col("VariantData")].trim_end_matches(';');
    let parts: Vec<&str> = vd.split(':').collect();
    assert_eq!(parts[0], "E3", "allele label");
    let variants: HashSet<&str> = parts[1].split(',').collect();
    let expected: HashSet<&str> = ["19-45411941-T-C", "19-45412079-C-T"].into_iter().collect();
    assert_eq!(variants, expected, "VariantData variant set must match PyPGx");
    assert!(parts[2].split(',').all(|f| f == "nan"), "fractions are nan (no AF): {:?}", parts[2]);

    // Rust's own order is deterministic: position-sorted (45411941 < 45412079).
    assert_eq!(parts[1], "19-45411941-T-C,19-45412079-C-T", "Rust order is position-sorted");
}
