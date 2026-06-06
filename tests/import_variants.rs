//! Byte-parity test for `import_variants` (in-memory VcfFrame path) vs the
//! Python reference. The input VCF exercises region slicing (one variant out of
//! the CYP4F2 region) and duplicate dropping (a repeated variant). Ground truth
//! generated from pypgx 0.26.0 in `.refenv`.

use pypgx::fuc::VcfFrame;
use serde_json::Value;

const INPUT_VCF: &str = include_str!("fixtures/import_variants_input.vcf");
const TRUTH: &str = include_str!("fixtures/import_variants.json");
const PE_VCF: &str = include_str!("fixtures/phase_ext_input.vcf");
const PE_TRUTH: &str = include_str!("fixtures/phase_ext.json");

fn rows_of(v: &Value, key: &str) -> Vec<Vec<String>> {
    v[key]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|c| c.as_str().unwrap().to_string())
                .collect()
        })
        .collect()
}

fn cols_of(v: &Value, key: &str) -> Vec<String> {
    v[key]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn import_variants_wgs_matches_python() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();
    let vf = VcfFrame::from_string(INPUT_VCF);
    let archive = pypgx::import_variants("CYP4F2", &vf, "GRCh37", "WGS", None, false).unwrap();

    assert_eq!(archive.semantic_type(), t["semantic_type"].as_str().unwrap());

    let out = archive.as_vcf();
    let expected_cols: Vec<String> = t["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
    assert_eq!(out.columns, expected_cols, "columns");

    let expected_rows: Vec<Vec<String>> = t["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|c| c.as_str().unwrap().to_string())
                .collect()
        })
        .collect();
    assert_eq!(out.rows, expected_rows, "rows");

    // Metadata key order + values (Platform, Gene, Assembly, SemanticType).
    assert_eq!(
        archive.metadata,
        vec![
            ("Platform".to_string(), "WGS".to_string()),
            ("Gene".to_string(), "CYP4F2".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "VcfFrame[Imported]".to_string()),
        ]
    );
}

/// LongRead path → `_phase_extension`: a phased anchor (19-16008388) lets the
/// unphased het (19-15990431, a *4-defining variant) be oriented `1|0`, and the
/// hom-ref call is pseudo-phased. Exercises anchor scoring, the flip rule, the
/// per-row `PE` FORMAT, and `list_alleles`/`list_variants` variant/allele filters.
#[test]
fn import_variants_longread_phase_extension_matches_python() {
    let t: Value = serde_json::from_str(PE_TRUTH).unwrap();
    let vf = VcfFrame::from_string(PE_VCF);
    let archive = pypgx::import_variants("CYP4F2", &vf, "GRCh37", "LongRead", None, false).unwrap();

    assert_eq!(archive.semantic_type(), "VcfFrame[Consolidated]");
    let out = archive.as_vcf();
    assert_eq!(out.columns, cols_of(&t, "columns"), "columns");
    assert_eq!(out.rows, rows_of(&t, "rows"), "rows");
}
