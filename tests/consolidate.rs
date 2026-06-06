//! Byte-parity test for `create_consolidated_vcf` vs the Python reference.
//! Merges an Imported archive (built via `import_variants`) with a synthesized
//! VcfFrame[Phased] (standing in for Beagle output). One imported-only variant
//! exercises the `filter_vcf` merge + `_phase_extension` path.

use pypgx::fuc::VcfFrame;
use pypgx::sdk::{Archive, ArchiveData};
use serde_json::Value;

const IMPORTED_VCF: &str = include_str!("fixtures/consolidate_imported.vcf");
const PHASED_VCF: &str = include_str!("fixtures/consolidate_phased.vcf");
const TRUTH: &str = include_str!("fixtures/consolidate.json");

fn meta(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
    pairs
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

#[test]
fn create_consolidated_vcf_matches_python() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();

    let imp_vf = VcfFrame::from_string(IMPORTED_VCF);
    let imported = pypgx::import_variants("CYP4F2", &imp_vf, "GRCh37", "WGS", None, false).unwrap();

    let phased = Archive::new(
        meta(&[
            ("Platform", "WGS"),
            ("Gene", "CYP4F2"),
            ("Assembly", "GRCh37"),
            ("SemanticType", "VcfFrame[Phased]"),
        ]),
        ArchiveData::Vcf(VcfFrame::from_string(PHASED_VCF)),
    );

    let res = pypgx::create_consolidated_vcf(&imported, &phased).unwrap();
    assert_eq!(res.semantic_type(), "VcfFrame[Consolidated]");

    let out = res.as_vcf();
    let expected_cols: Vec<String> = t["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
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
    assert_eq!(out.columns, expected_cols, "columns");
    assert_eq!(out.rows, expected_rows, "rows");

    // Metadata carried from the phased archive, SemanticType updated.
    assert_eq!(
        res.metadata,
        meta(&[
            ("Platform", "WGS"),
            ("Gene", "CYP4F2"),
            ("Assembly", "GRCh37"),
            ("SemanticType", "VcfFrame[Consolidated]"),
        ])
    );
}
