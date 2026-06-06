//! Byte-parity test for `import_read_depth` (pure CovFrame slicing) vs Python.
//! Reads a `CovFrame[DepthOfCoverage]` fixture archive (positions in and out of
//! the CYP4F2 region) and checks the sliced `CovFrame[ReadDepth]` output.

use serde_json::Value;

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

const TRUTH: &str = include_str!("fixtures/import_read_depth.json");

#[test]
fn import_read_depth_matches_python() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();
    let doc = pypgx::Archive::from_file(&fixture("depth_of_coverage.zip")).unwrap();
    assert_eq!(doc.semantic_type(), "CovFrame[DepthOfCoverage]");

    let res = pypgx::import_read_depth("CYP4F2", &doc, None, false).unwrap();
    let cf = res.as_cov();

    let cols: Vec<String> = t["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
    let rows: Vec<Vec<String>> = t["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| {
            r.as_array()
                .unwrap()
                .iter()
                .map(|c| c.as_str().unwrap().to_string())
                .collect()
        })
        .collect();
    assert_eq!(cf.columns, cols, "columns");
    assert_eq!(cf.rows, rows, "rows (region slice)");

    assert_eq!(
        res.metadata,
        vec![
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "CovFrame[ReadDepth]".to_string()),
            ("Gene".to_string(), "CYP4F2".to_string()),
        ]
    );
}
