//! Byte-parity test for `compute_copy_number` (pure depth normalization) vs
//! Python. Reads a `CovFrame[ReadDepth]` + a `SampleTable[Statistics]` fixture
//! (chosen for clean divisions) and checks the `CovFrame[CopyNumber]` output.

use serde_json::Value;

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

const TRUTH: &str = include_str!("fixtures/copy_number.json");

#[test]
fn compute_copy_number_matches_python() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();
    let rd = pypgx::Archive::from_file(&fixture("cn_read_depth.zip")).unwrap();
    let stats = pypgx::Archive::from_file(&fixture("cn_stats.zip")).unwrap();

    let res = pypgx::compute_copy_number(&rd, &stats, None).unwrap();
    assert_eq!(res.semantic_type(), "CovFrame[CopyNumber]");
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
    assert_eq!(cf.rows, rows, "copy-number rows");

    assert_eq!(
        res.metadata,
        vec![
            ("Gene".to_string(), "CYP4F2".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "CovFrame[CopyNumber]".to_string()),
            ("Platform".to_string(), "WGS".to_string()),
            ("Control".to_string(), "EGFR".to_string()),
            ("Samples".to_string(), "None".to_string()),
        ]
    );
}
