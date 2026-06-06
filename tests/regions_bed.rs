//! Byte-parity test for `create_regions_bed` against the Python reference.
//! Ground truth (`tests/fixtures/regions_bed.json`) was generated from
//! pypgx 0.26.0 in `.refenv` (see tools/). Each entry is the full
//! `BedFrame.gr.df` rendered as string rows.

use pypgx::create_regions_bed;
use serde_json::Value;

const TRUTH: &str = include_str!("fixtures/regions_bed.json");

fn truth() -> Value {
    serde_json::from_str(TRUTH).expect("parse regions_bed.json")
}

/// Compare a BedFrame's string rows to a fixture key's expected rows.
fn assert_matches(key: &str, bf: &pypgx::BedFrame) {
    let t = truth();
    let expected_rows: Vec<Vec<String>> = t[key]
        .as_array()
        .unwrap_or_else(|| panic!("missing fixture key {key}"))
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|c| c.as_str().unwrap().to_string())
                .collect()
        })
        .collect();
    let expected_cols: Vec<String> = t[format!("{key}__cols")]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
    assert_eq!(bf.columns, expected_cols, "columns ({key})");
    assert_eq!(bf.rows, expected_rows, "rows ({key})");
}

fn g(genes: &[&str]) -> Vec<String> {
    genes.iter().map(|s| s.to_string()).collect()
}

#[test]
fn regions_bed_grch37() {
    let bf = create_regions_bed("GRCh37", false, false, false, false, false, None, false);
    assert_matches("grch37", &bf);
}

#[test]
fn regions_bed_grch38() {
    let bf = create_regions_bed("GRCh38", false, false, false, false, false, None, false);
    assert_matches("grch38", &bf);
}

#[test]
fn regions_bed_chr_prefix() {
    let bf = create_regions_bed("GRCh37", true, false, false, false, false, None, false);
    assert_matches("chr", &bf);
}

#[test]
fn regions_bed_merge() {
    let bf = create_regions_bed("GRCh37", false, true, false, false, false, None, false);
    assert_matches("merge", &bf);
}

#[test]
fn regions_bed_target_sv_var() {
    assert_matches(
        "target",
        &create_regions_bed("GRCh37", false, false, true, false, false, None, false),
    );
    assert_matches(
        "sv",
        &create_regions_bed("GRCh37", false, false, false, true, false, None, false),
    );
    assert_matches(
        "var",
        &create_regions_bed("GRCh37", false, false, false, false, true, None, false),
    );
}

#[test]
fn regions_bed_gene_filter_and_exclude() {
    assert_matches(
        "genes_cyp2d6_cyp2c9",
        &create_regions_bed(
            "GRCh37",
            false,
            false,
            false,
            false,
            false,
            Some(&g(&["CYP2D6", "CYP2C9"])),
            false,
        ),
    );
    assert_matches(
        "exclude_cyp2d6",
        &create_regions_bed(
            "GRCh37",
            false,
            false,
            false,
            false,
            false,
            Some(&g(&["CYP2D6"])),
            true,
        ),
    );
}
