//! Integration test for `estimate_phase_beagle` (feature `beagle`): shells out
//! to the beagle-rs binary to phase an Imported VcfFrame and reads back the
//! VcfFrame[Phased]. Requires the beagle-rs binary — built from the submodule
//! and pointed at via $BEAGLE_RS_BIN (set here to the release build).
//!
//! This verifies the *integration wiring* (invoke → phased VCF → archive), NOT
//! byte-parity with PyPGx: PyPGx bundles Beagle 22Jul22.46e while beagle-rs
//! targets 27Feb25.75f, and no 1KGP panel is present (pure gt= phasing).
#![cfg(feature = "beagle")]

use pypgx::fuc::VcfFrame;
use pypgx::sdk::{Archive, ArchiveData};

/// A small Imported VcfFrame: 4 samples, 5 sorted CYP4F2-region markers, het calls.
fn imported() -> Archive {
    let cols: Vec<String> = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2", "S3",
        "S4",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    let mk = |pos: &str, r: &str, a: &str, g: [&str; 4]| {
        let mut row = vec![
            "19".into(), pos.into(), ".".into(), r.into(), a.into(), ".".into(), ".".into(),
            ".".into(), "GT:AD:DP:AF".into(),
        ];
        for gg in g {
            row.push(format!("{gg}:10,10:20:0.5"));
        }
        row
    };
    // Sorted by POS (Beagle requires sorted markers).
    let rows = vec![
        mk("15990431", "C", "T", ["0/1", "0/0", "1/1", "0/1"]),
        mk("15995000", "T", "A", ["0/0", "0/1", "1/1", "0/1"]),
        mk("16000000", "G", "A", ["1/1", "0/1", "0/0", "0/1"]),
        mk("16008388", "A", "C", ["0/1", "1/1", "0/0", "0/1"]),
        mk("16010000", "G", "C", ["0/1", "0/1", "0/1", "1/1"]),
    ];
    Archive::new(
        vec![
            ("Platform".into(), "WGS".into()),
            ("Gene".into(), "CYP4F2".into()),
            ("Assembly".into(), "GRCh37".into()),
            ("SemanticType".into(), "VcfFrame[Imported]".into()),
        ],
        ArchiveData::Vcf(VcfFrame::new(Vec::new(), cols, rows)),
    )
}

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

/// Regression test for the GSTM1/HTR1A crash class: when 0 or 1 input markers
/// overlap the reference panel, the real Beagle (Java *and* beagle-rs) aborts
/// with a window/single-marker assert. PyPGx works around this in its wrapper by
/// short-circuiting those cases; `estimate_phase_beagle` mirrors that. These two
/// cases return *before* the beagle-rs binary is ever invoked, so the guarantee
/// holds for every gene regardless of which Beagle build is installed — and the
/// test needs no binary.
#[test]
fn overlap_precheck_never_crashes() {
    // 0 overlap: panel_zero has only 19-99999999-G-T, which none of the 5
    // imported markers match → statistical phasing is skipped, empty Phased.
    let phased = pypgx::external::estimate_phase_beagle(&imported(), Some(&fixture("panel_zero.vcf.gz")), false)
        .expect("0-overlap must not crash");
    assert_eq!(phased.semantic_type(), "VcfFrame[Phased]");
    assert!(
        phased.as_vcf().rows.is_empty(),
        "0-overlap should yield an empty Phased frame, got {} rows",
        phased.as_vcf().rows.len()
    );

    // 1 overlap: panel_one overlaps only at 19-16008388-A-C → that single row is
    // pseudo-phased (no statistical phasing), the rest dropped.
    let phased = pypgx::external::estimate_phase_beagle(&imported(), Some(&fixture("panel_one.vcf.gz")), false)
        .expect("1-overlap must not crash");
    assert_eq!(phased.semantic_type(), "VcfFrame[Phased]");
    let vf = phased.as_vcf();
    assert_eq!(vf.rows.len(), 1, "1-overlap should keep exactly the overlapping marker");
    assert_eq!(vf.rows[0][1], "16008388", "kept marker should be the overlapping POS");
    // Pseudo-phased + stripped to GT only.
    let gt = vf.rows[0][9].split(':').next().unwrap_or("");
    assert!(gt.contains('|'), "kept marker must be pseudo-phased, got {:?}", vf.rows[0][9]);
}

#[test]
fn estimate_phase_beagle_phases_via_binary() {
    let bin = format!(
        "{}/repos/beagle-rs/target/release/beagle-rs",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&bin).exists() {
        eprintln!("skip: beagle-rs binary not built at {bin}");
        return;
    }
    std::env::set_var("BEAGLE_RS_BIN", &bin);

    let phased = pypgx::external::estimate_phase_beagle(&imported(), None, false)
        .expect("beagle phasing");

    assert_eq!(phased.semantic_type(), "VcfFrame[Phased]");
    assert_eq!(phased.get("Program"), Some("Beagle"));
    let vf = phased.as_vcf();
    assert!(!vf.rows.is_empty(), "phased output should be non-empty");
    // Every sample call must now be phased (GT uses '|').
    for r in &vf.rows {
        for cell in &r[9..] {
            let gt = cell.split(':').next().unwrap_or("");
            assert!(gt.contains('|'), "expected phased GT, got {cell:?}");
        }
    }
}
