//! Smoke tests for the ruviz-backed plots (feature `plots`). Output is visual,
//! not byte-comparable, so we assert each function renders a non-empty PNG per
//! sample. Run with: `cargo test --features plots --test plots`.
#![cfg(feature = "plots")]

use pypgx::fuc::VcfFrame;

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn tmpdir(tag: &str) -> String {
    let d = format!("{}/pypgx_plots_{tag}", std::env::temp_dir().display());
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A PNG file exists, is non-empty, and starts with the PNG magic bytes.
fn assert_png(path: &str) {
    let bytes = std::fs::read(path).unwrap_or_else(|_| panic!("missing {path}"));
    assert!(bytes.len() > 100, "{path} too small ({} bytes)", bytes.len());
    assert_eq!(&bytes[..4], b"\x89PNG", "{path} is not a PNG");
}

const IV_VCF: &str = include_str!("fixtures/import_variants_input.vcf");

#[test]
fn all_plots_render_pngs() {
    let dir = tmpdir("render");

    // CovFrame[ReadDepth] → read-depth line.
    let rd = pypgx::Archive::from_file(&fixture("cn_read_depth.zip")).unwrap();
    let paths = pypgx::plot_bam_read_depth(&rd, Some(&dir), None).unwrap();
    assert_eq!(paths.len(), 2); // samples A, B
    paths.iter().for_each(|p| assert_png(p));

    // CovFrame[CopyNumber] (computed) → copy-number line.
    let stats = pypgx::Archive::from_file(&fixture("cn_stats.zip")).unwrap();
    let cn = pypgx::compute_copy_number(&rd, &stats, None).unwrap();
    pypgx::plot_bam_copy_number(&cn, Some(&dir), None)
        .unwrap()
        .iter()
        .for_each(|p| assert_png(p));

    // VcfFrame DP → read-depth scatter.
    let vf = VcfFrame::from_string(IV_VCF);
    pypgx::plot_vcf_read_depth("CYP4F2", &vf, "GRCh37", Some(&dir), None)
        .unwrap()
        .iter()
        .for_each(|p| assert_png(p));

    // Imported VcfFrame → allele-fraction scatter (REF/ALT).
    let imported = pypgx::import_variants("CYP4F2", &vf, "GRCh37", "WGS", None, false).unwrap();
    pypgx::plot_vcf_allele_fraction(&imported, Some(&dir), None)
        .unwrap()
        .iter()
        .for_each(|p| assert_png(p));

    // plot_cn_af: copy number + allele fraction, 2-panel. Samples differ
    // between cn (A,B) and imported (A); request sample A explicitly.
    let a = vec!["A".to_string()];
    pypgx::plot_cn_af(&cn, &imported, Some(&dir), Some(&a))
        .unwrap()
        .iter()
        .for_each(|p| assert_png(p));
}
