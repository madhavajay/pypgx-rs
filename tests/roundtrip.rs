//! Round-trip tests for `Archive` write/read (`to_file` / `from_file`),
//! exercising both the VcfFrame and SampleTable payload paths.

use pypgx::sdk::ArchiveData;

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn tmp(name: &str) -> String {
    format!("{}/{name}", std::env::temp_dir().display())
}

#[test]
fn vcf_archive_roundtrips() {
    // Read a VcfFrame[Consolidated] fixture, write it back, read again.
    let a = pypgx::Archive::from_file(&fixture("CYP4F2-GRCh37.zip")).unwrap();
    let out = tmp("pypgx_rt_vcf.zip");
    a.to_file(&out).unwrap();
    let b = pypgx::Archive::from_file(&out).unwrap();

    assert_eq!(a.metadata, b.metadata);
    let (va, vb) = (a.as_vcf(), b.as_vcf());
    assert_eq!(va.columns, vb.columns);
    assert_eq!(va.rows, vb.rows);
    assert_eq!(va.samples(), vb.samples());
    // Allele fraction survives the round-trip.
    assert_eq!(vb.get_af("A", "19-16008388-A-C"), Some(0.5));
    std::fs::remove_file(&out).ok();
}

#[test]
fn sample_table_archive_roundtrips() {
    // predict_alleles output is a SampleTable[Alleles]; write + read it back.
    let input = pypgx::Archive::from_file(&fixture("CYP4F2-GRCh37.zip")).unwrap();
    let result = pypgx::predict_alleles(&input).unwrap();
    assert!(matches!(result.data, ArchiveData::SampleTable(_)));

    let out = tmp("pypgx_rt_alleles.zip");
    result.to_file(&out).unwrap();
    let back = pypgx::Archive::from_file(&out).unwrap();

    assert_eq!(back.semantic_type(), "SampleTable[Alleles]");
    let t = back.as_sample_table();
    assert_eq!(
        t.loc("A"),
        &vec![
            "*1;".to_string(),
            "*2;".to_string(),
            ";".to_string(),
            "*2:19-16008388-A-C:0.5;*1:default;".to_string(),
        ]
    );
    std::fs::remove_file(&out).ok();
}
