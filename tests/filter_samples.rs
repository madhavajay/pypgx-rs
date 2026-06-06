//! Tests for `filter_samples` on VcfFrame and SampleTable archives, using the
//! single-sample ("A") CYP4F2 fixture.

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn samples(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn filter_vcf_include_and_exclude() {
    let a = pypgx::Archive::from_file(&fixture("CYP4F2-GRCh37.zip")).unwrap();
    let nrows = a.as_vcf().rows.len();

    // Include "A": keep the sample; 9 headers + 1 sample column; rows intact.
    let inc = pypgx::filter_samples(&a, &samples(&["A"]), false);
    assert_eq!(inc.as_vcf().samples(), vec!["A".to_string()]);
    assert_eq!(inc.as_vcf().columns.len(), 10);
    assert_eq!(inc.as_vcf().rows.len(), nrows);
    assert!(inc.as_vcf().rows.iter().all(|r| r.len() == 10));

    // Exclude "A": no samples remain; only the 9 fixed VCF columns.
    let exc = pypgx::filter_samples(&a, &samples(&["A"]), true);
    assert!(exc.as_vcf().samples().is_empty());
    assert_eq!(exc.as_vcf().columns.len(), 9);
    assert!(exc.as_vcf().rows.iter().all(|r| r.len() == 9));
    assert_eq!(exc.as_vcf().rows.len(), nrows);

    // Metadata is carried through unchanged.
    assert_eq!(inc.metadata, a.metadata);
}

#[test]
fn filter_sample_table_include_and_exclude() {
    let a = pypgx::Archive::from_file(&fixture("CYP4F2-GRCh37.zip")).unwrap();
    let alleles = pypgx::predict_alleles(&a).unwrap(); // SampleTable[Alleles]

    let inc = pypgx::filter_samples(&alleles, &samples(&["A"]), false);
    assert_eq!(inc.as_sample_table().index, vec!["A".to_string()]);
    assert_eq!(
        inc.as_sample_table().columns,
        alleles.as_sample_table().columns
    );

    let exc = pypgx::filter_samples(&alleles, &samples(&["A"]), true);
    assert!(exc.as_sample_table().index.is_empty());
    assert!(exc.as_sample_table().rows.is_empty());

    // Excluding a non-present sample keeps everyone.
    let keep = pypgx::filter_samples(&alleles, &samples(&["ZZZ"]), true);
    assert_eq!(keep.as_sample_table().index, vec!["A".to_string()]);
}
