//! End-to-end parity for the pure analytical pipeline:
//! `predict_alleles` → `call_genotypes` → `call_phenotypes` →
//! `combine_results` → `count_alleles`, verified against the Python reference
//! on the CYP4F2 fixture.

fn fixture(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn vertical_slice_matches_python() {
    let consolidated = pypgx::Archive::from_file(&fixture("CYP4F2-GRCh37.zip")).unwrap();
    let alleles = pypgx::predict_alleles(&consolidated).unwrap();
    let genotypes = pypgx::call_genotypes(Some(&alleles), None).unwrap();
    let phenotypes = pypgx::call_phenotypes(&genotypes).unwrap();
    let results =
        pypgx::combine_results(Some(&genotypes), Some(&phenotypes), Some(&alleles), None).unwrap();

    // Genotype + phenotype.
    assert_eq!(
        genotypes.as_sample_table().loc("A"),
        &vec!["*1/*2".to_string()]
    );
    assert_eq!(
        phenotypes.as_sample_table().loc("A"),
        &vec!["Indeterminate".to_string()]
    );

    // Combined results row (CNV is empty / NaN — no CNV calls provided).
    let rt = results.as_sample_table();
    assert_eq!(
        rt.columns,
        vec![
            "Genotype",
            "Phenotype",
            "Haplotype1",
            "Haplotype2",
            "AlternativePhase",
            "VariantData",
            "CNV"
        ]
    );
    assert_eq!(
        rt.loc("A"),
        &vec![
            "*1/*2".to_string(),
            "Indeterminate".to_string(),
            "*1;".to_string(),
            "*2;".to_string(),
            ";".to_string(),
            "*2:19-16008388-A-C:0.5;*1:default;".to_string(),
            String::new(),
        ]
    );

    // Allele counts, name-sorted.
    let counts = pypgx::count_alleles(&results);
    assert_eq!(counts, vec![("*1".to_string(), 1), ("*2".to_string(), 1)]);
}

#[test]
fn compare_genotypes_self_concordant() {
    // A results archive compared with itself is 100% concordant on Genotype;
    // CNV is empty so it has 0 comparable rows (Concordance: N/A).
    let consolidated = pypgx::Archive::from_file(&fixture("CYP4F2-GRCh37.zip")).unwrap();
    let alleles = pypgx::predict_alleles(&consolidated).unwrap();
    let genotypes = pypgx::call_genotypes(Some(&alleles), None).unwrap();
    let phenotypes = pypgx::call_phenotypes(&genotypes).unwrap();
    let results =
        pypgx::combine_results(Some(&genotypes), Some(&phenotypes), Some(&alleles), None).unwrap();

    let report = pypgx::compare_genotypes(&results, &results, false);
    assert!(report.contains("# Genotype"));
    assert!(report.contains("Concordance: 1.000 (1/1)"));
    assert!(report.contains("# CNV"));
    assert!(report.contains("Concordance: N/A"));
}
