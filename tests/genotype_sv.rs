//! SV/CNV genotyper coverage. Until now the 13 SV genotypers + the alleles↔CNV
//! merge in `call_genotypes` were exercised by no test (no `SampleTable[CNVCalls]`
//! fixture existed). These cases assert the Rust output matches PyPGx 0.26.0's
//! `call_genotypes` truth, captured in `.refenv` for the same inputs.

use pypgx::sdk::{Archive, ArchiveData, SampleTable};

fn table(meta: &[(&str, &str)], cols: &[&str], index: &[&str], rows: Vec<Vec<&str>>) -> Archive {
    Archive::new(
        meta.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect(),
        ArchiveData::SampleTable(SampleTable {
            index: index.iter().map(|s| s.to_string()).collect(),
            columns: cols.iter().map(|s| s.to_string()).collect(),
            rows: rows.into_iter().map(|r| r.into_iter().map(|s| s.to_string()).collect()).collect(),
        }),
    )
}

/// CNV-only path `(None, Some(cnv))` + the GSTT1 genotyper (fixed diplotype per
/// CNV state). PyPGx truth: Normal→*A/*A, WholeDel1→*0/*A, WholeDel1Hom→*0/*0.
#[test]
fn gsst1_cnv_only_matches_pypgx() {
    let cnv = table(
        &[("Gene", "GSTT1"), ("Assembly", "GRCh37"), ("SemanticType", "SampleTable[CNVCalls]")],
        &["CNV"],
        &["Normal", "WholeDel1", "WholeDel1Hom"],
        vec![vec!["Normal"], vec!["WholeDel1"], vec!["WholeDel1Hom"]],
    );
    let gt = pypgx::call_genotypes(None, Some(&cnv)).expect("call_genotypes");
    let t = gt.as_sample_table();
    assert_eq!(t.loc("Normal"), &vec!["*A/*A".to_string()]);
    assert_eq!(t.loc("WholeDel1"), &vec!["*0/*A".to_string()]);
    assert_eq!(t.loc("WholeDel1Hom"), &vec!["*0/*0".to_string()]);
}

/// Alleles + CNV merge path `(Some, Some)` + the CYP4F2 genotyper (`WholeDel1`
/// pairs the priority-first allele with `*DEL`). PyPGx truth: *1/*3 + WholeDel1
/// → *3/*DEL; alleles-only → *1/*3.
#[test]
fn cyp4f2_alleles_plus_cnv_matches_pypgx() {
    let meta = [("Gene", "CYP4F2"), ("Assembly", "GRCh37"), ("SemanticType", "SampleTable[Alleles]")];
    let alleles = table(
        &meta,
        &["Haplotype1", "Haplotype2", "AlternativePhase", "VariantData"],
        &["S1"],
        vec![vec!["*3;*1;", "*1;", ";", "*1:default;*3:19-15990431-C-T:1.0;"]],
    );

    // alleles-only → *1/*3
    let gt = pypgx::call_genotypes(Some(&alleles), None).expect("call_genotypes alleles-only");
    assert_eq!(gt.as_sample_table().loc("S1"), &vec!["*1/*3".to_string()]);

    // alleles + WholeDel1 → *3/*DEL
    let cnv = table(
        &[("Gene", "CYP4F2"), ("Assembly", "GRCh37"), ("SemanticType", "SampleTable[CNVCalls]")],
        &["CNV"],
        &["S1"],
        vec![vec!["WholeDel1"]],
    );
    let gt = pypgx::call_genotypes(Some(&alleles), Some(&cnv)).expect("call_genotypes merge");
    assert_eq!(gt.as_sample_table().loc("S1"), &vec!["*3/*DEL".to_string()]);
}

/// Broad SV-genotyper coverage exercising the shared helpers (`call_duplication`,
/// `call_multiplication`, `call_linked`), the dedicated `cyp2d6` genotyper, and
/// several per-gene match arms. Each expected value is PyPGx 0.26.0 truth for the
/// same (gene, Haplotype1/2, VariantData, CNV) inputs, captured in `.refenv`.
#[test]
fn sv_genotyper_helpers_match_pypgx() {
    // (gene, haplotype1, haplotype2, variant_data, cnv, expected genotype)
    let cases: &[(&str, &str, &str, &str, &str, &str)] = &[
        ("CYP2A6", "*2;*1;", "*1;", "*1:default;*2:19-x-A-G:0.8;", "WholeDup1", "*1/*2x2"),
        ("CYP2E1", "*2;*1;", "*1;", "*1:default;*2:10-x-C-T:0.8;", "WholeMultip1", "*1/*2x3"),
        ("UGT1A4", "*1;", "*1;", "*1:default;", "NoncodingDel1", "*1/*S1"),
        ("SLC22A2", "*K432Q;*1;", "*1;", "*1:default;*K432Q:6-x-G-A:0.9;", "NoncodingDel1", "*1/*S1"),
        ("GSTM1", "*A;", "*A;", "*A:default;", "WholeDel1", "*0/*A"),
        ("GSTM1", "*A;", "*A;", "*A:default;", "WholeDel1Hom", "*0/*0"),
        ("SULT1A1", "*1;", "*1;", "*1:default;", "WholeDel1Hom", "*DEL/*DEL"),
        ("CYP2D6", "*4;*1;", "*1;", "*1:default;*4:22-x-C-T:0.9;", "WholeDel1", "*4/*5"),
        ("CYP2D6", "*4;*1;", "*1;", "*1:default;*4:22-x-C-T:0.9;", "Normal", "*1/*4"),
        ("UGT2B17", "*1;", "*1;", "*1:default;", "WholeDel1", "*1/*2"),
        ("CYP2E1", "*7;*1;", "*1;", "*1:default;*7:10-x-G-C:0.9;", "PartialDup1", "*1/*S1"),
    ];
    for &(gene, h1, h2, vd, cnv_call, expected) in cases {
        let alleles = table(
            &[("Gene", gene), ("Assembly", "GRCh37"), ("SemanticType", "SampleTable[Alleles]")],
            &["Haplotype1", "Haplotype2", "AlternativePhase", "VariantData"],
            &["S1"],
            vec![vec![h1, h2, ";", vd]],
        );
        let cnv = table(
            &[("Gene", gene), ("Assembly", "GRCh37"), ("SemanticType", "SampleTable[CNVCalls]")],
            &["CNV"],
            &["S1"],
            vec![vec![cnv_call]],
        );
        let gt = pypgx::call_genotypes(Some(&alleles), Some(&cnv))
            .unwrap_or_else(|e| panic!("{gene} {cnv_call}: {e}"));
        assert_eq!(
            gt.as_sample_table().loc("S1"),
            &vec![expected.to_string()],
            "{gene} {cnv_call} should match PyPGx",
        );
    }
}
