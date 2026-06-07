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
