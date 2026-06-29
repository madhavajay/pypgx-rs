use pypgx::fuc::VcfFrame;
use pypgx::sdk::{Archive, ArchiveData};

#[test]
fn g6pd_hg01621_haplotype_columns_match_pypgx_phase_orientation() {
    let columns: Vec<String> = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "HG01621",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    let rows = vec![vec![
        "X".to_string(),
        "154533596".to_string(),
        ".".to_string(),
        "C".to_string(),
        "G".to_string(),
        ".".to_string(),
        ".".to_string(),
        "Phased".to_string(),
        "GT:AD:DP:AF".to_string(),
        "1|0:10,10:20:0.5".to_string(),
    ]];
    let consolidated = Archive::new(
        vec![
            ("Platform".to_string(), "WGS".to_string()),
            ("Gene".to_string(), "G6PD".to_string()),
            ("Assembly".to_string(), "GRCh38".to_string()),
            (
                "SemanticType".to_string(),
                "VcfFrame[Consolidated]".to_string(),
            ),
        ],
        ArchiveData::Vcf(VcfFrame::new(Vec::new(), columns, rows)),
    );

    let alleles = pypgx::predict_alleles(&consolidated).unwrap();
    let table = alleles.as_sample_table();
    let row = table.loc("HG01621");
    let h1 = table.columns.iter().position(|c| c == "Haplotype1").unwrap();
    let h2 = table.columns.iter().position(|c| c == "Haplotype2").unwrap();

    assert_eq!(
        row[h1],
        "Seattle, Lodi, Modena, Ferrara II, Athens-like;"
    );
    assert_eq!(row[h2], "B (reference);");
}
