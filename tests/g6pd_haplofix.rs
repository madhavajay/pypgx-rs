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

#[test]
fn g6pd_hg02614_variant_data_is_coordinate_sorted_and_semantically_pypgx_equivalent() {
    let columns: Vec<String> = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "HG02614",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    let rows = vec![
        vec![
            "X".to_string(),
            "154535342".to_string(),
            ".".to_string(),
            "C".to_string(),
            "T".to_string(),
            ".".to_string(),
            ".".to_string(),
            "Phased".to_string(),
            "GT:AD:DP:AF".to_string(),
            "1|0:10,10:20:0.0,0.34".to_string(),
        ],
        vec![
            "X".to_string(),
            "154535277".to_string(),
            ".".to_string(),
            "T".to_string(),
            "C".to_string(),
            ".".to_string(),
            ".".to_string(),
            "Phased".to_string(),
            "GT:AD:DP:AF".to_string(),
            "1|1:10,10:20:0.0,0.27".to_string(),
        ],
    ];
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
    let row = table.loc("HG02614");
    let col = |name: &str| table.columns.iter().position(|c| c == name).unwrap();

    assert_eq!(row[col("Haplotype1")], "Sierra Leone;");
    assert_eq!(row[col("Haplotype2")], "A;");
    assert_eq!(row[col("AlternativePhase")], ";");
    assert_eq!(
        row[col("VariantData")],
        "Sierra Leone:X-154535277-T-C,X-154535342-C-T:0.27,0.34;A:X-154535277-T-C:0.27;"
    );

    let python_observed_set_order =
        "Sierra Leone:X-154535342-C-T,X-154535277-T-C:0.34,0.27;A:X-154535277-T-C:0.27;";
    assert_eq!(
        canonical_variant_data(&row[col("VariantData")]),
        canonical_variant_data(python_observed_set_order),
        "VariantData differs only by Python set iteration order"
    );
}

fn canonical_variant_data(s: &str) -> Vec<(String, Vec<(String, String)>)> {
    let mut out = Vec::new();
    for segment in s.split(';').filter(|x| !x.is_empty()) {
        let parts: Vec<&str> = segment.split(':').collect();
        let allele = parts[0].to_string();
        let mut pairs: Vec<(String, String)> = if parts.len() == 3 {
            parts[1]
                .split(',')
                .zip(parts[2].split(','))
                .map(|(variant, fraction)| (variant.to_string(), fraction.to_string()))
                .collect()
        } else {
            vec![(
                parts.get(1).copied().unwrap_or("").to_string(),
                String::new(),
            )]
        };
        pairs.sort();
        out.push((allele, pairs));
    }
    out.sort();
    out
}
