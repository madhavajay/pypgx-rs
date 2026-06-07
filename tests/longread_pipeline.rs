//! End-to-end test for `run_long_read_pipeline` vs the Python reference.
//! Runs the full chain (import_variants[LongRead] → predict_alleles →
//! call_genotypes → call_phenotypes → combine_results), writes the archives,
//! reads back `results.zip`, and compares the Results table. Also stress-tests
//! multi-variant `predict_alleles` VariantData ordering (the *4 allele).

use pypgx::fuc::VcfFrame;
use serde_json::Value;

const PE_VCF: &str = include_str!("fixtures/phase_ext_input.vcf");
const TRUTH: &str = include_str!("fixtures/longread_pipeline.json");

#[test]
fn run_long_read_pipeline_matches_python() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();
    let vf = VcfFrame::from_string(PE_VCF);

    let out = format!("{}/pypgx_lrp_test", std::env::temp_dir().display());
    std::fs::remove_dir_all(&out).ok();
    pypgx::run_long_read_pipeline("CYP4F2", &out, &vf, "GRCh37", false, None, false).unwrap();

    // All five intermediate archives written.
    for f in [
        "consolidated-variants.zip",
        "alleles.zip",
        "genotypes.zip",
        "phenotypes.zip",
        "results.zip",
    ] {
        assert!(
            std::path::Path::new(&format!("{out}/{f}")).exists(),
            "missing {f}"
        );
    }

    let results = pypgx::Archive::from_file(&format!("{out}/results.zip")).unwrap();
    let st = results.as_sample_table();

    let cols: Vec<String> = t["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
    assert_eq!(st.columns, cols, "columns");
    assert_eq!(st.index, vec!["A".to_string()], "index");

    let expected: Vec<String> = t["rows"][0]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
    let got = st.loc("A");

    // Columns 0..5 (Genotype, Phenotype, Haplotype1/2, AlternativePhase) match
    // byte-for-byte.
    assert_eq!(&got[..5], &expected[..5], "results row (Genotype..AlternativePhase)");

    // Column 5 (VariantData): for multi-variant alleles (e.g. *4), Python joins
    // a CPython `set` (hash order) while Rust uses a deterministic order — the
    // variant↔fraction PAIRS are identical, only the comma-order differs (a
    // documented, unavoidable divergence; see TODO §12). Compare semantically.
    assert_eq!(
        canonical_variant_data(&got[5]),
        canonical_variant_data(&expected[5]),
        "VariantData (set-equal): got {:?} vs {:?}",
        got[5],
        expected[5]
    );

    // CNV: empty here. Python re-reads the empty TSV cell as NaN ("nan"); the
    // Rust reader keeps it as "". Same serialized archive, different rendering.
    assert!(got[6].is_empty() || got[6] == "nan", "CNV cell: {:?}", got[6]);
}

/// Canonicalize a VariantData string ("allele:variants:fractions;…") by sorting
/// each allele's (variant, fraction) pairs, so a set-order difference between
/// CPython and Rust does not cause a spurious mismatch.
fn canonical_variant_data(s: &str) -> Vec<(String, Vec<(String, String)>)> {
    let mut out = Vec::new();
    for seg in s.split(';').filter(|x| !x.is_empty()) {
        let parts: Vec<&str> = seg.split(':').collect();
        let allele = parts[0].to_string();
        let mut pairs: Vec<(String, String)> = if parts.len() == 3 {
            parts[1]
                .split(',')
                .zip(parts[2].split(','))
                .map(|(v, f)| (v.to_string(), f.to_string()))
                .collect()
        } else {
            // e.g. "allele:default"
            vec![(parts.get(1).copied().unwrap_or("").to_string(), String::new())]
        };
        pairs.sort();
        out.push((allele, pairs));
    }
    out.sort();
    out
}
