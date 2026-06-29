//! End-to-end test for `run_ngs_pipeline` on the native path: pre-phased input
//! (→ no Beagle) and no depth-of-coverage (→ no CNV/sklearn). Verifies the
//! Results table vs Python. Without the `beagle` feature, the Beagle arm is
//! exercised structurally and surfaces `NotPorted`.

use pypgx::fuc::VcfFrame;
use serde_json::Value;

const PHASED_VCF: &str = include_str!("fixtures/ngs_phased.vcf");
const TRUTH: &str = include_str!("fixtures/ngs_pipeline.json");

fn canonical_variant_data(s: &str) -> Vec<(String, Vec<(String, String)>)> {
    let mut out = Vec::new();
    for seg in s.split(';').filter(|x| !x.is_empty()) {
        let parts: Vec<&str> = seg.split(':').collect();
        let mut pairs: Vec<(String, String)> = if parts.len() == 3 {
            parts[1]
                .split(',')
                .zip(parts[2].split(','))
                .map(|(v, f)| (v.to_string(), f.to_string()))
                .collect()
        } else {
            vec![(parts.get(1).copied().unwrap_or("").to_string(), String::new())]
        };
        pairs.sort();
        out.push((parts[0].to_string(), pairs));
    }
    out.sort();
    out
}

#[test]
fn run_ngs_pipeline_native_path_matches_python() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();
    let vf = VcfFrame::from_string(PHASED_VCF);
    let out = format!("{}/pypgx_ngs_test", std::env::temp_dir().display());
    std::fs::remove_dir_all(&out).ok();

    pypgx::run_ngs_pipeline(
        "CYP4F2", &out, Some(&vf), None, None, "WGS", "GRCh37", None, false, None, false, None,
        None,
    )
    .unwrap();

    let results = pypgx::Archive::from_file(&format!("{out}/results.zip")).unwrap();
    let st = results.as_sample_table();
    let expected: Vec<String> = t["rows"][0]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_string())
        .collect();
    let got = st.loc("A");

    assert_eq!(&got[..5], &expected[..5], "Genotype..AlternativePhase");
    assert_eq!(
        canonical_variant_data(&got[5]),
        canonical_variant_data(&expected[5]),
        "VariantData (set-equal)"
    );
    assert!(got[6].is_empty() || got[6] == "nan");
}

/// Without the `beagle` feature, unphased NGS input surfaces `NotPorted`.
#[cfg(not(feature = "beagle"))]
#[test]
fn run_ngs_pipeline_beagle_arm_is_notported() {
    // An unphased input forces the statistical-phasing branch.
    let unphased = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
                    19\t16008388\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,20:30\n";
    let vf = VcfFrame::from_string(unphased);
    let out = format!("{}/pypgx_ngs_beagle", std::env::temp_dir().display());
    std::fs::remove_dir_all(&out).ok();

    let err = pypgx::run_ngs_pipeline(
        "CYP4F2", &out, Some(&vf), None, None, "WGS", "GRCh37", None, false, None, false, None,
        None,
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("not yet ported") && err.to_string().contains("beagle"),
        "expected Beagle NotPorted, got: {err}"
    );
}
