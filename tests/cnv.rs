//! Verifies the native RBF OvR-SVM decision function reproduces sklearn's
//! `OneVsRestClassifier(SVC).predict` exactly. Reference model + predictions
//! generated from scikit-learn 1.7.2 in `.refenv` (`tests/fixtures/cnv_svm.json`).

use pypgx::cnv::{median_filter, CnvEstimator, CnvModel};
use serde_json::Value;

const TRUTH: &str = include_str!("fixtures/cnv_svm.json");
const MEDIAN: &str = include_str!("fixtures/cnv_median.json");

#[test]
fn median_filter_matches_scipy() {
    let cases: Value = serde_json::from_str(MEDIAN).unwrap();
    for case in cases.as_array().unwrap() {
        let input = f64s(&case["input"]);
        let size = case["size"].as_u64().unwrap() as usize;
        let expected = f64s(&case["out"]);
        let got = median_filter(&input, size);
        assert_eq!(got.len(), expected.len());
        for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
            assert!((g - e).abs() < 1e-9, "size {size} idx {i}: {g} vs scipy {e}");
        }
    }
}

fn f64s(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}

fn model(t: &Value) -> CnvModel {
    let classes = t["classes"].as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect();
    let estimators = t["estimators"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| CnvEstimator {
            label: e["label"].as_i64().unwrap(),
            gamma: e["gamma"].as_f64().unwrap(),
            intercept: e["intercept"].as_f64().unwrap(),
            dual_coef: f64s(&e["dual_coef"]),
            support_vectors: e["support_vectors"]
                .as_array()
                .unwrap()
                .iter()
                .map(f64s)
                .collect(),
        })
        .collect();
    CnvModel { classes, estimators }
}

/// End-to-end: `predict_cnv` (process_copy_number at the real CYP2A6 region
/// dimension → RBF predict → code→name), `Model[CNV]` archive round-trip
/// (serde JSON), and `test_cnv_caller` accuracy.
#[test]
fn predict_cnv_and_test_end_to_end() {
    use pypgx::cnv::{CnvEstimator, CnvModel};
    use pypgx::fuc::CovFrame;
    use pypgx::sdk::{Archive, ArchiveData, SampleTable};

    let (start, end) = (41_339_442i64, 41_396_352i64); // CYP2A6 GRCh37
    let n = (end - start + 1) as usize;

    // Dense copy-number profile (sample A ≈ CN 2 everywhere).
    let rows: Vec<Vec<String>> = (start..=end)
        .map(|p| vec!["19".to_string(), p.to_string(), "2.0".to_string()])
        .collect();
    let cn = Archive::new(
        vec![
            ("Gene".to_string(), "CYP2A6".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "CovFrame[CopyNumber]".to_string()),
            ("Platform".to_string(), "WGS".to_string()),
        ],
        ArchiveData::Cov(CovFrame {
            columns: vec!["Chromosome".into(), "Position".into(), "A".into()],
            rows,
        }),
    );

    // Model that scores class 0 ("Normal") far above class 1 for x≈2.
    let model = CnvModel {
        classes: vec![0, 1],
        estimators: vec![
            CnvEstimator { label: 0, gamma: 0.001, intercept: 0.5, dual_coef: vec![1.0], support_vectors: vec![vec![2.0; n]] },
            CnvEstimator { label: 1, gamma: 0.001, intercept: 0.0, dual_coef: vec![1.0], support_vectors: vec![vec![4.0; n]] },
        ],
    };
    let model_archive = Archive::new(
        vec![
            ("Gene".to_string(), "CYP2A6".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "Model[CNV]".to_string()),
        ],
        ArchiveData::Model(model),
    );

    // Model[CNV] archive round-trips through the zip (data.json).
    let tmp = format!("{}/pypgx_cnv_model.zip", std::env::temp_dir().display());
    model_archive.to_file(&tmp).unwrap();
    let model_back = Archive::from_file(&tmp).unwrap();
    assert_eq!(model_back.semantic_type(), "Model[CNV]");

    let calls = pypgx::predict_cnv(&cn, Some(&model_back)).unwrap();
    assert_eq!(calls.semantic_type(), "SampleTable[CNVCalls]");
    assert_eq!(calls.as_sample_table().loc("A"), &vec!["Normal".to_string()]);

    let cnv_calls = Archive::new(
        vec![
            ("Gene".to_string(), "CYP2A6".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "SampleTable[CNVCalls]".to_string()),
        ],
        ArchiveData::SampleTable(SampleTable {
            index: vec!["A".to_string()],
            columns: vec!["CNV".to_string()],
            rows: vec![vec!["Normal".to_string()]],
        }),
    );
    let report = pypgx::test_cnv_caller(&model_back, &cn, &cnv_calls).unwrap();
    assert_eq!((report.accuracy, report.correct, report.total), (1.0, 1, 1));
}

/// `predict_cnv(cnv_caller=None)` resolves the default model from the bundle
/// (`$PYPGX_BUNDLE/cnv/{assembly}/{gene}.zip`) — the path the NGS pipeline uses.
#[test]
fn predict_cnv_resolves_default_model_from_bundle() {
    use pypgx::cnv::{CnvEstimator, CnvModel};
    use pypgx::fuc::CovFrame;
    use pypgx::sdk::{Archive, ArchiveData};

    let (start, end) = (41_339_442i64, 41_396_352i64); // CYP2A6 GRCh37
    let n = (end - start + 1) as usize;

    // A bundle laid out like pypgx-bundle: {bundle}/cnv/GRCh37/CYP2A6.zip.
    let bundle = std::env::temp_dir().join(format!("pypgx_bundle_test_{}", std::process::id()));
    let cnv_dir = bundle.join("cnv").join("GRCh37");
    std::fs::create_dir_all(&cnv_dir).unwrap();
    let model = CnvModel {
        classes: vec![0, 1],
        estimators: vec![
            CnvEstimator { label: 0, gamma: 0.001, intercept: 0.5, dual_coef: vec![1.0], support_vectors: vec![vec![2.0; n]] },
            CnvEstimator { label: 1, gamma: 0.001, intercept: 0.0, dual_coef: vec![1.0], support_vectors: vec![vec![4.0; n]] },
        ],
    };
    Archive::new(
        vec![
            ("Gene".to_string(), "CYP2A6".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "Model[CNV]".to_string()),
        ],
        ArchiveData::Model(model),
    )
    .to_file(cnv_dir.join("CYP2A6.zip").to_str().unwrap())
    .unwrap();

    let cn = Archive::new(
        vec![
            ("Gene".to_string(), "CYP2A6".to_string()),
            ("Assembly".to_string(), "GRCh37".to_string()),
            ("SemanticType".to_string(), "CovFrame[CopyNumber]".to_string()),
            ("Platform".to_string(), "WGS".to_string()),
        ],
        ArchiveData::Cov(CovFrame {
            columns: vec!["Chromosome".into(), "Position".into(), "A".into()],
            rows: (start..=end).map(|p| vec!["19".into(), p.to_string(), "2.0".into()]).collect(),
        }),
    );

    std::env::set_var("PYPGX_BUNDLE", &bundle);
    // No explicit caller → must load {bundle}/cnv/GRCh37/CYP2A6.zip and predict.
    let calls = pypgx::predict_cnv(&cn, None).expect("default model resolved from bundle");
    assert_eq!(calls.as_sample_table().loc("A"), &vec!["Normal".to_string()]);
    std::env::remove_var("PYPGX_BUNDLE");
    let _ = std::fs::remove_dir_all(&bundle);
}

#[test]
fn rbf_ovr_predict_matches_sklearn() {
    let t: Value = serde_json::from_str(TRUTH).unwrap();
    let m = model(&t);
    let x_test: Vec<Vec<f64>> = t["X_test"].as_array().unwrap().iter().map(f64s).collect();
    let expected: Vec<i64> = t["predictions"].as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect();

    let got: Vec<i64> = x_test.iter().map(|x| m.predict(x)).collect();
    assert_eq!(got, expected, "OvR-SVM predictions must match sklearn");

    // Per-estimator decision_function scores must match sklearn too (tight tol).
    for (k, est) in m.estimators.iter().enumerate() {
        let sk_scores = f64s(&t["scores"][k]);
        for (i, x) in x_test.iter().enumerate() {
            let d = est.decision(x);
            assert!(
                (d - sk_scores[i]).abs() < 1e-9,
                "estimator {k} sample {i}: {d} vs sklearn {}",
                sk_scores[i]
            );
        }
    }
}
