//! Breadth parity: assert the additional pure `core` functions match the
//! Python reference (captured in `tests/fixtures/truth2.json`).

use pypgx::core;
use pypgx::fuc::python_float_str;
use serde_json::Value;

const TRUTH2: &str = include_str!("fixtures/truth2.json");

fn truth() -> Value {
    serde_json::from_str(TRUTH2).expect("parse truth2.json")
}

fn ints(v: &Value) -> Vec<i64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_i64().unwrap())
        .collect()
}

fn strs(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn predict_score_matches() {
    let t = truth();
    for (k, v) in t["predict_score"].as_object().unwrap() {
        let (gene, allele) = k.split_once('|').unwrap();
        let got = python_float_str(core::predict_score(gene, allele));
        assert_eq!(got, v.as_str().unwrap(), "predict_score {k}");
    }
}

#[test]
fn predict_phenotype_matches() {
    let t = truth();
    for (k, v) in t["predict_phenotype"].as_object().unwrap() {
        let parts: Vec<&str> = k.split('|').collect();
        let got = core::predict_phenotype(parts[0], parts[1], parts[2]);
        assert_eq!(got, v.as_str().unwrap(), "predict_phenotype {k}");
    }
}

#[test]
fn get_priority_matches() {
    let t = truth();
    for (k, v) in t["get_priority"].as_object().unwrap() {
        let (gene, ph) = k.split_once('|').unwrap();
        let got = core::get_priority(gene, ph).unwrap();
        assert_eq!(got, v.as_str().unwrap(), "get_priority {k}");
    }
}

#[test]
fn get_region_matches() {
    let t = truth();
    for (k, v) in t["get_region"].as_object().unwrap() {
        let (gene, asm) = k.split_once('|').unwrap();
        assert_eq!(
            core::get_region(gene, asm).unwrap(),
            v.as_str().unwrap(),
            "get_region {k}"
        );
    }
}

#[test]
fn get_exons_match() {
    let t = truth();
    for (k, v) in t["get_exon_starts"].as_object().unwrap() {
        let (gene, asm) = k.split_once('|').unwrap();
        assert_eq!(
            core::get_exon_starts(gene, asm).unwrap(),
            ints(v),
            "exon_starts {k}"
        );
    }
    for (k, v) in t["get_exon_ends"].as_object().unwrap() {
        let (gene, asm) = k.split_once('|').unwrap();
        assert_eq!(
            core::get_exon_ends(gene, asm).unwrap(),
            ints(v),
            "exon_ends {k}"
        );
    }
}

#[test]
fn get_strand_and_paralog_match() {
    let t = truth();
    for (g, v) in t["get_strand"].as_object().unwrap() {
        assert_eq!(
            core::get_strand(g).unwrap(),
            v.as_str().unwrap(),
            "strand {g}"
        );
    }
    for (g, v) in t["get_paralog"].as_object().unwrap() {
        assert_eq!(core::get_paralog(g), v.as_str().unwrap(), "paralog {g}");
    }
}

#[test]
fn list_functions_match() {
    let t = truth();
    for (g, v) in t["list_functions"].as_object().unwrap() {
        let got: Vec<String> = core::list_functions(Some(g))
            .into_iter()
            .map(|x| x.unwrap_or_else(|| "nan".to_string()))
            .collect();
        assert_eq!(got, strs(v), "list_functions {g}");
    }
    let got_all: Vec<String> = core::list_functions(None)
        .into_iter()
        .map(|x| x.unwrap_or_else(|| "nan".to_string()))
        .collect();
    assert_eq!(got_all, strs(&t["list_functions_all"]));
}

#[test]
fn list_phenotypes_match() {
    let t = truth();
    for (g, v) in t["list_phenotypes"].as_object().unwrap() {
        assert_eq!(
            core::list_phenotypes(Some(g)),
            strs(v),
            "list_phenotypes {g}"
        );
    }
    assert_eq!(core::list_phenotypes(None), strs(&t["list_phenotypes_all"]));
}

#[test]
fn scalar_predicates_match() {
    let t = truth();
    for (k, v) in t["is_legit_allele"].as_object().unwrap() {
        let (gene, allele) = k.split_once('|').unwrap();
        assert_eq!(
            core::is_legit_allele(gene, allele),
            v.as_bool().unwrap(),
            "is_legit {k}"
        );
    }
    for (g, v) in t["has_score"].as_object().unwrap() {
        assert_eq!(
            core::has_score(g).unwrap(),
            v.as_bool().unwrap(),
            "has_score {g}"
        );
    }
    for (g, v) in t["has_sv"].as_object().unwrap() {
        assert_eq!(core::has_sv(g).unwrap(), v.as_bool().unwrap(), "has_sv {g}");
    }
    for (k, v) in t["get_function"].as_object().unwrap() {
        let (gene, allele) = k.split_once('|').unwrap();
        let got = core::get_function(gene, allele)
            .unwrap()
            .unwrap_or_else(|| "nan".to_string());
        assert_eq!(got, v.as_str().unwrap(), "get_function {k}");
    }
    for (k, v) in t["get_score"].as_object().unwrap() {
        let (gene, allele) = k.split_once('|').unwrap();
        let got = match core::get_score(gene, allele).unwrap() {
            Some(f) => python_float_str(f),
            None => "nan".to_string(),
        };
        assert_eq!(got, v.as_str().unwrap(), "get_score {k}");
    }
}

#[test]
fn get_recommendation_matches() {
    let t = truth();
    for (k, v) in t["get_recommendation"].as_object().unwrap() {
        let p: Vec<&str> = k.split('|').collect();
        let (gene2, phenotype2) = if p.len() == 5 {
            (Some(p[3]), Some(p[4]))
        } else {
            (None, None)
        };
        let got = core::get_recommendation(p[0], p[1], p[2], gene2, phenotype2).unwrap();
        assert_eq!(got, v.as_str().unwrap(), "get_recommendation {k}");
    }
}
