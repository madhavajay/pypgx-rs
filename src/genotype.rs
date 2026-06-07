//! Port of `pypgx.api.genotype` — final diplotype calling from candidate star
//! alleles and (optionally) CNV calls. `call_genotypes` dispatches to a
//! per-gene genotyper; genes without SV fall through to the simple path.
//!
//! The genotypers mirror PyPGx's explicit `if/elif` ladders branch-for-branch,
//! so several arms are deliberately identical (e.g. different CNV states that
//! happen to yield the same diplotype). We keep them separate for a faithful
//! 1-for-1 port rather than collapsing them.
#![allow(clippy::if_same_then_else)]

use std::collections::HashMap;

use crate::core;
use crate::sdk::{Archive, ArchiveData, PgxError, SampleTable};

/// Parsed `VariantData` entry for one allele: `Empty` mirrors the Python `[]`
/// (a default allele, falsy), `Data` holds the allele-fraction list (PyPGx's
/// `VariantData[allele][1]`; the variant names at `[0]` are unused downstream).
#[derive(Clone, Debug)]
enum VData {
    Empty,
    Data(Vec<f64>),
}

impl VData {
    fn truthy(&self) -> bool {
        matches!(self, VData::Data(..))
    }
    fn fractions(&self) -> &[f64] {
        match self {
            VData::Data(f) => f,
            VData::Empty => &[],
        }
    }
}

/// One parsed sample row.
struct GtRow {
    h1: Vec<String>,
    h2: Vec<String>,
    vd: HashMap<String, VData>,
    cnv: Option<String>,
}

impl GtRow {
    fn a1(&self) -> &str {
        &self.h1[0]
    }
    fn a2(&self) -> &str {
        &self.h2[0]
    }
    fn vd(&self, allele: &str) -> &VData {
        self.vd.get(allele).unwrap_or(&VData::Empty)
    }
    fn h1_has(&self, allele: &str) -> bool {
        self.h1.iter().any(|x| x == allele)
    }
    fn h2_has(&self, allele: &str) -> bool {
        self.h2.iter().any(|x| x == allele)
    }
}

fn split_semis(s: &str) -> Vec<String> {
    s.trim_matches(';')
        .split(';')
        .map(|x| x.to_string())
        .collect()
}

fn parse_variant_data(s: &str) -> HashMap<String, VData> {
    let mut d = HashMap::new();
    for allele in split_semis(s) {
        if allele.is_empty() {
            continue;
        }
        let fields: Vec<&str> = allele.split(':').collect();
        if allele.contains("default") {
            d.insert(fields[0].to_string(), VData::Empty);
        } else {
            // fields = [allele, variants, fractions]; only fractions are used.
            let fracs: Vec<f64> = fields[2]
                .split(',')
                .map(|x| x.parse::<f64>().expect("fraction"))
                .collect();
            d.insert(fields[0].to_string(), VData::Data(fracs));
        }
    }
    d
}

/// Genes with a dedicated SV-aware genotyper (`sv_genotypers` in PyPGx).
const SV_GENES: &[&str] = &[
    "CYP2A6", "CYP2B6", "CYP2D6", "CYP2E1", "CYP4F2", "G6PD", "GSTM1", "GSTT1", "SLC22A2",
    "SULT1A1", "UGT1A4", "UGT2B15", "UGT2B17",
];

/// `call_genotypes(alleles, cnv_calls)` → `SampleTable[Genotypes]`.
pub fn call_genotypes(
    alleles: Option<&Archive>,
    cnv_calls: Option<&Archive>,
) -> Result<Archive, PgxError> {
    if let Some(a) = alleles {
        a.check_type(&["SampleTable[Alleles]"])?;
    }
    if let Some(c) = cnv_calls {
        c.check_type(&["SampleTable[CNVCalls]"])?;
    }

    let (index, rows, gene, assembly): (Vec<String>, Vec<GtRow>, String, String) =
        match (alleles, cnv_calls) {
            (Some(a), Some(c)) => {
                let at = a.as_sample_table();
                let ct = c.as_sample_table();
                let aset: std::collections::HashSet<&String> = at.index.iter().collect();
                let cset: std::collections::HashSet<&String> = ct.index.iter().collect();
                if aset != cset {
                    return Err(PgxError::IncorrectMetadata(
                        "SampleTable[Alleles] and SampleTable[CNVCalls] have different samples".into(),
                    ));
                }
                if a.get("Gene") != c.get("Gene") {
                    return Err(PgxError::IncorrectMetadata(
                        "Found two different target genes".into(),
                    ));
                }
                let cnv_idx = ct.columns.iter().position(|x| x == "CNV").ok_or_else(|| {
                    PgxError::IncorrectMetadata("SampleTable[CNVCalls] missing 'CNV' column".into())
                })?;
                let mut rows = Vec::new();
                for (i, sample) in at.index.iter().enumerate() {
                    let cj = ct.index.iter().position(|x| x == sample).ok_or_else(|| {
                        PgxError::External(format!("sample {sample} missing from CNV calls"))
                    })?;
                    rows.push(parse_alleles_row(at, i, Some(ct.rows[cj][cnv_idx].clone())));
                }
                (
                    at.index.clone(),
                    rows,
                    a.get("Gene").unwrap().to_string(),
                    a.get("Assembly").unwrap().to_string(),
                )
            }
            (Some(a), None) => {
                let at = a.as_sample_table();
                let rows = (0..at.index.len())
                    .map(|i| parse_alleles_row(at, i, None))
                    .collect();
                (
                    at.index.clone(),
                    rows,
                    a.get("Gene").unwrap().to_string(),
                    a.get("Assembly").unwrap().to_string(),
                )
            }
            (None, Some(c)) => {
                let ct = c.as_sample_table();
                let cnv_idx = ct.columns.iter().position(|x| x == "CNV").ok_or_else(|| {
                    PgxError::IncorrectMetadata("SampleTable[CNVCalls] missing 'CNV' column".into())
                })?;
                let rows = (0..ct.index.len())
                    .map(|i| GtRow {
                        h1: Vec::new(),
                        h2: Vec::new(),
                        vd: HashMap::new(),
                        cnv: Some(ct.rows[i][cnv_idx].clone()),
                    })
                    .collect();
                (
                    ct.index.clone(),
                    rows,
                    c.get("Gene").unwrap().to_string(),
                    c.get("Assembly").unwrap().to_string(),
                )
            }
            (None, None) => {
                return Err(PgxError::IncorrectMetadata(
                    "Either SampleTable[Alleles] or SampleTable[CNVCalls] must be provided".into(),
                ));
            }
        };

    let is_sv = SV_GENES.contains(&gene.as_str());
    let mut out_rows = Vec::new();
    for mut row in rows {
        if is_sv && row.cnv.is_none() {
            // PyPGx warns and assumes no SV.
            row.cnv = Some("AssumeNormal".to_string());
        }
        let result = genotype_row(&gene, &assembly, &row);
        let sorted = core::sort_alleles(&result, "name", None, &assembly);
        out_rows.push(vec![sorted.join("/")]);
    }

    let metadata = vec![
        ("Gene".to_string(), gene),
        ("Assembly".to_string(), assembly),
        (
            "SemanticType".to_string(),
            "SampleTable[Genotypes]".to_string(),
        ),
    ];
    let table = SampleTable {
        index,
        columns: vec!["Genotype".to_string()],
        rows: out_rows,
    };
    Ok(Archive::new(metadata, ArchiveData::SampleTable(table)))
}

fn parse_alleles_row(at: &SampleTable, i: usize, cnv: Option<String>) -> GtRow {
    let col = |name: &str| at.columns.iter().position(|c| c == name).unwrap();
    let h1 = split_semis(&at.rows[i][col("Haplotype1")]);
    let h2 = split_semis(&at.rows[i][col("Haplotype2")]);
    let vd = parse_variant_data(&at.rows[i][col("VariantData")]);
    GtRow { h1, h2, vd, cnv }
}

fn s(x: &str) -> String {
    x.to_string()
}

/// Dispatch to the per-gene genotyper, returning the (pre-name-sort) result.
fn genotype_row(gene: &str, assembly: &str, r: &GtRow) -> Vec<String> {
    let cnv = r.cnv.as_deref().unwrap_or("");
    let priority =
        |alleles: &[String]| core::sort_alleles(alleles, "priority", Some(gene), assembly);
    match gene {
        "CYP2A6" => {
            let (a1, a2) = (r.a1(), r.a2());
            let sp = priority(&[s(a1), s(a2)]);
            let s1 = &sp[0];
            match cnv {
                "Normal" | "AssumeNormal" | "ParalogWholeDup1" | "ParalogWholeDel1" => {
                    vec![s(a1), s(a2)]
                }
                "WholeDel1" | "WholeDel2" | "WholeDel3" => vec![s1.clone(), s("*4")],
                "WholeDel1Hom" | "WholeDel2Hom" => vec![s("*4"), s("*4")],
                "Hybrid2" => vec![s1.clone(), s("*12")],
                "Hybrid2Hom" => vec![s("*12"), s("*12")],
                "Hybrid3" => vec![s1.clone(), s("*34")],
                "WholeDup1" | "WholeDup2" | "WholeDup3" => call_duplication(r),
                _ => vec![s("Indeterminate")],
            }
        }
        "CYP2B6" => {
            let (a1, a2) = (r.a1(), r.a2());
            let p = priority(&[s(a1), s(a2)])[0].clone();
            match cnv {
                "Normal" | "AssumeNormal" => vec![s(a1), s(a2)],
                "Hybrid1" => vec![p, s("*29")],
                "WholeDup1" => call_duplication(r),
                _ => vec![s("Indeterminate")],
            }
        }
        "CYP2D6" => cyp2d6(assembly, r),
        "CYP2E1" => {
            let (a1, a2) = (r.a1(), r.a2());
            match cnv {
                "Normal" | "AssumeNormal" => vec![s(a1), s(a2)],
                "PartialDup1" => call_linked(r, "*7", "*S1"),
                "PartialDup1Hom" => vec![s("*S1"), s("*S1")],
                "WholeDup1" | "WholeDup2" => call_duplication(r),
                "WholeMultip1" => call_multiplication(r),
                _ => vec![s("Indeterminate")],
            }
        }
        "CYP4F2" => {
            let (a1, a2) = (r.a1(), r.a2());
            let s1 = priority(&[s(a1), s(a2)])[0].clone();
            match cnv {
                "Normal" | "AssumeNormal" => vec![s(a1), s(a2)],
                "WholeDel1" => vec![s1, s("*DEL")],
                _ => vec![s("Indeterminate")],
            }
        }
        "G6PD" => {
            let (a1, a2) = (r.a1(), r.a2());
            let s1 = priority(&[s(a1), s(a2)])[0].clone();
            match cnv {
                "Female" | "AssumeNormal" => vec![s(a1), s(a2)],
                "Male" => vec![s1, s("MALE")],
                _ => vec![s("Indeterminate")],
            }
        }
        "GSTM1" => {
            let (a1, a2) = (r.a1(), r.a2());
            let s1 = priority(&[s(a1), s(a2)])[0].clone();
            match cnv {
                "Normal" | "AssumeNormal" | "NoncodingDel1" => vec![s(a1), s(a2)],
                "WholeDel1" | "WholeDel1+NoncodingDel1" | "WholeDel2" => vec![s1, s("*0")],
                "WholeDel1Hom" | "WholeDel1+WholeDel2" => vec![s("*0"), s("*0")],
                "WholeDup1" => call_duplication(r),
                _ => vec![s("Indeterminate")],
            }
        }
        "GSTT1" => match cnv {
            "Normal" | "AssumeNormal" => vec![s("*A"), s("*A")],
            "WholeDel1" => vec![s("*A"), s("*0")],
            "WholeDel1Hom" => vec![s("*0"), s("*0")],
            _ => vec![s("Indeterminate")],
        },
        "SLC22A2" => slc22a2(r),
        "SULT1A1" => {
            let (a1, a2) = (r.a1(), r.a2());
            let s1 = priority(&[s(a1), s(a2)])[0].clone();
            match cnv {
                "Normal" | "AssumeNormal" => vec![s(a1), s(a2)],
                "WholeDel1" => vec![s1, s("*DEL")],
                "WholeDel1Hom" => vec![s("*DEL"), s("*DEL")],
                "WholeDup1" => call_duplication(r),
                "WholeMultip1" => call_multiplication(r),
                "WholeMultip2" => vec![s("Indeterminate")],
                _ => vec![s("Indeterminate")],
            }
        }
        "UGT1A4" => {
            let (a1, a2) = (r.a1(), r.a2());
            match cnv {
                "Normal" | "AssumeNormal" => vec![s(a1), s(a2)],
                "NoncodingDel1" => call_linked(r, "*1", "*S1"),
                "NoncodingDel1Hom" => vec![s("*S1"), s("*S1")],
                "NoncodingDel2" => call_linked(r, "*1", "*S2"),
                "NoncodingDup1" => call_linked(r, "*1", "*S3"),
                _ => vec![s("Indeterminate")],
            }
        }
        "UGT2B15" => {
            let (a1, _a2) = (r.a1(), r.a2());
            match cnv {
                "Normal" | "AssumeNormal" => vec![s(a1), s(r.a2())],
                "WholeDel1" => vec![s(a1), s("*S4")],
                "PartialDel1" => vec![s(a1), s("*S1")],
                "PartialDel2" => vec![s(a1), s("*S2")],
                "PartialDel3" => vec![s(a1), s("*S3")],
                _ => vec![s("Indeterminate")],
            }
        }
        "UGT2B17" => match cnv {
            "Normal" | "AssumeNormal" => vec![s("*1"), s("*1")],
            "WholeDel1" => vec![s("*1"), s("*2")],
            "WholeDel1Hom" => vec![s("*2"), s("*2")],
            "PartialDel2" => vec![s("*1"), s("*S2")],
            "PartialDel3" => vec![s("*1"), s("*S3")],
            "WholeDel1+PartialDel1" => vec![s("*2"), s("*S1")],
            "WholeDel1+PartialDel2" => vec![s("*2"), s("*S2")],
            "WholeDel1+PartialDel3" => vec![s("*2"), s("*S3")],
            _ => vec![s("Indeterminate")],
        },
        // SimpleGenotyper: genes without SV.
        _ => vec![s(r.a1()), s(r.a2())],
    }
}

fn cyp2d6(assembly: &str, r: &GtRow) -> Vec<String> {
    let (a1, a2) = (r.a1(), r.a2());
    let s1 = core::sort_alleles(&[s(a1), s(a2)], "priority", Some("CYP2D6"), assembly)[0].clone();
    let cnv = r.cnv.as_deref().unwrap_or("");
    let tandem = |linked: &str, target: &str, alt: &str| -> Vec<String> {
        let h1 = r.h1_has(linked);
        let h2 = r.h2_has(linked);
        if h1 && h2 {
            vec![s(a1), s(target)]
        } else if h1 && !h2 {
            vec![s(a2), s(alt)]
        } else if !h1 && h2 {
            vec![s(a1), s(alt)]
        } else {
            vec![s("Indeterminate")]
        }
    };
    match cnv {
        "Normal" | "AssumeNormal" | "ParalogPartialDel1" => vec![s(a1), s(a2)],
        "WholeDel1" => vec![s1, s("*5")],
        "WholeDel1Hom" => vec![s("*5"), s("*5")],
        "WholeDup1" => call_duplication(r),
        "WholeMultip1" => call_multiplication(r),
        "WholeDel1+Tandem3" => vec![s("*5"), s("*13+*1")],
        "Tandem1A" => tandem("*4", "*68+*4", "*68+*4"),
        "Tandem1B" => {
            let (h1, h2) = (r.h1_has("*4"), r.h2_has("*4"));
            if h1 && h2 {
                vec![s("*68+*4"), s("*68+*4")]
            } else if h1 && !h2 {
                vec![s(a2), s("*68x2+*4")]
            } else if !h1 && h2 {
                vec![s(a1), s("*68x2+*4")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        "Tandem2A" => tandem("*10", "*36+*10", "*36+*10"),
        "Tandem2B" => {
            let (h1, h2) = (r.h1_has("*10"), r.h2_has("*10"));
            if h1 && h2 {
                vec![s("*36+*10"), s("*36+*10")]
            } else if h1 && !h2 {
                vec![s(a2), s("*36x2+*10")]
            } else if !h1 && h2 {
                vec![s(a1), s("*36x2+*10")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        "Tandem2C" => {
            let (h1, h2) = (r.h1_has("*10"), r.h2_has("*10"));
            if h1 && h2 {
                vec![s("*36+*10"), s("*36x2+*10")]
            } else if h1 && !h2 {
                vec![s(a2), s("*36x3+*10")]
            } else if !h1 && h2 {
                vec![s(a1), s("*36x3+*10")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        "Tandem3" => {
            let (h1, h2) = (r.h1_has("*1"), r.h2_has("*1"));
            if h1 && h2 {
                vec![s("*1"), s("*13+*1")]
            } else if h1 && !h2 {
                vec![s(a2), s("*13+*1")]
            } else if !h1 && h2 {
                vec![s(a1), s("*13+*1")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        "WholeDel1+Tandem1A" => {
            if a1.contains("*4") || a2.contains("*4") {
                vec![s("*5"), s("*68+*4")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        "WholeDup1+Tandem1A" => {
            let (h1, h2) = (r.h1_has("*4"), r.h2_has("*4"));
            if h1 && h2 {
                vec![s("*4x2"), s("*68+*4")]
            } else if h1 && !h2 {
                vec![format!("{a2}x2"), s("*68+*4")]
            } else if !h1 && h2 {
                vec![format!("{a1}x2"), s("*68+*4")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        _ => vec![s("Indeterminate")],
    }
}

fn slc22a2(r: &GtRow) -> Vec<String> {
    let (a1, a2) = (r.a1(), r.a2());
    let cnv = r.cnv.as_deref().unwrap_or("");
    match cnv {
        "Normal" | "AssumeNormal" => vec![s(a1), s(a2)],
        "NoncodingDel1" => call_linked(r, "*K432Q", "*S1"),
        "NoncodingDel1Hom" => vec![s("*S1"), s("*S1")],
        "PartialDel1" => call_linked(r, "*3", "*S2"),
        "NoncodingDel1+PartialDel1" => {
            if (r.h1_has("*3") || r.h2_has("*3")) && (r.h1_has("*K432Q") || r.h2_has("*K432Q")) {
                vec![s("*S1"), s("*S2")]
            } else {
                vec![s("Indeterminate")]
            }
        }
        _ => vec![s("Indeterminate")],
    }
}

/// `_call_duplication` — whole-gene duplication (3 copies).
fn call_duplication(r: &GtRow) -> Vec<String> {
    let (a1, a2) = (r.a1(), r.a2());
    let h1 = if r.vd(a1).truthy() {
        r.vd(a1).fractions().iter().all(|&x| x > 0.5)
    } else {
        true
    };
    let h2 = if r.vd(a2).truthy() {
        r.vd(a2).fractions().iter().all(|&x| x > 0.5)
    } else {
        true
    };
    if h1 && h2 {
        if a1 == a2 {
            vec![s(a2), format!("{a1}x2")]
        } else if r.vd(a1).truthy() && !r.vd(a2).truthy() {
            vec![s(a2), format!("{a1}x2")]
        } else if !r.vd(a1).truthy() && r.vd(a2).truthy() {
            vec![s(a1), format!("{a2}x2")]
        } else if r.h2_has(a1) {
            vec![s(a1), format!("{a2}x2")]
        } else if r.h1_has(a2) {
            vec![s(a2), format!("{a1}x2")]
        } else {
            vec![s("Indeterminate")]
        }
    } else if h1 && !h2 {
        vec![s(a2), format!("{a1}x2")]
    } else if !h1 && h2 {
        vec![s(a1), format!("{a2}x2")]
    } else {
        vec![s("Indeterminate")]
    }
}

/// `_call_multiplication` — whole-gene multiplication (4 copies).
fn call_multiplication(r: &GtRow) -> Vec<String> {
    let (a1, a2) = (r.a1(), r.a2());
    let h1 = if r.vd(a1).truthy() {
        r.vd(a1).fractions().iter().all(|&x| x > 0.6)
    } else {
        true
    };
    let h2 = if r.vd(a2).truthy() {
        r.vd(a2).fractions().iter().all(|&x| x > 0.6)
    } else {
        true
    };
    if h1 && h2 {
        if a1 == a2 {
            vec![s(a2), format!("{a1}x3")]
        } else if r.vd(a1).truthy() && !r.vd(a2).truthy() {
            vec![s(a2), format!("{a1}x3")]
        } else if !r.vd(a1).truthy() && r.vd(a2).truthy() {
            vec![s(a1), format!("{a2}x3")]
        } else {
            vec![s("Indeterminate")]
        }
    } else if h1 && !h2 {
        vec![s(a2), format!("{a1}x3")]
    } else if !h1 && h2 {
        vec![s(a1), format!("{a2}x3")]
    } else {
        vec![s("Indeterminate")]
    }
}

/// `_call_linked_allele`.
fn call_linked(r: &GtRow, linked: &str, target: &str) -> Vec<String> {
    let (a1, a2) = (r.a1(), r.a2());
    let h1 = r.h1_has(linked);
    let h2 = r.h2_has(linked);
    if h1 && h2 {
        vec![s(a1), s(target)]
    } else if h1 && !h2 {
        vec![s(a2), s(target)]
    } else if !h1 && h2 {
        vec![s(a1), s(target)]
    } else {
        vec![s("Indeterminate")]
    }
}
