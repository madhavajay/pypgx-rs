//! Port of `pypgx.api.core` — reference tables plus the allele / phenotype /
//! score logic. Tables are embedded from `data/*.csv` (mirroring PyPGx's
//! `package_data`) and parsed with the same NaN semantics as `pandas.read_csv`.

use std::collections::{HashMap, HashSet};

use crate::fuc::{self, VcfFrame};
use crate::sdk::PgxError;
use crate::table::{Cell, Frame};

/// Allele-function priority order (`core.FUNCTION_ORDER`), used by
/// `sort_alleles(by='priority')`.
pub const FUNCTION_ORDER: &[&str] = &[
    "No Function",
    "Severely Decreased Function",
    "Decreased Function",
    "Possible Decreased Function",
    "Increased Function",
    "Possible Increased Function",
    "Class I (Deficient with CNSHA)",
    "Class II (Deficient)",
    "Class III (Deficient)",
    "Unfavorable Response",
    "Malignant Hyperthermia Associated",
    "Uncertain Function",
    "Unknown Function",
    "Normal Function",
    "Favorable Response",
    "Class IV (Normal)",
];

// ---- embedded reference tables -------------------------------------------

const ALLELE_CSV: &str = include_str!("../data/allele-table.csv");
const CNV_CSV: &str = include_str!("../data/cnv-table.csv");
const CPIC_CSV: &str = include_str!("../data/cpic-table.csv");
const DIPLOTYPE_CSV: &str = include_str!("../data/diplotype-table.csv");
const EQUATION_CSV: &str = include_str!("../data/equation-table.csv");
const GENE_CSV: &str = include_str!("../data/gene-table.csv");
const PHENOTYPE_CSV: &str = include_str!("../data/phenotype-table.csv");
const RECOMMENDATION_CSV: &str = include_str!("../data/recommendation-table.csv");
const VARIANT_CSV: &str = include_str!("../data/variant-table.csv");

pub fn load_allele_table() -> Frame {
    Frame::from_csv(ALLELE_CSV, true)
}
pub fn load_cnv_table() -> Frame {
    Frame::from_csv(CNV_CSV, true)
}
pub fn load_cpic_table() -> Frame {
    Frame::from_csv(CPIC_CSV, true)
}
pub fn load_diplotype_table() -> Frame {
    Frame::from_csv(DIPLOTYPE_CSV, true)
}
pub fn load_equation_table() -> Frame {
    Frame::from_csv(EQUATION_CSV, true)
}
pub fn load_gene_table() -> Frame {
    Frame::from_csv(GENE_CSV, true)
}
pub fn load_phenotype_table() -> Frame {
    Frame::from_csv(PHENOTYPE_CSV, true)
}
pub fn load_recommendation_table() -> Frame {
    // PyPGx loads this one with na_filter=False.
    Frame::from_csv(RECOMMENDATION_CSV, false)
}
pub fn load_variant_table() -> Frame {
    // PyPGx casts Chromosome to str; our cells are already strings.
    Frame::from_csv(VARIANT_CSV, true)
}

// ---- gene listing --------------------------------------------------------

/// `list_genes(mode)` — genes filtered by `target`/`control`/`all`.
pub fn list_genes(mode: &str) -> Vec<String> {
    let df = load_gene_table();
    let gene_c = df.col("Gene");
    type RowPred = Box<dyn Fn(&Vec<Cell>) -> bool>;
    let keep: RowPred = match mode {
        "target" => {
            let c = df.col("Target");
            Box::new(move |r| r[c].is_true())
        }
        "control" => {
            let c = df.col("Control");
            Box::new(move |r| r[c].is_true())
        }
        _ => Box::new(|_| true),
    };
    df.rows
        .iter()
        .filter(|r| keep(r))
        .filter_map(|r| r[gene_c].as_str().map(|s| s.to_string()))
        .collect()
}

/// `is_target_gene(gene)`.
pub fn is_target_gene(gene: &str) -> bool {
    list_genes("target").iter().any(|g| g == gene)
}

/// `has_phenotype(gene)` — whether the gene has a PhenotypeMethod.
pub fn has_phenotype(gene: &str) -> Result<bool, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    let df = load_gene_table();
    let gene_c = df.col("Gene");
    let pm_c = df.col("PhenotypeMethod");
    let genes: Vec<&str> = df
        .rows
        .iter()
        .filter(|r| !r[pm_c].is_null())
        .filter_map(|r| r[gene_c].as_str())
        .collect();
    Ok(genes.contains(&gene))
}

// ---- simple accessors ----------------------------------------------------

fn gene_field(gene: &str, field: &str) -> Option<String> {
    let df = load_gene_table();
    let rows = df.filter_eq("Gene", gene);
    let row = rows.first().expect("gene not in gene table");
    row[df.col(field)].as_str().map(|s| s.to_string())
}

/// `get_ref_allele(gene)`.
pub fn get_ref_allele(gene: &str) -> String {
    gene_field(gene, "RefAllele").expect("RefAllele")
}

/// `get_default_allele(gene, assembly)`.
pub fn get_default_allele(gene: &str, assembly: &str) -> String {
    gene_field(gene, &format!("{assembly}Default")).expect("Default allele")
}

/// `get_strand(gene)`.
pub fn get_strand(gene: &str) -> Result<String, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    Ok(gene_field(gene, "Strand").expect("Strand"))
}

/// `get_paralog(gene)` — empty string when none.
pub fn get_paralog(gene: &str) -> String {
    gene_field(gene, "Paralog").unwrap_or_default()
}

/// `get_function(gene, allele)` — `None` mirrors a NaN function.
pub fn get_function(gene: &str, allele: &str) -> Result<Option<String>, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    let df = load_allele_table();
    let gene_c = df.col("Gene");
    let sa_c = df.col("StarAllele");
    let fn_c = df.col("Function");
    let row = df
        .rows
        .iter()
        .find(|r| r[gene_c].as_str() == Some(gene) && r[sa_c].as_str() == Some(allele));
    match row {
        None => Err(PgxError::AlleleNotFound {
            gene: gene.to_string(),
            allele: allele.to_string(),
        }),
        Some(r) => Ok(r[fn_c].as_str().map(|s| s.to_string())),
    }
}

/// `get_variant_impact(variant)` — empty string for a NaN impact.
pub fn get_variant_impact(variant: &str) -> Result<String, PgxError> {
    let df = load_variant_table();
    let g37 = df.col("GRCh37Name");
    let g38 = df.col("GRCh38Name");
    let imp = df.col("Impact");
    let row = df
        .rows
        .iter()
        .find(|r| r[g37].as_str() == Some(variant) || r[g38].as_str() == Some(variant));
    match row {
        None => Err(PgxError::VariantNotFound(variant.to_string())),
        Some(r) => Ok(r[imp].as_str().unwrap_or("").to_string()),
    }
}

/// `get_variant_synonyms(gene, assembly)`.
pub fn get_variant_synonyms(gene: &str, assembly: &str) -> HashMap<String, String> {
    let df = load_variant_table();
    let gene_c = df.col("Gene");
    let syn_c = df.col(&format!("{assembly}Synonym"));
    let name_c = df.col(&format!("{assembly}Name"));
    let mut synonyms = HashMap::new();
    for r in df.rows.iter().filter(|r| r[gene_c].as_str() == Some(gene)) {
        if let Some(syn) = r[syn_c].as_str() {
            for v in syn.split(',') {
                if let Some(name) = r[name_c].as_str() {
                    synonyms.insert(v.to_string(), name.to_string());
                }
            }
        }
    }
    synonyms
}

// ---- allele / variant listing -------------------------------------------

/// `list_alleles(gene, variants, assembly)`.
pub fn list_alleles(gene: &str, variants: Option<&[String]>, assembly: &str) -> Vec<String> {
    if !is_target_gene(gene) {
        panic!("NotTargetGeneError: {gene}");
    }
    let df = load_allele_table();
    let gene_c = df.col("Gene");
    let sa_c = df.col("StarAllele");
    let core_c = df.col(&format!("{assembly}Core"));
    let tag_c = df.col(&format!("{assembly}Tag"));
    df.rows
        .iter()
        .filter(|r| r[gene_c].as_str() == Some(gene))
        .filter(|r| match variants {
            None => true,
            Some(vs) => {
                let mut l: Vec<&str> = Vec::new();
                if let Some(c) = r[core_c].as_str() {
                    l.extend(c.split(','));
                }
                if let Some(t) = r[tag_c].as_str() {
                    l.extend(t.split(','));
                }
                vs.iter().all(|x| l.contains(&x.as_str()))
            }
        })
        .filter_map(|r| r[sa_c].as_str().map(|s| s.to_string()))
        .collect()
}

/// `list_variants(gene, alleles, mode, assembly)` — coordinate-sorted, unique.
pub fn list_variants(
    gene: &str,
    alleles: Option<&[String]>,
    mode: &str,
    assembly: &str,
) -> Vec<String> {
    if !is_target_gene(gene) {
        panic!("NotTargetGeneError: {gene}");
    }
    let df = load_allele_table();
    let gene_c = df.col("Gene");
    let sa_c = df.col("StarAllele");
    let core_c = df.col(&format!("{assembly}Core"));
    let tag_c = df.col(&format!("{assembly}Tag"));

    let owned;
    let alleles: &[String] = match alleles {
        Some(a) => a,
        None => {
            owned = list_alleles(gene, None, assembly);
            &owned
        }
    };

    let mut core_variants: Vec<String> = Vec::new();
    let mut tag_variants: Vec<String> = Vec::new();
    for allele in alleles {
        let row = df
            .rows
            .iter()
            .find(|r| r[gene_c].as_str() == Some(gene) && r[sa_c].as_str() == Some(allele));
        let row = row.unwrap_or_else(|| panic!("AlleleNotFoundError: {gene}/{allele}"));
        if let Some(c) = row[core_c].as_str() {
            core_variants.extend(c.split(',').map(|s| s.to_string()));
        }
        if let Some(t) = row[tag_c].as_str() {
            tag_variants.extend(t.split(',').map(|s| s.to_string()));
        }
    }

    let results: Vec<String> = match mode {
        "all" => core_variants.into_iter().chain(tag_variants).collect(),
        "core" => core_variants,
        "tag" => tag_variants,
        _ => panic!("Incorrect mode: {mode}"),
    };
    let set: HashSet<String> = results.into_iter().collect();
    fuc::sort_variants(set)
}

// ---- definition table ----------------------------------------------------

/// `build_definition_table(gene, assembly)` — a VcfFrame of SNV/indel-defined
/// star alleles (sorted by coordinate).
pub fn build_definition_table(gene: &str, assembly: &str) -> VcfFrame {
    if !is_target_gene(gene) {
        panic!("NotTargetGeneError: {gene}");
    }
    let at = load_allele_table();
    let gene_c = at.col("Gene");
    let sa_c = at.col("StarAllele");
    let sv_c = at.col("SV");
    let core_c = at.col(&format!("{assembly}Core"));
    let gene_rows: Vec<&Vec<Cell>> = at
        .rows
        .iter()
        .filter(|r| r[gene_c].as_str() == Some(gene))
        .collect();

    // Ordered unique list of core variants across this gene's alleles.
    let mut variants: Vec<String> = Vec::new();
    for r in &gene_rows {
        if let Some(core) = r[core_c].as_str() {
            for v in core.split(',') {
                if !variants.iter().any(|x| x == v) {
                    variants.push(v.to_string());
                }
            }
        }
    }

    // Allele genotype columns (skip SV alleles and alleles with no core).
    let mut allele_names: Vec<String> = Vec::new();
    let mut allele_cols: Vec<Vec<String>> = Vec::new();
    for r in &gene_rows {
        let is_sv = r[sv_c].is_true();
        let core = r[core_c].as_str();
        if is_sv || core.is_none() {
            continue;
        }
        let core_set: Vec<&str> = core.unwrap().split(',').collect();
        let col: Vec<String> = variants
            .iter()
            .map(|x| {
                if core_set.contains(&x.as_str()) {
                    "1"
                } else {
                    "0"
                }
                .to_string()
            })
            .collect();
        allele_names.push(r[sa_c].as_str().unwrap().to_string());
        allele_cols.push(col);
    }

    // Per-variant VCF fields from the variant table.
    let vt = load_variant_table();
    let vt_gene_c = vt.col("Gene");
    let vt_name_c = vt.col(&format!("{assembly}Name"));
    let vt_chrom_c = vt.col("Chromosome");
    let vt_rsid_c = vt.col("rsID");
    let vt_impact_c = vt.col("Impact");
    let vt_gene_rows: Vec<&Vec<Cell>> = vt
        .rows
        .iter()
        .filter(|r| r[vt_gene_c].as_str() == Some(gene))
        .collect();

    let mut columns: Vec<String> = fuc::HEADERS.iter().map(|s| s.to_string()).collect();
    columns.extend(allele_names.iter().cloned());

    let mut rows: Vec<Vec<String>> = Vec::new();
    for (i, variant) in variants.iter().enumerate() {
        let parts: Vec<&str> = variant.split('-').collect();
        let pos = parts[1];
        let rf = parts[2];
        let alt = parts[3];
        let s = vt_gene_rows
            .iter()
            .find(|r| r[vt_name_c].as_str() == Some(variant.as_str()))
            .unwrap_or_else(|| panic!("variant {variant} not in variant table"));
        let chrom = s[vt_chrom_c].as_str().unwrap_or("nan").to_string();
        let id = s[vt_rsid_c].as_str().unwrap_or("nan").to_string();
        let impact = s[vt_impact_c].as_str().unwrap_or("nan").to_string();
        let mut row = vec![
            chrom,
            pos.to_string(),
            id,
            rf.to_string(),
            alt.to_string(),
            ".".to_string(),
            ".".to_string(),
            format!("VI={impact}"),
            "GT".to_string(),
        ];
        for col in &allele_cols {
            row.push(col[i].clone());
        }
        rows.push(row);
    }

    let meta = vec![
        "##fileformat=VCFv4.1".to_string(),
        "##INFO=<ID=VI,Number=1,Type=String,Description=\"Variant impact\">".to_string(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
    ];
    VcfFrame::new(meta, columns, rows).sort()
}

// ---- collapse / sort -----------------------------------------------------

/// `collapse_alleles(gene, alleles, assembly)` — drop alleles whose core
/// variants are a subset of another candidate's.
pub fn collapse_alleles(gene: &str, alleles: &[String], assembly: &str) -> Vec<String> {
    let mut redundant = vec![false; alleles.len()];
    for (i, a) in alleles.iter().enumerate() {
        for b in alleles.iter() {
            if a == b {
                continue;
            }
            let v1: HashSet<String> =
                list_variants(gene, Some(std::slice::from_ref(a)), "core", assembly)
                    .into_iter()
                    .collect();
            let v2: HashSet<String> =
                list_variants(gene, Some(std::slice::from_ref(b)), "core", assembly)
                    .into_iter()
                    .collect();
            if v1.is_subset(&v2) {
                redundant[i] = true;
                break;
            }
        }
    }
    alleles
        .iter()
        .enumerate()
        .filter(|(i, _)| !redundant[*i])
        .map(|(_, a)| a.clone())
        .collect()
}

/// `sort_alleles(alleles, by, gene, assembly)`.
pub fn sort_alleles(
    alleles: &[String],
    by: &str,
    gene: Option<&str>,
    assembly: &str,
) -> Vec<String> {
    let mut out = alleles.to_vec();
    match by {
        "priority" => {
            let gene = gene.expect("Gene is required when sorting by priority");
            out.sort_by_cached_key(|allele| priority_key(allele, gene, assembly));
        }
        "name" => {
            out.sort_by_cached_key(|allele| name_key(allele));
        }
        _ => panic!("unknown sort mode: {by}"),
    }
    out
}

fn priority_key(allele: &str, gene: &str, assembly: &str) -> (i64, i64, i64, bool) {
    let function = get_function(gene, allele)
        .expect("get_function")
        .expect("function is NaN (cannot index FUNCTION_ORDER)");
    let a = FUNCTION_ORDER
        .iter()
        .position(|f| *f == function)
        .expect("function in FUNCTION_ORDER") as i64;
    let core_variants = list_variants(
        gene,
        Some(std::slice::from_ref(&allele.to_string())),
        "core",
        assembly,
    );
    let b = -(core_variants.len() as i64);
    let impacts: Vec<String> = core_variants
        .iter()
        .map(|x| get_variant_impact(x).expect("impact"))
        .filter(|x| !x.is_empty())
        .collect();
    let c = -(impacts.len() as i64);
    let d = allele == get_ref_allele(gene);
    (a, b, c, d)
}

fn name_key(allele: &str) -> (i64, i64, usize, String) {
    let mut n: i64 = 99999;
    let mut cn: i64 = 1;
    if allele == "Reference" {
        n = 0;
    } else if allele == "Indeterminate" {
        n += 1;
    } else if allele.contains("c.") {
        // DPYD-style names: digits before the first '>'.
        let head = allele.split('>').next().unwrap_or(allele);
        let digits: String = head.chars().filter(|c| c.is_ascii_digit()).collect();
        if let Ok(v) = digits.parse::<i64>() {
            n = v;
        }
    } else if !allele.contains('*') {
        // leave n at default
    } else {
        let first = allele.split('+').next().unwrap_or(allele);
        let base = first.split('x').next().unwrap_or(first).replace('*', "");
        let first_char = base.chars().next();
        if first_char.map(|c| c.is_ascii_digit()).unwrap_or(false) {
            let digits: String = base.chars().filter(|c| c.is_ascii_digit()).collect();
            if let Ok(v) = digits.parse::<i64>() {
                n = v;
            }
        }
        if first.contains('x') {
            if let Some(times) = first.split('x').nth(1) {
                if let Ok(v) = times.parse::<i64>() {
                    cn = v;
                }
            }
        }
    }
    (n, cn, allele.chars().count(), allele.to_string())
}

// ---- phenotype / score ---------------------------------------------------

/// `has_score(gene)`.
pub fn has_score(gene: &str) -> Result<bool, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    let df = load_gene_table();
    Ok(df
        .filter_eq("Gene", gene)
        .first()
        .map(|r| r[df.col("PhenotypeMethod")].as_str() == Some("Score"))
        .unwrap_or(false))
}

/// `has_sv(gene)` (gene-level form).
pub fn has_sv(gene: &str) -> Result<bool, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    let df = load_gene_table();
    Ok(df
        .filter_eq("Gene", gene)
        .first()
        .map(|r| r[df.col("SV")].is_true())
        .unwrap_or(false))
}

/// `get_score(gene, allele)` — `None` mirrors NaN.
pub fn get_score(gene: &str, allele: &str) -> Result<Option<f64>, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    if !has_score(gene)? {
        return Ok(None);
    }
    let df = load_allele_table();
    let gene_c = df.col("Gene");
    let sa_c = df.col("StarAllele");
    let as_c = df.col("ActivityScore");
    match df
        .rows
        .iter()
        .find(|r| r[gene_c].as_str() == Some(gene) && r[sa_c].as_str() == Some(allele))
    {
        None => Err(PgxError::AlleleNotFound {
            gene: gene.to_string(),
            allele: allele.to_string(),
        }),
        Some(r) => Ok(r[as_c].as_str().and_then(|s| s.parse::<f64>().ok())),
    }
}

/// `is_legit_allele(gene, allele)`.
pub fn is_legit_allele(gene: &str, allele: &str) -> bool {
    list_alleles(gene, None, "GRCh37")
        .iter()
        .any(|a| a == allele)
}

/// `list_functions(gene)` — unique Function values (NaN kept as `None`).
pub fn list_functions(gene: Option<&str>) -> Vec<Option<String>> {
    let df = load_allele_table();
    if let Some(g) = gene {
        if !is_target_gene(g) {
            panic!("NotTargetGeneError: {g}");
        }
    }
    let gene_c = df.col("Gene");
    let fn_c = df.col("Function");
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for r in &df.rows {
        if let Some(g) = gene {
            if r[gene_c].as_str() != Some(g) {
                continue;
            }
        }
        let key = r[fn_c].as_str().map(|s| s.to_string());
        if seen.insert(key.clone()) {
            out.push(key);
        }
    }
    out
}

/// `get_region(gene, assembly)`.
pub fn get_region(gene: &str, assembly: &str) -> Result<String, PgxError> {
    if !list_genes("all").iter().any(|g| g == gene) {
        return Err(PgxError::GeneNotFound(gene.to_string()));
    }
    Ok(gene_field(gene, &format!("{assembly}Region")).expect("Region"))
}

fn exon_positions(gene: &str, field: &str) -> Result<Vec<i64>, PgxError> {
    if !list_genes("all").iter().any(|g| g == gene) {
        return Err(PgxError::GeneNotFound(gene.to_string()));
    }
    let s = gene_field(gene, field).expect("exon field");
    Ok(s.trim_matches(',')
        .split(',')
        .map(|x| x.parse::<i64>().expect("exon position"))
        .collect())
}

/// `get_exon_starts(gene, assembly)`.
pub fn get_exon_starts(gene: &str, assembly: &str) -> Result<Vec<i64>, PgxError> {
    exon_positions(gene, &format!("{assembly}ExonStarts"))
}

/// `get_exon_ends(gene, assembly)`.
pub fn get_exon_ends(gene: &str, assembly: &str) -> Result<Vec<i64>, PgxError> {
    exon_positions(gene, &format!("{assembly}ExonEnds"))
}

/// `get_priority(gene, phenotype)`.
pub fn get_priority(gene: &str, phenotype: &str) -> Result<String, PgxError> {
    if !is_target_gene(gene) {
        return Err(PgxError::NotTargetGene(gene.to_string()));
    }
    if !list_phenotypes(None).iter().any(|p| p == phenotype) {
        return Err(PgxError::PhenotypeNotFound(phenotype.to_string()));
    }
    let df = load_phenotype_table();
    let gene_c = df.col("Gene");
    let ph_c = df.col("Phenotype");
    let pr_c = df.col("Priority");
    let row = df
        .rows
        .iter()
        .find(|r| r[gene_c].as_str() == Some(gene) && r[ph_c].as_str() == Some(phenotype))
        .expect("phenotype row");
    Ok(row[pr_c].as_str().unwrap_or("").to_string())
}

/// Score as `f64` (NaN mirrors a pandas NaN); panics on AlleleNotFound, as
/// PyPGx raises.
fn get_score_f64(gene: &str, allele: &str) -> f64 {
    get_score(gene, allele)
        .expect("get_score")
        .unwrap_or(f64::NAN)
}

/// `predict_score(gene, allele)` — NaN when the gene has no activity score or
/// the allele has uncertain function. Handles SV (`x`, `+`).
pub fn predict_score(gene: &str, allele: &str) -> f64 {
    if !is_target_gene(gene) {
        panic!("NotTargetGeneError: {gene}");
    }
    if !has_score(gene).unwrap() {
        return f64::NAN;
    }
    if has_sv(gene).unwrap() {
        allele
            .split('+')
            .map(|x| {
                if x.contains('x') {
                    let mut it = x.split('x');
                    let base = it.next().unwrap();
                    let times: f64 = it.next().unwrap().parse().expect("cn");
                    get_score_f64(gene, base) * times
                } else {
                    get_score_f64(gene, x)
                }
            })
            .sum()
    } else {
        get_score_f64(gene, allele)
    }
}

/// Evaluate a chained comparison from the equation table (e.g.
/// `0 <= score < 0.25`) with `score` substituted, mirroring Python's `eval`.
fn eval_equation(equation: &str, score: f64) -> bool {
    let tokens: Vec<&str> = equation.split_whitespace().collect();
    let val = |tok: &str| -> f64 {
        if tok == "score" {
            score
        } else {
            tok.parse().expect("equation literal")
        }
    };
    // Chained comparison: operands at even indices, operators at odd indices.
    let mut i = 0;
    while i + 2 < tokens.len() {
        let l = val(tokens[i]);
        let r = val(tokens[i + 2]);
        let ok = match tokens[i + 1] {
            "<=" => l <= r,
            "<" => l < r,
            ">=" => l >= r,
            ">" => l > r,
            "==" => l == r,
            "!=" => l != r,
            op => panic!("unknown operator: {op}"),
        };
        if !ok {
            return false;
        }
        i += 2;
    }
    true
}

/// `predict_phenotype(gene, a, b)`.
pub fn predict_phenotype(gene: &str, a: &str, b: &str) -> String {
    if !is_target_gene(gene) {
        panic!("NotTargetGeneError: {gene}");
    }
    let method = gene_field(gene, "PhenotypeMethod");
    match method.as_deref() {
        Some("Score") => {
            let score = predict_score(gene, a) + predict_score(gene, b);
            if score.is_nan() {
                return "Indeterminate".to_string();
            }
            let df = load_equation_table();
            let gene_c = df.col("Gene");
            let eq_c = df.col("Equation");
            let ph_c = df.col("Phenotype");
            for r in df.rows.iter().filter(|r| r[gene_c].as_str() == Some(gene)) {
                if eval_equation(r[eq_c].as_str().unwrap(), score) {
                    return r[ph_c].as_str().unwrap().to_string();
                }
            }
            "Indeterminate".to_string()
        }
        Some("Diplotype") => {
            let df = load_diplotype_table();
            let gene_c = df.col("Gene");
            let dp_c = df.col("Diplotype");
            let ph_c = df.col("Phenotype");
            let candidates = [format!("{a}/{b}"), format!("{b}/{a}")];
            df.rows
                .iter()
                .filter(|r| r[gene_c].as_str() == Some(gene))
                .find(|r| {
                    r[dp_c]
                        .as_str()
                        .map(|d| candidates.iter().any(|c| c == d))
                        .unwrap_or(false)
                })
                .map(|r| r[ph_c].as_str().unwrap().to_string())
                .unwrap_or_else(|| "Indeterminate".to_string())
        }
        _ => "Indeterminate".to_string(),
    }
}

/// `list_phenotypes(gene)` — sorted unique phenotypes.
pub fn list_phenotypes(gene: Option<&str>) -> Vec<String> {
    let df = load_phenotype_table();
    if let Some(g) = gene {
        if !is_target_gene(g) {
            panic!("NotTargetGeneError: {g}");
        }
    }
    let gene_c = df.col("Gene");
    let ph_c = df.col("Phenotype");
    let mut set: Vec<String> = Vec::new();
    for r in &df.rows {
        if let Some(g) = gene {
            if r[gene_c].as_str() != Some(g) {
                continue;
            }
        }
        if let Some(p) = r[ph_c].as_str() {
            if !set.contains(&p.to_string()) {
                set.push(p.to_string());
            }
        }
    }
    set.sort();
    set
}

/// `get_recommendation(drug, gene1, phenotype1, gene2, phenotype2)` — drug
/// recommendation for a phenotype (or phenotype pair).
pub fn get_recommendation(
    drug: &str,
    gene1: &str,
    phenotype1: &str,
    gene2: Option<&str>,
    phenotype2: Option<&str>,
) -> Result<String, PgxError> {
    let all = list_genes("all");
    if !all.iter().any(|g| g == gene1) {
        return Err(PgxError::GeneNotFound(gene1.to_string()));
    }
    if let Some(g2) = gene2 {
        if !all.iter().any(|g| g == g2) {
            return Err(PgxError::GeneNotFound(g2.to_string()));
        }
    }
    if !list_phenotypes(Some(gene1)).iter().any(|p| p == phenotype1) {
        return Err(PgxError::PhenotypeNotFound(format!(
            "{phenotype1} in {gene1}"
        )));
    }
    if let (Some(g2), Some(p2)) = (gene2, phenotype2) {
        if !list_phenotypes(Some(g2)).iter().any(|p| p == p2) {
            return Err(PgxError::PhenotypeNotFound(format!("{p2} in {g2}")));
        }
    }

    let df = load_recommendation_table();
    let drug_c = df.col("Drug");
    let g1_c = df.col("Gene1");
    let p1_c = df.col("Phenotype1");
    let g2_c = df.col("Gene2");
    let p2_c = df.col("Phenotype2");
    let rec_c = df.col("Recommendation");

    if !df.rows.iter().any(|r| r[drug_c].as_str() == Some(drug)) {
        panic!("Drug not found: {drug}");
    }
    let rows: Vec<&Vec<Cell>> = df
        .rows
        .iter()
        .filter(|r| r[drug_c].as_str() == Some(drug))
        .collect();

    let mut target_genes: HashSet<&str> = HashSet::new();
    for r in &rows {
        if let Some(g) = r[g1_c].as_str() {
            target_genes.insert(g);
        }
        if let Some(g) = r[g2_c].as_str() {
            target_genes.insert(g);
        }
    }
    if !target_genes.contains(gene1) {
        panic!("{gene1} does not have any recommendations for {drug}");
    }
    if let Some(g2) = gene2 {
        if !target_genes.contains(g2) {
            panic!("{g2} does not have any recommendations for {drug}");
        }
    }

    let uniq = |col: usize| -> Vec<String> {
        let mut seen = Vec::new();
        for r in &rows {
            if let Some(v) = r[col].as_str() {
                if !seen.iter().any(|x| x == v) {
                    seen.push(v.to_string());
                }
            }
        }
        seen
    };
    let rec_of = |row: Option<&&Vec<Cell>>| -> String {
        row.expect("recommendation row")[rec_c]
            .as_str()
            .unwrap_or("")
            .to_string()
    };

    // Single-gene drug: every Gene2 is the literal string 'None'.
    if uniq(g2_c) == ["None"] {
        let row = rows
            .iter()
            .find(|r| r[g1_c].as_str() == Some(gene1) && r[p1_c].as_str() == Some(phenotype1));
        return Ok(rec_of(row));
    }

    let gene1_in_g1 = uniq(g1_c).iter().any(|g| g == gene1);

    if gene2.is_none() {
        // Multi-gene drug queried with a single gene (PyPGx warns here).
        let row = if gene1_in_g1 {
            rows.iter()
                .find(|r| r[g1_c].as_str() == Some(gene1) && r[p1_c].as_str() == Some(phenotype1))
        } else {
            rows.iter()
                .find(|r| r[g2_c].as_str() == Some(gene1) && r[p2_c].as_str() == Some(phenotype1))
        };
        return Ok(rec_of(row));
    }

    let gene2 = gene2.unwrap();
    let phenotype2 = phenotype2.unwrap();
    let row = if gene1_in_g1 {
        rows.iter().find(|r| {
            r[g1_c].as_str() == Some(gene1)
                && r[p1_c].as_str() == Some(phenotype1)
                && r[g2_c].as_str() == Some(gene2)
                && r[p2_c].as_str() == Some(phenotype2)
        })
    } else {
        rows.iter().find(|r| {
            r[g2_c].as_str() == Some(gene1)
                && r[p2_c].as_str() == Some(phenotype1)
                && r[g1_c].as_str() == Some(gene2)
                && r[p1_c].as_str() == Some(phenotype2)
        })
    };
    Ok(rec_of(row))
}
