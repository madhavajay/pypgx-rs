//! Port of the `pypgx.api.utils` functions. Currently implements
//! `predict_alleles`; external-tool wrappers (Beagle phasing, depth, CNV) are
//! tracked in TODO.md and will shell out to the same programs as PyPGx.

use std::collections::{HashMap, HashSet};

use crate::bed::BedFrame;
use crate::core;
use crate::fuc::{python_float_str, sort_variants, VcfFrame};
use crate::sdk::{Archive, ArchiveData, PgxError, SampleTable};

/// `predict_alleles(consolidated_variants)` — predict candidate star alleles
/// from observed SNVs/indels. Input is a `VcfFrame[Consolidated]` archive;
/// output is a `SampleTable[Alleles]` archive.
pub fn predict_alleles(consolidated_variants: &Archive) -> Result<Archive, PgxError> {
    consolidated_variants.check_type(&["VcfFrame[Consolidated]"])?;

    let gene = consolidated_variants
        .get("Gene")
        .expect("Gene metadata")
        .to_string();
    let assembly = consolidated_variants
        .get("Assembly")
        .expect("Assembly metadata")
        .to_string();

    let definition_table = core::build_definition_table(&gene, &assembly);
    let ref_allele = core::get_ref_allele(&gene);
    let default_allele = core::get_default_allele(&gene, &assembly);
    let defining_variants: HashSet<String> = core::list_variants(&gene, None, "all", &assembly)
        .into_iter()
        .collect();
    let variant_synonyms = core::get_variant_synonyms(&gene, &assembly);

    let vcf = consolidated_variants.as_vcf();

    // reformatted_variants: synonym name -> observed variant.
    let mut reformatted_variants: HashMap<String, String> = HashMap::new();
    for x in vcf.to_variants() {
        if let Some(y) = variant_synonyms.get(&x) {
            reformatted_variants.insert(y.clone(), x.clone());
        }
    }

    // star_alleles in definition-table sample order; values are the (ordered,
    // de-duplicated) defining variants of each allele.
    let mut star_alleles: Vec<(String, Vec<String>)> = Vec::new();
    let def_cols = &definition_table.columns;
    for allele in definition_table.samples() {
        let ac = def_cols.iter().position(|c| *c == allele).unwrap();
        let chrom_c = def_cols.iter().position(|c| c == "CHROM").unwrap();
        let pos_c = def_cols.iter().position(|c| c == "POS").unwrap();
        let ref_c = def_cols.iter().position(|c| c == "REF").unwrap();
        let alt_c = def_cols.iter().position(|c| c == "ALT").unwrap();
        let mut vars: Vec<String> = Vec::new();
        for r in &definition_table.rows {
            if r[ac] == "1" {
                let v = format!("{}-{}-{}-{}", r[chrom_c], r[pos_c], r[ref_c], r[alt_c]);
                if !vars.contains(&v) {
                    vars.push(v);
                }
            }
        }
        star_alleles.push((allele, sort_variants(vars)));
    }
    let star_sets: Vec<HashSet<String>> = star_alleles
        .iter()
        .map(|(_, v)| v.iter().cloned().collect())
        .collect();

    // one_haplotype: candidate alleles whose defining variants ⊆ observed.
    let one_haplotype = |observed: &HashSet<String>| -> Vec<String> {
        let mut candidates: Vec<String> = Vec::new();
        for (idx, (allele, _)) in star_alleles.iter().enumerate() {
            if star_sets[idx].is_subset(observed) {
                candidates.push(allele.clone());
            }
        }
        let mut candidates = core::collapse_alleles(&gene, &candidates, &assembly);
        if ref_allele != default_allele
            && !candidates.contains(&ref_allele)
            && !candidates.contains(&default_allele)
        {
            candidates.push(default_allele.clone());
        }
        if candidates.is_empty() {
            candidates.push(default_allele.clone());
        }
        core::sort_alleles(&candidates, "priority", Some(&gene), &assembly)
    };

    // one_row: the observed defining variant for sample/haplotype, or "".
    let chrom_c = vcf.columns.iter().position(|c| c == "CHROM").unwrap();
    let pos_c = vcf.columns.iter().position(|c| c == "POS").unwrap();
    let ref_c = vcf.columns.iter().position(|c| c == "REF").unwrap();
    let alt_c = vcf.columns.iter().position(|c| c == "ALT").unwrap();
    let one_row = |row: &[String], sample_c: usize, i: usize| -> String {
        let cell = &row[sample_c];
        let gt = cell.split(':').next().unwrap_or("");
        if gt.contains('.') {
            return String::new();
        }
        let j: usize = match gt.split('|').nth(i).and_then(|x| x.parse().ok()) {
            Some(j) => j,
            None => return String::new(),
        };
        if j == 0 {
            return String::new();
        }
        let alt = match row[alt_c].split(',').nth(j - 1) {
            Some(a) => a,
            None => return String::new(),
        };
        let mut variant = format!("{}-{}-{}-{}", row[chrom_c], row[pos_c], row[ref_c], alt);
        if let Some(syn) = variant_synonyms.get(&variant) {
            variant = syn.clone();
        }
        if !defining_variants.contains(&variant) {
            return String::new();
        }
        variant
    };

    let mut index: Vec<String> = Vec::new();
    let mut out_rows: Vec<Vec<String>> = Vec::new();

    for sample in vcf.samples() {
        let sample_c = vcf.columns.iter().position(|c| *c == sample).unwrap();
        let mut results: Vec<String> = Vec::new();
        let mut alt_phase: Vec<String> = Vec::new();
        let mut all_alleles: Vec<String> = Vec::new();

        for i in 0..3usize {
            let candidates: Vec<String> = if i == 2 {
                let phase_set: HashSet<String> = alt_phase.iter().cloned().collect();
                let cands = one_haplotype(&phase_set);
                let new: Vec<String> = cands
                    .into_iter()
                    .filter(|x| !all_alleles.contains(x))
                    .collect();
                for x in &new {
                    if !all_alleles.contains(x) {
                        all_alleles.push(x.clone());
                    }
                }
                all_alleles = core::sort_alleles(&all_alleles, "priority", Some(&gene), &assembly);
                new
            } else {
                let observed: Vec<String> = vcf
                    .rows
                    .iter()
                    .map(|r| one_row(r, sample_c, i))
                    .filter(|x| !x.is_empty())
                    .collect();
                for x in &observed {
                    if !alt_phase.contains(x) {
                        alt_phase.push(x.clone());
                    }
                }
                let observed_set: HashSet<String> = observed.into_iter().collect();
                let cands = one_haplotype(&observed_set);
                for x in &cands {
                    if !all_alleles.contains(x) {
                        all_alleles.push(x.clone());
                    }
                }
                cands
            };
            results.push(format!("{};", candidates.join(";")));
        }

        // VariantData column.
        let mut af_list: Vec<String> = Vec::new();
        for allele in &all_alleles {
            if *allele == default_allele {
                af_list.push(format!("{allele}:default"));
            } else {
                let idx = star_alleles.iter().position(|(a, _)| a == allele).unwrap();
                let vars = &star_alleles[idx].1;
                let variants = vars.join(",");
                let fractions: Vec<String> = vars
                    .iter()
                    .map(|x| {
                        let target = reformatted_variants.get(x).map(|s| s.as_str()).unwrap_or(x);
                        match vcf.get_af(&sample, target) {
                            Some(f) => python_float_str(f),
                            None => "nan".to_string(),
                        }
                    })
                    .collect();
                af_list.push(format!("{allele}:{variants}:{}", fractions.join(",")));
            }
        }
        results.push(format!("{};", af_list.join(";")));

        index.push(sample);
        out_rows.push(results);
    }

    let mut metadata = consolidated_variants.copy_metadata();
    for kv in metadata.iter_mut() {
        if kv.0 == "SemanticType" {
            kv.1 = "SampleTable[Alleles]".to_string();
        }
    }
    let table = SampleTable {
        index,
        columns: vec![
            "Haplotype1".to_string(),
            "Haplotype2".to_string(),
            "AlternativePhase".to_string(),
            "VariantData".to_string(),
        ],
        rows: out_rows,
    };
    Ok(Archive::new(metadata, ArchiveData::SampleTable(table)))
}

/// `call_phenotypes(genotypes)` → `SampleTable[Phenotypes]`.
pub fn call_phenotypes(genotypes: &Archive) -> Result<Archive, PgxError> {
    genotypes.check_type(&["SampleTable[Genotypes]"])?;
    let gene = genotypes.get("Gene").expect("Gene metadata").to_string();
    let t = genotypes.as_sample_table();
    let gt_c = t.columns.iter().position(|c| c == "Genotype").ok_or_else(|| {
        PgxError::IncorrectMetadata("SampleTable[Genotypes] missing 'Genotype' column".into())
    })?;

    let mut rows = Vec::new();
    for r in &t.rows {
        let genotype = &r[gt_c];
        let phenotype = if genotype == "Indeterminate" {
            "Indeterminate".to_string()
        } else {
            let mut it = genotype.split('/');
            let a1 = it.next().unwrap_or("");
            let a2 = it.next().unwrap_or("");
            core::predict_phenotype(&gene, a1, a2)
        };
        rows.push(vec![phenotype]);
    }

    let metadata = vec![
        ("Gene".to_string(), gene),
        (
            "SemanticType".to_string(),
            "SampleTable[Phenotypes]".to_string(),
        ),
    ];
    let table = SampleTable {
        index: t.index.clone(),
        columns: vec!["Phenotype".to_string()],
        rows,
    };
    Ok(Archive::new(metadata, ArchiveData::SampleTable(table)))
}

/// `combine_results(genotypes, phenotypes, alleles, cnv_calls)` →
/// `SampleTable[Results]`. Absent columns become empty (pandas NaN).
pub fn combine_results(
    genotypes: Option<&Archive>,
    phenotypes: Option<&Archive>,
    alleles: Option<&Archive>,
    cnv_calls: Option<&Archive>,
) -> Result<Archive, PgxError> {
    if let Some(a) = genotypes {
        a.check_type(&["SampleTable[Genotypes]"])?;
    }
    if let Some(a) = phenotypes {
        a.check_type(&["SampleTable[Phenotypes]"])?;
    }
    if let Some(a) = alleles {
        a.check_type(&["SampleTable[Alleles]"])?;
    }
    if let Some(a) = cnv_calls {
        a.check_type(&["SampleTable[CNVCalls]"])?;
    }

    let tables: Vec<&Archive> = [genotypes, phenotypes, alleles, cnv_calls]
        .into_iter()
        .flatten()
        .collect();
    if tables.is_empty() {
        return Err(PgxError::IncorrectMetadata("No input data detected".into()));
    }

    let mut metadata = Vec::new();
    for k in ["Gene", "Assembly"] {
        let vals: Vec<&str> = tables.iter().filter_map(|t| t.get(k)).collect();
        let uniq: HashSet<&str> = vals.iter().copied().collect();
        if uniq.len() > 1 {
            return Err(PgxError::IncorrectMetadata(format!(
                "Found incompatible inputs: {vals:?}"
            )));
        }
        if let Some(v) = vals.first() {
            metadata.push((k.to_string(), v.to_string()));
        }
    }

    // Use the first table's sample order as the row index.
    let index = tables[0].as_sample_table().index.clone();

    // For each output column, find the source table that provides it.
    const COLS: [&str; 7] = [
        "Genotype",
        "Phenotype",
        "Haplotype1",
        "Haplotype2",
        "AlternativePhase",
        "VariantData",
        "CNV",
    ];
    let lookup = |col: &str, sample: &str| -> String {
        for t in &tables {
            let st = t.as_sample_table();
            if let Some(ci) = st.columns.iter().position(|c| c == col) {
                if let Some(ri) = st.index.iter().position(|s| s == sample) {
                    return st.rows[ri][ci].clone();
                }
            }
        }
        String::new()
    };

    let rows: Vec<Vec<String>> = index
        .iter()
        .map(|sample| COLS.iter().map(|col| lookup(col, sample)).collect())
        .collect();

    metadata.push((
        "SemanticType".to_string(),
        "SampleTable[Results]".to_string(),
    ));
    let table = SampleTable {
        index,
        columns: COLS.iter().map(|s| s.to_string()).collect(),
        rows,
    };
    Ok(Archive::new(metadata, ArchiveData::SampleTable(table)))
}

/// `compare_genotypes(first, second, verbose)` — concordance report for the
/// `Genotype` and `CNV` columns. Returns the report text PyPGx would print.
pub fn compare_genotypes(first: &Archive, second: &Archive, verbose: bool) -> String {
    first
        .check_type(&["SampleTable[Results]"])
        .expect("first type");
    second
        .check_type(&["SampleTable[Results]"])
        .expect("second type");
    let (a, b) = (first.as_sample_table(), second.as_sample_table());

    let mut out = String::new();
    for col in ["Genotype", "CNV"] {
        out.push_str(&format!("# {col}\n"));
        let total = a.rows.len();
        out.push_str(&format!("Total: {total}\n"));

        let ac = a.columns.iter().position(|c| c == col);
        let bc = b.columns.iter().position(|c| c == col);
        // Align on shared samples; drop rows with a NaN (empty) in either.
        let mut pairs: Vec<(String, String, String)> = Vec::new();
        for (i, sample) in a.index.iter().enumerate() {
            let Some(bi) = b.index.iter().position(|s| s == sample) else {
                continue;
            };
            let va = ac.map(|c| a.rows[i][c].clone()).unwrap_or_default();
            let vb = bc.map(|c| b.rows[bi][c].clone()).unwrap_or_default();
            if va.is_empty() || vb.is_empty() {
                continue;
            }
            pairs.push((sample.clone(), va, vb));
        }
        out.push_str(&format!("Compared: {}\n", pairs.len()));
        let concordant = pairs.iter().filter(|(_, x, y)| x == y).count();
        if !pairs.is_empty() {
            out.push_str(&format!(
                "Concordance: {:.3} ({}/{})\n",
                concordant as f64 / pairs.len() as f64,
                concordant,
                pairs.len()
            ));
        } else {
            out.push_str("Concordance: N/A\n");
        }
        if verbose {
            out.push_str("Discordant:\n");
            let discordant: Vec<&(String, String, String)> =
                pairs.iter().filter(|(_, x, y)| x != y).collect();
            if discordant.is_empty() {
                out.push_str("None\n");
            } else {
                for (sample, x, y) in discordant {
                    out.push_str(&format!("{sample}\t{x}\t{y}\n"));
                }
            }
        }
    }
    out
}

/// `count_alleles(results)` — star-allele counts, ordered by name.
pub fn count_alleles(results: &Archive) -> Vec<(String, usize)> {
    results
        .check_type(&["SampleTable[Results]"])
        .expect("results type");
    let t = results.as_sample_table();
    let gt_c = t.columns.iter().position(|c| c == "Genotype").unwrap();

    let mut counts: HashMap<String, usize> = HashMap::new();
    for r in &t.rows {
        let alleles: Vec<String> = if r[gt_c] == "Indeterminate" {
            vec!["Indeterminate".to_string(), "Indeterminate".to_string()]
        } else {
            r[gt_c].split('/').map(|x| x.to_string()).collect()
        };
        for a in alleles {
            *counts.entry(a).or_insert(0) += 1;
        }
    }
    let names: Vec<String> = counts.keys().cloned().collect();
    let ordered = core::sort_alleles(&names, "name", None, "GRCh37");
    ordered
        .into_iter()
        .map(|a| (a.clone(), counts[&a]))
        .collect()
}

/// `import_variants(gene, vcf, ...)` — slice an input VcfFrame to the target
/// gene, drop duplicate variants, strip to `GT:AD:DP`, add `AF`, optionally
/// subset samples, and return a `VcfFrame[Imported]` (or `[Consolidated]` when
/// already fully phased) archive. Mirrors `utils.import_variants`.
///
/// Takes an in-memory `VcfFrame` (the `vcf=VcfFrame` path). The file path with
/// bgzf+tabix region access will route through noodles. `platform="LongRead"`
/// needs `_phase_extension` (not yet ported) — returns `NotPorted`.
pub fn import_variants(
    gene: &str,
    vcf: &VcfFrame,
    assembly: &str,
    platform: &str,
    samples: Option<&[String]>,
    exclude: bool,
) -> Result<Archive, PgxError> {
    let region = core::get_region(gene, assembly)?;
    let mut vf = vcf.slice(&region);
    // Drop duplicate variant records, keeping the first (PyPGx warns; we don't).
    vf = vf.drop_duplicates(&["CHROM", "POS", "REF", "ALT"]);
    vf = vf.update_chr_prefix("remove");
    vf = vf.strip("GT:AD:DP");
    vf = vf.add_af();
    if let Some(s) = samples {
        vf = vf.subset(s, exclude);
    }

    let semantic_type = if platform == "LongRead" {
        vf = phase_extension(&vf, gene, assembly);
        "VcfFrame[Consolidated]"
    } else if vf.phased() {
        "VcfFrame[Consolidated]"
    } else {
        vf = vf.unphase();
        "VcfFrame[Imported]"
    };

    // X-chromosome haploid male calls interfere downstream for G6PD.
    if gene == "G6PD" {
        vf = vf.diploidize();
    }

    let metadata = vec![
        ("Platform".to_string(), platform.to_string()),
        ("Gene".to_string(), gene.to_string()),
        ("Assembly".to_string(), assembly.to_string()),
        ("SemanticType".to_string(), semantic_type.to_string()),
    ];
    Ok(Archive::new(metadata, ArchiveData::Vcf(vf)))
}

/// `create_consolidated_vcf(imported, phased)` — merge a `VcfFrame[Imported]`
/// (genotype data: AD/DP/AF) with a `VcfFrame[Phased]` (Beagle output) into a
/// `VcfFrame[Consolidated]`: phased variants get their imported data appended,
/// imported-only variants are kept, and the union is phase-extended. Faithful
/// port of `utils.create_consolidated_vcf`.
pub fn create_consolidated_vcf(
    imported_variants: &Archive,
    phased_variants: &Archive,
) -> Result<Archive, PgxError> {
    imported_variants.check_type(&["VcfFrame[Imported]"])?;
    phased_variants.check_type(&["VcfFrame[Phased]"])?;
    let gene = imported_variants.get("Gene").expect("Gene");
    let assembly = imported_variants.get("Assembly").expect("Assembly");
    let platform = imported_variants.get("Platform").expect("Platform");
    assert_eq!(Some(gene), phased_variants.get("Gene"), "different genes");
    assert_eq!(
        Some(assembly),
        phased_variants.get("Assembly"),
        "different assemblies"
    );
    assert_eq!(
        Some(platform),
        phased_variants.get("Platform"),
        "different platforms"
    );

    let format = if platform == "WGS" || platform == "Targeted" {
        "GT:AD:DP:AF"
    } else {
        "GT"
    };
    let vf1 = imported_variants.as_vcf().strip(format);
    let vf2 = phased_variants.as_vcf().strip("GT");

    // For each phased variant, append its imported genotype data (minus the GT
    // field). INFO becomes 'Phased', FORMAT becomes `format`.
    let mut vf3_rows = Vec::with_capacity(vf2.rows.len());
    for r in &vf2.rows {
        let variant = format!("{}-{}-{}-{}", r[0], r[1], r[3], r[4]);
        let mut nr = r.clone();
        if let Some(s) = vf1.fetch(&variant) {
            for (col, cell) in nr.iter_mut().enumerate().skip(9) {
                let imported_minus_gt = s[col].split(':').skip(1).collect::<Vec<_>>().join(":");
                *cell = format!("{cell}:{imported_minus_gt}");
            }
        }
        nr[7] = "Phased".to_string();
        nr[8] = format.to_string();
        vf3_rows.push(nr);
    }

    // Imported variants absent from the phased set, then merge + sort.
    let vf4 = vf1.filter_vcf(&vf2, true);
    let mut vf5_rows = vf3_rows;
    vf5_rows.extend(vf4.rows);
    let vf5 = VcfFrame::new(Vec::new(), vf2.columns.clone(), vf5_rows).sort();
    let vf6 = phase_extension(&vf5, gene, assembly);

    let mut metadata = phased_variants.copy_metadata();
    if let Some(e) = metadata.iter_mut().find(|(k, _)| k == "SemanticType") {
        e.1 = "VcfFrame[Consolidated]".to_string();
    }
    Ok(Archive::new(metadata, ArchiveData::Vcf(vf6)))
}

/// `_phase_extension(vf, gene, assembly)` — estimate haplotype phase of variants
/// that read-backed phasing left unphased, by scoring each het call's two
/// orientations against already-phased "anchor" variants grouped by haplotype.
/// Faithful port of `utils._phase_extension`. Adds a `PE` FORMAT field carrying
/// the four anchor scores.
fn phase_extension(vf: &VcfFrame, gene: &str, assembly: &str) -> VcfFrame {
    let samples = vf.samples();

    // anchors[sample] = [haplotype-0 variants, haplotype-1 variants] from phased calls.
    let mut anchors: HashMap<String, [Vec<String>; 2]> = HashMap::new();
    for s in &samples {
        anchors.insert(s.clone(), [Vec::new(), Vec::new()]);
    }
    for r in &vf.rows {
        for allele in r[4].split(',') {
            let variant = format!("{}-{}-{}-{}", r[0], r[1], r[3], allele);
            for (si, s) in samples.iter().enumerate() {
                let gt = r[9 + si].split(':').next().unwrap_or("");
                if !gt.contains('|') {
                    continue;
                }
                let h: Vec<&str> = gt.split('|').collect();
                let a = anchors.get_mut(s).unwrap();
                if h[0] != "0" {
                    a[0].push(variant.clone());
                }
                if h[1] != "0" {
                    a[1].push(variant.clone());
                }
            }
        }
    }

    let variant_synonyms = core::get_variant_synonyms(gene, assembly);

    let mut new_rows = Vec::with_capacity(vf.rows.len());
    for r in &vf.rows {
        if crate::fuc::row_phased(&r[9..]) {
            new_rows.push(r.clone());
            continue;
        }
        let mut nr = r.clone();
        nr[8] = format!("{}:PE", r[8]);
        let alts: Vec<&str> = r[4].split(',').collect();
        for (si, s) in samples.iter().enumerate() {
            let cell = &r[9 + si];
            if !crate::fuc::gt_het(cell) {
                nr[9 + si] = format!("{}:0,0,0,0", crate::fuc::gt_pseudophase(cell));
                continue;
            }
            // scores[i][j]: best anchor overlap of called allele i with haplotype j.
            let mut scores = [[0i64, 0], [0i64, 0]];
            // gt_het accepts both '/' and '|', so split on either; a non-diploid
            // call can't be phase-extended, so pseudo-phase it and move on.
            let gt: Vec<&str> = cell.split(':').next().unwrap_or("").split(['/', '|']).collect();
            if gt.len() < 2 {
                nr[9 + si] = format!("{}:0,0,0,0", crate::fuc::gt_pseudophase(cell));
                continue;
            }
            for i in 0..2 {
                if gt[i] == "0" {
                    continue;
                }
                let alt_allele = alts[gt[i].parse::<usize>().unwrap() - 1];
                let mut variant = format!("{}-{}-{}-{}", r[0], r[1], r[3], alt_allele);
                if let Some(syn) = variant_synonyms.get(&variant) {
                    variant = syn.clone();
                }
                let star_alleles =
                    core::list_alleles(gene, Some(&[variant.clone()]), assembly);
                for (j, anchor_hap) in anchors[s].iter().enumerate() {
                    for star_allele in &star_alleles {
                        let sv =
                            core::list_variants(gene, Some(&[star_allele.clone()]), "all", assembly);
                        let score = anchor_hap.iter().filter(|x| sv.contains(x)).count() as i64;
                        if score > scores[i][j] {
                            scores[i][j] = score;
                        }
                    }
                }
            }
            let (a, b, c, d) = (scores[0][0], scores[0][1], scores[1][0], scores[1][1]);
            let flip = if a.max(b) == c.max(d) {
                (a < b && c > d) || (a == b && c > d) || (a < b && c == d)
            } else if a.max(b) > c.max(d) {
                !(a > b)
            } else {
                c > d
            };
            let result_gt = if flip {
                format!("{}|{}", gt[1], gt[0])
            } else {
                format!("{}|{}", gt[0], gt[1])
            };
            let rest: Vec<&str> = cell.split(':').skip(1).collect();
            nr[9 + si] = format!("{}:{}:{},{},{},{}", result_gt, rest.join(":"), a, b, c, d);
        }
        new_rows.push(nr);
    }
    VcfFrame::new(Vec::new(), vf.columns.clone(), new_rows)
}

/// `compute_copy_number(read_depth, control_statistics, samples_without_sv)` —
/// normalize per-position read depth into copy number. Pure (no BAM): intra-
/// sample normalization divides each sample's depth by its control median × 2;
/// for `Targeted` data an inter-sample (per-position) normalization follows.
/// Mirrors `utils.compute_copy_number`.
pub fn compute_copy_number(
    read_depth: &Archive,
    control_statistics: &Archive,
    samples_without_sv: Option<&[String]>,
) -> Result<Archive, PgxError> {
    read_depth.check_type(&["CovFrame[ReadDepth]"])?;
    control_statistics.check_type(&["SampleTable[Statistics]"])?;
    let cf = read_depth.as_cov();
    let stats = control_statistics.as_sample_table();
    let samples = cf.samples();
    if samples.iter().collect::<HashSet<_>>() != stats.index.iter().collect::<HashSet<_>>() {
        return Err(PgxError::External(
            "read-depth and control-statistics have different sample sets".into(),
        ));
    }

    let num = |v: &str| -> Result<f64, PgxError> {
        v.parse::<f64>()
            .map_err(|_| PgxError::External(format!("non-numeric depth value {v:?}")))
    };
    // Control median ('50%') per sample.
    let pct50 = stats.columns.iter().position(|c| c == "50%").ok_or_else(|| {
        PgxError::IncorrectMetadata("SampleTable[Statistics] missing '50%' column".into())
    })?;
    let mut medians: std::collections::HashMap<&str, f64> = std::collections::HashMap::new();
    for s in &samples {
        let i = stats.index.iter().position(|x| x == s).ok_or_else(|| {
            PgxError::External(format!("sample {s} missing from control statistics"))
        })?;
        medians.insert(s.as_str(), num(&stats.rows[i][pct50])?);
    }

    // Intra-sample normalization → a float matrix (rows × samples).
    let mut mat: Vec<Vec<f64>> = Vec::with_capacity(cf.rows.len());
    for r in &cf.rows {
        let mut out = Vec::with_capacity(samples.len());
        for (k, s) in samples.iter().enumerate() {
            out.push(num(&r[2 + k])? / medians[s.as_str()] * 2.0);
        }
        mat.push(out);
    }

    // Inter-sample (per-position) normalization for targeted sequencing.
    if read_depth.get("Platform") == Some("Targeted") {
        let cols: Vec<usize> = match samples_without_sv {
            Some(sw) => samples
                .iter()
                .enumerate()
                .filter(|(_, s)| sw.contains(s))
                .map(|(i, _)| i)
                .collect(),
            None => (0..samples.len()).collect(),
        };
        for row in &mut mat {
            let mut vals: Vec<f64> = cols.iter().map(|&i| row[i]).collect();
            let med = median_of(&mut vals);
            // `.replace(0, np.nan)`: a zero median makes the whole row NaN.
            for v in row.iter_mut() {
                *v = if med == 0.0 { f64::NAN } else { *v / med * 2.0 };
            }
        }
    }

    let rows = cf
        .rows
        .iter()
        .zip(&mat)
        .map(|(r, m)| {
            let mut nr = vec![r[0].clone(), r[1].clone()];
            nr.extend(m.iter().map(|&v| python_float_str(v)));
            nr
        })
        .collect();
    let cf_out = crate::fuc::CovFrame {
        columns: cf.columns.clone(),
        rows,
    };

    let mut metadata = read_depth.copy_metadata();
    if let Some(e) = metadata.iter_mut().find(|(k, _)| k == "SemanticType") {
        e.1 = "CovFrame[CopyNumber]".to_string();
    }
    let control = control_statistics.get("Control").unwrap_or("").to_string();
    metadata.push(("Control".to_string(), control));
    metadata.push((
        "Samples".to_string(),
        match samples_without_sv {
            None => "None".to_string(),
            Some(sw) => sw.join(","),
        },
    ));
    Ok(Archive::new(metadata, ArchiveData::Cov(cf_out)))
}

/// Median of a slice (pandas semantics: mean of the two middle values for an
/// even count). Mutates `vals` (sorts it).
fn median_of(vals: &mut [f64]) -> f64 {
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = vals.len();
    if n == 0 {
        f64::NAN
    } else if n % 2 == 1 {
        vals[n / 2]
    } else {
        (vals[n / 2 - 1] + vals[n / 2]) / 2.0
    }
}

// ---- CNV calling (RBF OvR-SVM) -------------------------------------------
// predict/test are pure (verified decision function + median_filter); only
// train needs sklears. PyPGx's pickled `Model[CNV]` is converted to the Rust
// `data.json` form once via tools/convert_cnv_model.py.

/// `_process_copy_number` — gap-fill missing positions (ffill/bfill) then apply
/// a 1000-wide median filter per sample, matching `utils._process_copy_number`.
fn process_copy_number(copy_number: &Archive) -> Archive {
    let cf = copy_number.as_cov();
    let gene = copy_number.get("Gene").expect("Gene");
    let assembly = copy_number.get("Assembly").expect("Assembly");
    let region = core::get_region(gene, assembly).expect("region");
    let (_, start, end) = crate::fuc::parse_region(&region);
    let (start, end) = (start.expect("start"), end.expect("end"));
    let nsamples = cf.samples().len();

    // Densify positions if the frame is sparse vs the region span. An empty
    // frame has nothing to densify (and no rows[0]), so skip it.
    let mut rows = cf.rows.clone();
    if !rows.is_empty() && end - start + 1 > rows.len() as i64 {
        let first: i64 = rows[0][1].parse().unwrap();
        let last: i64 = rows.last().unwrap()[1].parse().unwrap();
        let by_pos: std::collections::HashMap<i64, &Vec<String>> =
            cf.rows.iter().map(|r| (r[1].parse::<i64>().unwrap(), r)).collect();
        let mut filled: Vec<Vec<String>> = Vec::new();
        // ffill/bfill seed: first known row.
        let mut last_known: Option<Vec<String>> = None;
        for pos in (first - 1)..=last {
            match by_pos.get(&pos) {
                Some(r) => {
                    last_known = Some((*r).clone());
                    filled.push((*r).clone());
                }
                None => {
                    // forward-fill; back-fill handled in a second pass.
                    let mut nr = vec![String::new(), pos.to_string()];
                    nr.resize(2 + nsamples, String::new());
                    filled.push(nr);
                }
            }
        }
        // forward fill (chrom + sample cols), then back fill leading gaps.
        let _ = last_known;
        for col in (0..2 + nsamples).filter(|&c| c != 1) {
            let mut prev = String::new();
            for r in filled.iter_mut() {
                if r[col].is_empty() {
                    r[col] = prev.clone();
                } else {
                    prev = r[col].clone();
                }
            }
            let mut next = String::new();
            for r in filled.iter_mut().rev() {
                if r[col].is_empty() {
                    r[col] = next.clone();
                } else {
                    next = r[col].clone();
                }
            }
        }
        rows = filled;
    }

    // Per-sample 1000-wide median filter.
    let filtered: Vec<Vec<f64>> = (0..nsamples)
        .map(|s| {
            let col: Vec<f64> = rows.iter().map(|r| r[2 + s].parse::<f64>().unwrap()).collect();
            crate::cnv::median_filter(&col, 1000)
        })
        .collect();
    let new_rows = rows
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut nr = vec![r[0].clone(), r[1].clone()];
            nr.extend((0..nsamples).map(|s| python_float_str(filtered[s][i])));
            nr
        })
        .collect();
    Archive::new(
        copy_number.copy_metadata(),
        ArchiveData::Cov(crate::fuc::CovFrame {
            columns: cf.columns.clone(),
            rows: new_rows,
        }),
    )
}

/// CNV code → star-allele Name for `gene`, ordered as in the cnv-table.
fn cnv_names(gene: &str) -> Vec<String> {
    let t = core::load_cnv_table();
    let (gc, nc) = (t.col("Gene"), t.col("Name"));
    t.rows
        .iter()
        .filter(|r| r[gc].as_str() == Some(gene))
        .filter_map(|r| r[nc].as_str().map(|s| s.to_string()))
        .collect()
}

/// Resolve the `pypgx-bundle` root: `$PYPGX_BUNDLE`, else `$HOME/pypgx-bundle`
/// (mirrors PyPGx's `sdk.get_bundle_path`).
pub fn bundle_path() -> Result<String, PgxError> {
    if let Ok(p) = std::env::var("PYPGX_BUNDLE") {
        return Ok(p);
    }
    if let Ok(home) = std::env::var("HOME") {
        return Ok(format!("{home}/pypgx-bundle"));
    }
    Err(PgxError::BundleNotFound(
        "set $PYPGX_BUNDLE (or pass an explicit cnv_caller)".into(),
    ))
}

/// `predict_cnv(copy_number, cnv_caller)` — predict per-sample CNV calls with a
/// `Model[CNV]`, producing a `SampleTable[CNVCalls]`. Pure (RBF decision
/// function). When `cnv_caller=None`, the default model for this gene/assembly is
/// loaded from the bundle (`{bundle}/cnv/{assembly}/{gene}.zip`, a Rust-converted
/// `Model[CNV]` — see `tools/convert_cnv_models_all.py`), mirroring PyPGx.
pub fn predict_cnv(copy_number: &Archive, cnv_caller: Option<&Archive>) -> Result<Archive, PgxError> {
    copy_number.check_type(&["CovFrame[CopyNumber]"])?;
    let gene = copy_number.get("Gene").expect("Gene").to_string();

    // Resolve the caller: an explicit model, else the gene/assembly default from
    // the bundle. `loaded` outlives the borrow taken in the `None` arm.
    let loaded;
    let caller = match cnv_caller {
        Some(c) => c,
        None => {
            let assembly = copy_number.get("Assembly").expect("Assembly");
            let path = format!("{}/cnv/{assembly}/{gene}.zip", bundle_path()?);
            loaded = Archive::from_file(&path).map_err(|e| {
                PgxError::BundleNotFound(format!("default CNV model {path}: {e}"))
            })?;
            &loaded
        }
    };
    caller.check_type(&["Model[CNV]"])?;
    let model = caller.as_model();

    let processed = process_copy_number(copy_number);
    let cf = processed.as_cov();
    let names = cnv_names(&gene);
    let samples = cf.samples();
    let mut calls: Vec<Vec<String>> = Vec::with_capacity(samples.len());
    for si in 0..samples.len() {
        let mut x: Vec<f64> = Vec::with_capacity(cf.rows.len());
        for r in &cf.rows {
            x.push(r[2 + si].parse::<f64>().map_err(|_| {
                PgxError::External(format!("non-numeric copy-number value {:?}", r[2 + si]))
            })?);
        }
        let label = model.predict(&x) as usize;
        let name = names.get(label).ok_or_else(|| {
            PgxError::External(format!("CNV model predicted label {label} with no name for {gene}"))
        })?;
        calls.push(vec![name.clone()]);
    }

    let mut metadata = copy_number.copy_metadata();
    if let Some(e) = metadata.iter_mut().find(|(k, _)| k == "SemanticType") {
        e.1 = "SampleTable[CNVCalls]".to_string();
    }
    Ok(Archive::new(
        metadata,
        ArchiveData::SampleTable(SampleTable {
            index: samples,
            columns: vec!["CNV".to_string()],
            rows: calls,
        }),
    ))
}

/// Result of `test_cnv_caller`: accuracy plus the label-ordered confusion matrix.
#[derive(Clone, Debug)]
pub struct CnvTestReport {
    pub accuracy: f64,
    pub correct: usize,
    pub total: usize,
    pub labels: Vec<String>,
    pub confusion: Vec<Vec<usize>>,
}

/// `test_cnv_caller(cnv_caller, copy_number, cnv_calls)` — predict with the model
/// and score against known calls. Pure (`predict` + a confusion matrix).
pub fn test_cnv_caller(
    cnv_caller: &Archive,
    copy_number: &Archive,
    cnv_calls: &Archive,
) -> Result<CnvTestReport, PgxError> {
    cnv_caller.check_type(&["Model[CNV]"])?;
    cnv_calls.check_type(&["SampleTable[CNVCalls]"])?;
    let gene = copy_number.get("Gene").expect("Gene").to_string();
    let predicted = predict_cnv(copy_number, Some(cnv_caller))?;
    let pt = predicted.as_sample_table();
    let at = cnv_calls.as_sample_table();

    let labels = cnv_names(&gene);
    let label_idx = |n: &str| labels.iter().position(|l| l == n).unwrap();
    let mut confusion = vec![vec![0usize; labels.len()]; labels.len()];
    let (mut correct, mut total) = (0usize, 0usize);
    for (i, sample) in pt.index.iter().enumerate() {
        let pred = &pt.rows[i][0];
        let actual = &at.loc(sample)[0];
        confusion[label_idx(actual)][label_idx(pred)] += 1;
        total += 1;
        if pred == actual {
            correct += 1;
        }
    }
    Ok(CnvTestReport {
        accuracy: if total == 0 { 0.0 } else { correct as f64 / total as f64 },
        correct,
        total,
        labels,
        confusion,
    })
}

// ---- BAM depth (feature `bam`; depth engine = samtools-rs `native::depth`) ----
// NOTE: byte-parity with PyPGx (pysam) is unverified here — no BAM fixtures /
// samtools to diff against. The orchestration + the pure `describe` are faithful.

/// `compute_target_depth(gene, bams, assembly, bed)` — per-position read depth
/// over the target gene from BAMs → `CovFrame[ReadDepth]`.
#[cfg(feature = "bam")]
pub fn compute_target_depth(
    gene: &str,
    bams: &[String],
    assembly: &str,
    bed: Option<&str>,
) -> Result<Archive, Box<dyn std::error::Error>> {
    let region = core::get_region(gene, assembly)?;
    let cf = crate::fuc::CovFrame::from_bam(bams, &region, true)?;
    let platform = if bed.is_some() { "Targeted" } else { "WGS" };
    let metadata = vec![
        ("Gene".to_string(), gene.to_string()),
        ("Assembly".to_string(), assembly.to_string()),
        ("SemanticType".to_string(), "CovFrame[ReadDepth]".to_string()),
        ("Platform".to_string(), platform.to_string()),
    ];
    Ok(Archive::new(metadata, ArchiveData::Cov(cf)))
}

/// `compute_control_statistics(gene, bams, assembly, bed)` — per-sample depth
/// summary over a control region → `SampleTable[Statistics]`. `gene` may be a
/// known gene name or a raw `chrom:start-end` region.
#[cfg(feature = "bam")]
pub fn compute_control_statistics(
    gene: &str,
    bams: &[String],
    assembly: &str,
    bed: Option<&str>,
) -> Result<Archive, Box<dyn std::error::Error>> {
    let region = if core::list_genes("all").iter().any(|g| g == gene) {
        core::get_region(gene, assembly)?
    } else {
        gene.to_string()
    };
    let cf = crate::fuc::CovFrame::from_bam(bams, &region, false)?;
    let stats = describe_cov(&cf);
    let platform = if bed.is_some() { "Targeted" } else { "WGS" };
    let metadata = vec![
        ("Control".to_string(), gene.to_string()),
        ("Assembly".to_string(), assembly.to_string()),
        ("SemanticType".to_string(), "SampleTable[Statistics]".to_string()),
        ("Platform".to_string(), platform.to_string()),
    ];
    Ok(Archive::new(metadata, ArchiveData::SampleTable(stats)))
}

/// `prepare_depth_of_coverage(bams, ...)` — genome-wide (merged SV-gene regions)
/// per-position depth → `CovFrame[DepthOfCoverage]`.
#[cfg(feature = "bam")]
pub fn prepare_depth_of_coverage(
    bams: &[String],
    assembly: &str,
    bed: Option<&str>,
    genes: Option<&[String]>,
    exclude: bool,
) -> Result<Archive, Box<dyn std::error::Error>> {
    let regions =
        create_regions_bed(assembly, false, true, false, true, false, genes, exclude).to_regions();
    let mut columns = vec!["Chromosome".to_string(), "Position".to_string()];
    let mut rows = Vec::new();
    for (i, region) in regions.iter().enumerate() {
        let cf = crate::fuc::CovFrame::from_bam(bams, region, true)?;
        if i == 0 {
            columns = cf.columns;
        }
        rows.extend(cf.rows);
    }
    let platform = if bed.is_some() { "Targeted" } else { "WGS" };
    let metadata = vec![
        ("Assembly".to_string(), assembly.to_string()),
        (
            "SemanticType".to_string(),
            "CovFrame[DepthOfCoverage]".to_string(),
        ),
        ("Platform".to_string(), platform.to_string()),
    ];
    Ok(Archive::new(
        metadata,
        ArchiveData::Cov(crate::fuc::CovFrame { columns, rows }),
    ))
}

/// `cf.df.iloc[:, 2:].describe().T` — per-sample summary stats (count, mean,
/// std[ddof=1], min, 25/50/75%, max) as a `SampleTable[Statistics]`.
#[cfg(feature = "bam")]
fn describe_cov(cf: &crate::fuc::CovFrame) -> SampleTable {
    let columns: Vec<String> = ["count", "mean", "std", "min", "25%", "50%", "75%", "max"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let samples = cf.samples();
    let mut rows = Vec::new();
    for s in &samples {
        let ci = cf.columns.iter().position(|c| c == s).unwrap();
        let mut vals: Vec<f64> = cf.rows.iter().filter_map(|r| r[ci].parse().ok()).collect();
        rows.push(describe_stats(&mut vals).iter().map(|v| python_float_str(*v)).collect());
    }
    SampleTable {
        index: samples,
        columns,
        rows,
    }
}

/// pandas `Series.describe()` values: [count, mean, std(ddof=1), min, 25%, 50%,
/// 75%, max] with linear-interpolation percentiles. Mutates `vals` (sorts it).
#[cfg(feature = "bam")]
fn describe_stats(vals: &mut [f64]) -> [f64; 8] {
    let n = vals.len();
    if n == 0 {
        return [0.0; 8];
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let count = n as f64;
    let mean = vals.iter().sum::<f64>() / count;
    let std = if n > 1 {
        (vals.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (count - 1.0)).sqrt()
    } else {
        f64::NAN
    };
    let pct = |q: f64| -> f64 {
        let idx = q * (n as f64 - 1.0);
        let lo = idx.floor() as usize;
        let hi = idx.ceil() as usize;
        vals[lo] + (idx - lo as f64) * (vals[hi] - vals[lo])
    };
    [count, mean, std, vals[0], pct(0.25), pct(0.5), pct(0.75), vals[n - 1]]
}

/// `import_read_depth(gene, depth_of_coverage, ...)` — slice a
/// `CovFrame[DepthOfCoverage]` to the target gene and (optionally) subset
/// samples, producing a `CovFrame[ReadDepth]`. Pure (no BAM reading); mirrors
/// `utils.import_read_depth`.
pub fn import_read_depth(
    gene: &str,
    depth_of_coverage: &Archive,
    samples: Option<&[String]>,
    exclude: bool,
) -> Result<Archive, PgxError> {
    depth_of_coverage.check_type(&["CovFrame[DepthOfCoverage]"])?;

    let mut metadata = depth_of_coverage.copy_metadata();
    metadata.push(("Gene".to_string(), gene.to_string()));
    if let Some(e) = metadata.iter_mut().find(|(k, _)| k == "SemanticType") {
        e.1 = "CovFrame[ReadDepth]".to_string();
    }
    let assembly = metadata
        .iter()
        .find(|(k, _)| k == "Assembly")
        .map(|(_, v)| v.clone())
        .expect("Assembly");

    let region = core::get_region(gene, &assembly)?;
    let mut cf = depth_of_coverage.as_cov().update_chr_prefix("remove").slice(&region);
    if let Some(s) = samples {
        cf = cf.subset(s, exclude);
    }
    Ok(Archive::new(metadata, ArchiveData::Cov(cf)))
}

/// `filter_samples(archive, samples, exclude)` — subset an archive to the given
/// samples. Pure. Handles `VcfFrame[*]` (column subset) and `SampleTable[*]`
/// (row subset); `CovFrame[*]` is handled once that payload is ported (until
/// then the data is returned unchanged).
pub fn filter_samples(archive: &Archive, samples: &[String], exclude: bool) -> Archive {
    let st = archive.semantic_type();
    let data = if st.contains("VcfFrame") {
        ArchiveData::Vcf(archive.as_vcf().subset(samples, exclude))
    } else if st.contains("CovFrame") {
        ArchiveData::Cov(archive.as_cov().subset(samples, exclude))
    } else if st.contains("SampleTable") {
        let t = archive.as_sample_table();
        let idxs: Vec<usize> = if exclude {
            (0..t.index.len())
                .filter(|&i| !samples.contains(&t.index[i]))
                .collect()
        } else {
            // Include in the requested order (matching pandas `.loc[samples]`).
            samples
                .iter()
                .filter_map(|s| t.index.iter().position(|x| x == s))
                .collect()
        };
        ArchiveData::SampleTable(SampleTable {
            index: idxs.iter().map(|&i| t.index[i].clone()).collect(),
            columns: t.columns.clone(),
            rows: idxs.iter().map(|&i| t.rows[i].clone()).collect(),
        })
    } else {
        archive.data.clone()
    };
    Archive::new(archive.copy_metadata(), data)
}

/// `create_regions_bed(...)` — the BED of all regions PyPGx uses, built from the
/// gene table. Pure (no external tool); mirrors `utils.create_regions_bed`.
///
/// Filters are applied in PyPGx's order: `genes`/`exclude`, then `target_genes`,
/// `sv_genes`, `var_genes`. The result is chromosome-ordered (gene-table order
/// preserved within a chromosome); `add_chr_prefix` and `merge` post-process it.
#[allow(clippy::too_many_arguments)]
pub fn create_regions_bed(
    assembly: &str,
    add_chr_prefix: bool,
    merge: bool,
    target_genes: bool,
    sv_genes: bool,
    var_genes: bool,
    genes: Option<&[String]>,
    exclude: bool,
) -> BedFrame {
    let df = core::load_gene_table();
    let region_col = format!("{assembly}Region");
    let (gene_c, region_c) = (df.col("Gene"), df.col(&region_col));
    let (target_c, sv_c, var_c) = (df.col("Target"), df.col("SV"), df.col("Variants"));

    let mut data: Vec<(String, i64, i64, String)> = Vec::new();
    for r in &df.rows {
        let gene = r[gene_c].as_str().expect("gene name").to_string();
        if let Some(gs) = genes {
            let in_list = gs.iter().any(|g| g == &gene);
            // `exclude` is ignored when `genes` is None (matched above).
            if exclude == in_list {
                continue;
            }
        }
        if (target_genes && !r[target_c].is_true())
            || (sv_genes && !r[sv_c].is_true())
            || (var_genes && !r[var_c].is_true())
        {
            continue;
        }
        let region = r[region_c].as_str().expect("region");
        let (chrom, start, end) = crate::fuc::parse_region(region);
        data.push((chrom, start.expect("start"), end.expect("end"), gene));
    }

    let mut bf = BedFrame::from_regions(data);
    if add_chr_prefix {
        bf.add_chr_prefix();
    }
    if merge {
        bf = bf.merge();
    }
    bf
}
