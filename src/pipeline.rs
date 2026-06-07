//! Port of `pypgx.api.pipeline` — the end-to-end pipelines that chain the
//! per-step API functions and write each intermediate archive to an output dir.
//!
//! `run_long_read_pipeline` is fully native (read-backed phasing, no externals).
//! `run_chip_pipeline` is native when the input is already fully phased; the
//! unphased branch calls `estimate_phase_beagle`, which is deferred to
//! `beagle-rs` (returns `PgxError::NotPorted`). `run_ngs_pipeline` additionally
//! needs depth/CNV + Beagle.

use std::path::Path;

use crate::fuc::VcfFrame;
use crate::sdk::{Archive, ArchiveData, PgxError};
use crate::{api, core, external, genotype};

type PipelineResult = Result<(), Box<dyn std::error::Error>>;

/// Run a pipeline body, converting any panic into a clean `Err` so a single bad
/// gene (or any unexpected failure) never aborts the caller. Mirrors how PyPGx
/// raises rather than crashes the process; here we additionally catch panics
/// from not-yet-graceful internal paths.
fn guard<F: FnOnce() -> PipelineResult>(body: F) -> PipelineResult {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(body)) {
        Ok(result) => result,
        Err(payload) => {
            let msg = payload
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| payload.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "unknown panic".to_string());
            Err(Box::new(PgxError::External(format!("pipeline panicked: {msg}"))))
        }
    }
}

/// `os.mkdir(output)` with PyPGx's `force` semantics (rmtree first if forcing).
fn ensure_output(output: &str, force: bool) -> std::io::Result<()> {
    let p = Path::new(output);
    if p.exists() && force {
        std::fs::remove_dir_all(p)?;
    }
    std::fs::create_dir(p)
}

/// The shared tail: consolidated variants → alleles → genotypes → phenotypes →
/// results, writing each archive. Mirrors the identical suffix of every pipeline.
fn finish(output: &str, consolidated: &Archive) -> PipelineResult {
    let alleles = api::predict_alleles(consolidated)?;
    alleles.to_file(&format!("{output}/alleles.zip"))?;
    let genotypes = genotype::call_genotypes(Some(&alleles), None)?;
    genotypes.to_file(&format!("{output}/genotypes.zip"))?;
    let phenotypes = api::call_phenotypes(&genotypes)?;
    phenotypes.to_file(&format!("{output}/phenotypes.zip"))?;
    let results = api::combine_results(Some(&genotypes), Some(&phenotypes), Some(&alleles), None)?;
    results.to_file(&format!("{output}/results.zip"))?;
    Ok(())
}

/// `run_long_read_pipeline` — read-backed phasing path; fully native.
#[allow(clippy::too_many_arguments)]
pub fn run_long_read_pipeline(
    gene: &str,
    output: &str,
    variants: &VcfFrame,
    assembly: &str,
    force: bool,
    samples: Option<&[String]>,
    exclude: bool,
) -> PipelineResult {
    guard(move || {
        if !core::is_target_gene(gene) {
            return Err(Box::new(PgxError::NotTargetGene(gene.to_string())));
        }
        ensure_output(output, force)?;
        let consolidated =
            api::import_variants(gene, variants, assembly, "LongRead", samples, exclude)?;
        consolidated.to_file(&format!("{output}/consolidated-variants.zip"))?;
        finish(output, &consolidated)
    })
}

/// `run_chip_pipeline` — chip/array genotypes. Native when the input is already
/// fully phased; otherwise statistical phasing (`estimate_phase_beagle`) is
/// required and the call surfaces `PgxError::NotPorted` until `beagle-rs` lands.
#[allow(clippy::too_many_arguments)]
pub fn run_chip_pipeline(
    gene: &str,
    output: &str,
    variants: &VcfFrame,
    assembly: &str,
    panel: Option<&str>,
    impute: bool,
    force: bool,
    samples: Option<&[String]>,
    exclude: bool,
) -> PipelineResult {
    guard(move || {
        if !core::is_target_gene(gene) {
            return Err(Box::new(PgxError::NotTargetGene(gene.to_string())));
        }
        ensure_output(output, force)?;
        let imported = api::import_variants(gene, variants, assembly, "Chip", samples, exclude)?;
        imported.to_file(&format!("{output}/imported-variants.zip"))?;

        let consolidated = if imported.semantic_type() == "VcfFrame[Consolidated]" {
            imported
        } else {
            let phased = external::estimate_phase_beagle(&imported, panel, impute)?;
            phased.to_file(&format!("{output}/phased-variants.zip"))?;
            let consolidated = api::create_consolidated_vcf(&imported, &phased)?;
            consolidated.to_file(&format!("{output}/consolidated-variants.zip"))?;
            consolidated
        };
        finish(output, &consolidated)
    })
}

/// `run_ngs_pipeline` — the full NGS pipeline. The variant arm runs natively
/// when the input is pre-phased (`VcfFrame[Consolidated]`) or for haploid
/// MT-RNR1; the general unphased arm needs Beagle (`estimate_phase_beagle`).
/// The SV/CNV arm runs depth import + copy-number natively but `predict_cnv`
/// needs the sklearn model. So this orchestrates end-to-end and surfaces
/// `NotPorted` exactly where Beagle / the CNV model are still required.
#[allow(clippy::too_many_arguments)]
pub fn run_ngs_pipeline(
    gene: &str,
    output: &str,
    variants: Option<&VcfFrame>,
    depth_of_coverage: Option<&Archive>,
    control_statistics: Option<&Archive>,
    platform: &str,
    assembly: &str,
    panel: Option<&str>,
    force: bool,
    samples: Option<&[String]>,
    exclude: bool,
    samples_without_sv: Option<&[String]>,
    cnv_caller: Option<&Archive>,
) -> PipelineResult {
    guard(move || {
    if !core::is_target_gene(gene) {
        return Err(Box::new(PgxError::NotTargetGene(gene.to_string())));
    }
    let gt = core::load_gene_table();
    let row = gt
        .rows
        .iter()
        .find(|r| r[gt.col("Gene")].as_str() == Some(gene))
        .expect("gene row");
    let small_var = row[gt.col("Variants")].is_true();
    let large_var = row[gt.col("SV")].is_true();

    ensure_output(output, force)?;
    let mut alleles: Option<Archive> = None;
    let mut cnv_calls: Option<Archive> = None;

    // ---- SNV/indel arm ----
    if small_var && variants.is_some() {
        let imported =
            api::import_variants(gene, variants.unwrap(), assembly, platform, samples, exclude)?;
        imported.to_file(&format!("{output}/imported-variants.zip"))?;
        let consolidated = if imported.semantic_type() == "VcfFrame[Consolidated]" {
            imported.clone()
        } else if gene == "MT-RNR1" {
            // Haploid variants need no phasing — pseudophase in place.
            let mut metadata = imported.copy_metadata();
            if let Some(e) = metadata.iter_mut().find(|(k, _)| k == "SemanticType") {
                e.1 = "VcfFrame[Consolidated]".to_string();
            }
            let c = Archive::new(metadata, ArchiveData::Vcf(imported.as_vcf().pseudophase()));
            c.to_file(&format!("{output}/consolidated-variants.zip"))?;
            c
        } else {
            let phased = external::estimate_phase_beagle(&imported, panel, false)?;
            phased.to_file(&format!("{output}/phased-variants.zip"))?;
            let c = api::create_consolidated_vcf(&imported, &phased)?;
            c.to_file(&format!("{output}/consolidated-variants.zip"))?;
            c
        };
        let a = api::predict_alleles(&consolidated)?;
        a.to_file(&format!("{output}/alleles.zip"))?;
        #[cfg(feature = "plots")]
        {
            let dir = format!("{output}/allele-fraction-profile");
            std::fs::create_dir(&dir)?;
            crate::plot::plot_vcf_allele_fraction(&imported, Some(&dir), None)?;
        }
        alleles = Some(a);
    }

    // ---- SV/CNV arm ----
    if large_var && depth_of_coverage.is_some() {
        let doc = depth_of_coverage.unwrap();
        doc.check_type(&["CovFrame[DepthOfCoverage]"])?;
        let cs = control_statistics
            .ok_or_else(|| Box::<dyn std::error::Error>::from("SV detection requires SampleTable[Statistics]"))?;
        let cs = match samples {
            Some(s) => api::filter_samples(cs, s, exclude),
            None => cs.clone(),
        };
        let read_depth = api::import_read_depth(gene, doc, samples, exclude)?;
        read_depth.to_file(&format!("{output}/read-depth.zip"))?;
        let copy_number = api::compute_copy_number(&read_depth, &cs, samples_without_sv)?;
        copy_number.to_file(&format!("{output}/copy-number.zip"))?;
        let calls = api::predict_cnv(&copy_number, cnv_caller)?;
        calls.to_file(&format!("{output}/cnv-calls.zip"))?;
        #[cfg(feature = "plots")]
        {
            let dir = format!("{output}/copy-number-profile");
            std::fs::create_dir(&dir)?;
            crate::plot::plot_bam_copy_number(&copy_number, Some(&dir), None)?;
        }
        cnv_calls = Some(calls);
    }

    let genotypes = genotype::call_genotypes(alleles.as_ref(), cnv_calls.as_ref())?;
    genotypes.to_file(&format!("{output}/genotypes.zip"))?;
    let phenotypes = api::call_phenotypes(&genotypes)?;
    phenotypes.to_file(&format!("{output}/phenotypes.zip"))?;
    let results =
        api::combine_results(Some(&genotypes), Some(&phenotypes), alleles.as_ref(), cnv_calls.as_ref())?;
    results.to_file(&format!("{output}/results.zip"))?;
    Ok(())
    })
}
