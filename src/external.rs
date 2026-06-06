//! Functions from `pypgx.api.utils` / `pipeline` / `plot` that depend on
//! external programs (Beagle, samtools, bcftools, the `pypgx-bundle` panels),
//! scikit-learn (the pickled `Model[CNV]`), or matplotlib.
//!
//! Per the project's scope decision ("full parity, **defer externals**"), these
//! keep PyPGx's exact signatures but are not implemented in this environment.
//! Each returns [`PgxError::NotPorted`] naming the missing dependency, and the
//! doc comment records what it would do (and, where relevant, the exact command
//! PyPGx runs) so the wrappers can be filled in once the tools are available.

use crate::sdk::{Archive, PgxError};

macro_rules! deferred {
    ($name:ident ( $($arg:ident : $ty:ty),* $(,)? ) -> $ret:ty, $dep:expr, $doc:expr) => {
        #[doc = $doc]
        #[allow(unused_variables)]
        pub fn $name($($arg : $ty),*) -> Result<$ret, PgxError> {
            Err(PgxError::NotPorted($dep.to_string()))
        }
    };
}

// ---- Beagle phasing ------------------------------------------------------

deferred!(
    estimate_phase_beagle(imported_variants: &Archive, panel: Option<&str>, impute: bool) -> Archive,
    "java + beagle.22Jul22.46e.jar + pypgx-bundle 1KGP panel",
    "`estimate_phase_beagle` — statistical phasing of a `VcfFrame[Imported]` \
     archive into `VcfFrame[Phased]`. PyPGx runs: `java -Xmx2g -jar beagle.jar \
     gt=input.vcf chrom=<region> ref=<panel> out=output impute=<bool> em=<bool>` \
     with an EM-skip fallback. Needs the 1KGP reference panel from pypgx-bundle."
);

// ---- BAM / depth (samtools) ---------------------------------------------

deferred!(
    slice_bam(input: &str, output: &str, assembly: &str, genes: Option<&[String]>, exclude: bool) -> (),
    "samtools (via pysam/pybam)",
    "`slice_bam` — subset a BAM to PyPGx's gene regions."
);
deferred!(
    prepare_depth_of_coverage(bams: &[String], assembly: &str, bed: Option<&str>, genes: Option<&[String]>, exclude: bool) -> Archive,
    "samtools (CovFrame.from_bam)",
    "`prepare_depth_of_coverage` — per-base depth over SV gene regions \
     (`CovFrame[DepthOfCoverage]`)."
);
deferred!(
    compute_target_depth(gene: &str, bams: &[String], assembly: &str, bed: Option<&str>) -> Archive,
    "samtools (CovFrame.from_bam)",
    "`compute_target_depth` — read depth for one target gene (`CovFrame[ReadDepth]`)."
);
deferred!(
    compute_control_statistics(gene: &str, bams: &[String], assembly: &str, bed: Option<&str>) -> Archive,
    "samtools (CovFrame.from_bam)",
    "`compute_control_statistics` — per-sample depth statistics over a control gene \
     (`SampleTable[Statistics]`)."
);
deferred!(
    import_read_depth(depth_of_coverage: &Archive, gene: &str) -> Archive,
    "pycov.CovFrame",
    "`import_read_depth` — slice a depth-of-coverage archive to one gene."
);

// ---- VCF construction (bcftools / pyvcf) ---------------------------------

deferred!(
    import_variants(vcf: &str, gene: &str, assembly: &str, platform: &str) -> Archive,
    "bcftools / pyvcf VCF I/O",
    "`import_variants` — read a VCF region into a `VcfFrame[Imported]` archive."
);
deferred!(
    create_input_vcf(vcf: &str, fasta: &str, bams: &[String], assembly: &str) -> Archive,
    "bcftools (variant calling)",
    "`create_input_vcf` — call SNVs/indels from BAMs into a `VcfFrame[Imported]`."
);
deferred!(
    create_consolidated_vcf(imported_variants: &Archive, phased_variants: &Archive) -> Archive,
    "pyvcf VcfFrame machinery",
    "`create_consolidated_vcf` — merge imported + phased variants into \
     `VcfFrame[Consolidated]` (the input to `predict_alleles`)."
);
deferred!(
    filter_samples(archive: &Archive, samples: &[String], exclude: bool) -> Archive,
    "pyvcf/pycov subsetting",
    "`filter_samples` — subset an archive to the given samples."
);

// ---- CNV (scikit-learn) --------------------------------------------------

deferred!(
    compute_copy_number(read_depth: &Archive, control_statistics: &Archive, samples_without_sv: Option<&[String]>) -> Archive,
    "pycov/numpy depth normalization",
    "`compute_copy_number` — normalize read depth into `CovFrame[CopyNumber]`."
);
deferred!(
    predict_cnv(copy_number: &Archive, cnv_caller: Option<&Archive>) -> Archive,
    "scikit-learn (pickled Model[CNV])",
    "`predict_cnv` — predict CNV calls with the SVM model. The pickled sklearn \
     model cannot be loaded in Rust; needs retraining to a native format."
);
deferred!(
    train_cnv_caller(copy_number: &Archive, cnv_calls: &Archive, confusion_matrix: Option<&str>) -> Archive,
    "scikit-learn (OneVsRestClassifier/SVC)",
    "`train_cnv_caller` — train a `Model[CNV]` SVM classifier."
);
deferred!(
    test_cnv_caller(cnv_caller: &Archive, copy_number: &Archive, cnv_calls: &Archive) -> (),
    "scikit-learn",
    "`test_cnv_caller` — evaluate a `Model[CNV]` against known calls."
);

// ---- Pipelines (orchestrate the above) -----------------------------------

deferred!(
    run_ngs_pipeline(gene: &str, output: &str, variants: Option<&str>, depth_of_coverage: Option<&str>, control_statistics: Option<&str>) -> Archive,
    "Beagle + samtools + sklearn (full NGS pipeline)",
    "`run_ngs_pipeline` — end-to-end NGS genotyping; orchestrates phasing, \
     `predict_alleles`, CNV calling, and result combination."
);
deferred!(
    run_chip_pipeline(gene: &str, output: &str, variants: &str) -> Archive,
    "Beagle (chip pipeline)",
    "`run_chip_pipeline` — end-to-end SNP-array genotyping."
);
deferred!(
    run_long_read_pipeline(gene: &str, output: &str, variants: &str) -> Archive,
    "long-read tooling",
    "`run_long_read_pipeline` — end-to-end long-read genotyping."
);

// ---- Plotting (matplotlib) -----------------------------------------------

deferred!(
    plot_bam_copy_number(copy_number: &Archive, path: &str) -> (),
    "matplotlib (use plotters when ported)",
    "`plot_bam_copy_number` — copy-number profile plot."
);
deferred!(
    plot_bam_read_depth(read_depth: &Archive, path: &str) -> (),
    "matplotlib",
    "`plot_bam_read_depth` — read-depth profile plot."
);
deferred!(
    plot_cn_af(copy_number: &Archive, imported_variants: &Archive, path: &str) -> (),
    "matplotlib",
    "`plot_cn_af` — copy-number vs allele-fraction plot."
);
deferred!(
    plot_vcf_allele_fraction(imported_variants: &Archive, path: &str) -> (),
    "matplotlib",
    "`plot_vcf_allele_fraction` — allele-fraction plot."
);
deferred!(
    plot_vcf_read_depth(imported_variants: &Archive, path: &str) -> (),
    "matplotlib",
    "`plot_vcf_read_depth` — VCF read-depth plot."
);
