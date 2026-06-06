//! pypgx-rs — a 1-for-1 Rust port of [PyPGx](https://github.com/sbslee/pypgx),
//! a package for pharmacogenomics (PGx) research.
//!
//! Module layout mirrors the Python package:
//! - [`core`]  ← `pypgx.api.core`   (reference tables + allele/phenotype logic)
//! - [`api`]   ← `pypgx.api.utils`  (`predict_alleles`, ...)
//! - [`sdk`]   ← `pypgx.sdk`         (`Archive`, semantic types, errors)
//! - [`fuc`]   ← the slice of `fuc` PyPGx depends on (`VcfFrame`, variants)
//! - [`table`] ← a pandas-like frame used to load the CSV tables

pub mod api;
pub mod bed;
pub mod cnv;
pub mod core;
pub mod external;
pub mod fuc;
pub mod genotype;
pub mod pipeline;
#[cfg(feature = "plots")]
pub mod plot;
pub mod sdk;
pub mod table;

// Re-export the public surface, mirroring `pypgx/__init__.py`.
pub use crate::api::{
    call_phenotypes, combine_results, compare_genotypes, compute_copy_number, count_alleles,
    create_consolidated_vcf, create_regions_bed, filter_samples, import_read_depth, import_variants,
    predict_alleles, predict_cnv, test_cnv_caller,
};
pub use crate::cnv::CnvModel;
pub use crate::bed::BedFrame;
pub use crate::core::{
    build_definition_table, collapse_alleles, get_default_allele, get_exon_ends, get_exon_starts,
    get_function, get_paralog, get_priority, get_recommendation, get_ref_allele, get_region,
    get_score, get_strand, get_variant_impact, get_variant_synonyms, has_phenotype, has_score,
    has_sv, is_legit_allele, is_target_gene, list_alleles, list_functions, list_genes,
    list_phenotypes, list_variants, load_allele_table, load_cnv_table, load_cpic_table,
    load_diplotype_table, load_equation_table, load_gene_table, load_phenotype_table,
    load_recommendation_table, load_variant_table, predict_phenotype, predict_score, sort_alleles,
};
pub use crate::genotype::call_genotypes;
pub use crate::pipeline::{run_chip_pipeline, run_long_read_pipeline, run_ngs_pipeline};
pub use crate::sdk::{Archive, PgxError};
#[cfg(feature = "plots")]
pub use crate::plot::{
    plot_bam_copy_number, plot_bam_read_depth, plot_cn_af, plot_vcf_allele_fraction,
    plot_vcf_read_depth,
};
#[cfg(feature = "bam")]
pub use crate::api::{compute_control_statistics, compute_target_depth, prepare_depth_of_coverage};
