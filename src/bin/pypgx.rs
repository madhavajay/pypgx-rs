//! PyPGx CLI (port of `pypgx/__main__.py` + `pypgx/cli/*`).
//!
//! Implements the data-inspection and pure analytical commands that need no
//! external tools: the predict→genotype→phenotype→combine→report pipeline.
//! External-tool / CNV / plot subcommands are deferred (see TODO.md and the
//! `pypgx::external` module).

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "pypgx",
    version,
    about = "Pharmacogenomics (PGx) toolkit (Rust port)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// List genes in the gene table.
    ListGenes {
        /// Gene set: target, control, or all.
        #[arg(long, default_value = "target")]
        mode: String,
    },
    /// List star alleles for a target gene.
    ListAlleles {
        gene: String,
        #[arg(long, default_value = "GRCh37")]
        assembly: String,
    },
    /// List variants that define star alleles for a gene.
    ListVariants {
        gene: String,
        #[arg(long, default_value = "all")]
        mode: String,
        #[arg(long, default_value = "GRCh37")]
        assembly: String,
    },
    /// Predict candidate star alleles from a VcfFrame[Consolidated] archive.
    PredictAlleles {
        /// Input archive (.zip) with semantic type VcfFrame[Consolidated].
        input: String,
        /// Optional output archive (.zip) with semantic type SampleTable[Alleles].
        #[arg(long)]
        output: Option<String>,
    },
    /// Call diplotypes from an alleles (and optional CNV-calls) archive.
    CallGenotypes {
        /// SampleTable[Alleles] archive (.zip).
        #[arg(long)]
        alleles: Option<String>,
        /// SampleTable[CNVCalls] archive (.zip).
        #[arg(long)]
        cnv_calls: Option<String>,
        /// Output SampleTable[Genotypes] archive (.zip).
        #[arg(long)]
        output: Option<String>,
    },
    /// Call phenotypes from a genotypes archive.
    CallPhenotypes {
        /// SampleTable[Genotypes] archive (.zip).
        genotypes: String,
        #[arg(long)]
        output: Option<String>,
    },
    /// Combine genotypes/phenotypes/alleles/CNV into a results archive.
    CombineResults {
        #[arg(long)]
        genotypes: Option<String>,
        #[arg(long)]
        phenotypes: Option<String>,
        #[arg(long)]
        alleles: Option<String>,
        #[arg(long)]
        cnv_calls: Option<String>,
        #[arg(long)]
        output: Option<String>,
    },
    /// Compare two results archives (concordance report).
    CompareGenotypes {
        first: String,
        second: String,
        #[arg(long)]
        verbose: bool,
    },
    /// Count star alleles in a results archive.
    CountAlleles {
        /// SampleTable[Results] archive (.zip).
        results: String,
    },
    /// Print the metadata of an archive.
    PrintMetadata {
        /// Input archive (.zip).
        input: String,
    },
    /// Print the main data table of a SampleTable archive (TSV).
    PrintData {
        /// Input archive (.zip).
        input: String,
    },
}

fn load(path: &str) -> pypgx::Archive {
    pypgx::Archive::from_file(path).expect("read archive")
}

/// Print a SampleTable as TSV (empty header cell for the index column).
fn print_sample_table(archive: &pypgx::Archive) {
    let t = archive.as_sample_table();
    println!("\t{}", t.columns.join("\t"));
    for (i, idx) in t.index.iter().enumerate() {
        println!("{idx}\t{}", t.rows[i].join("\t"));
    }
}

fn save_or_print(archive: &pypgx::Archive, output: Option<String>) {
    match output {
        Some(path) => {
            archive.to_file(&path).expect("write archive");
            eprintln!("Saved {} to: {}", archive.semantic_type(), path);
        }
        None => print_sample_table(archive),
    }
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Command::ListGenes { mode } => {
            for g in pypgx::list_genes(&mode) {
                println!("{g}");
            }
        }
        Command::ListAlleles { gene, assembly } => {
            for a in pypgx::list_alleles(&gene, None, &assembly) {
                println!("{a}");
            }
        }
        Command::ListVariants {
            gene,
            mode,
            assembly,
        } => {
            for v in pypgx::list_variants(&gene, None, &mode, &assembly) {
                println!("{v}");
            }
        }
        Command::PredictAlleles { input, output } => {
            let archive = load(&input);
            let result = pypgx::predict_alleles(&archive).expect("predict_alleles");
            save_or_print(&result, output);
        }
        Command::CallGenotypes {
            alleles,
            cnv_calls,
            output,
        } => {
            let a = alleles.map(|p| load(&p));
            let c = cnv_calls.map(|p| load(&p));
            let result = pypgx::call_genotypes(a.as_ref(), c.as_ref()).expect("call_genotypes");
            save_or_print(&result, output);
        }
        Command::CallPhenotypes { genotypes, output } => {
            let gt = load(&genotypes);
            let result = pypgx::call_phenotypes(&gt).expect("call_phenotypes");
            save_or_print(&result, output);
        }
        Command::CombineResults {
            genotypes,
            phenotypes,
            alleles,
            cnv_calls,
            output,
        } => {
            let gt = genotypes.map(|p| load(&p));
            let ph = phenotypes.map(|p| load(&p));
            let al = alleles.map(|p| load(&p));
            let cn = cnv_calls.map(|p| load(&p));
            let result = pypgx::combine_results(gt.as_ref(), ph.as_ref(), al.as_ref(), cn.as_ref())
                .expect("combine_results");
            save_or_print(&result, output);
        }
        Command::CompareGenotypes {
            first,
            second,
            verbose,
        } => {
            let a = load(&first);
            let b = load(&second);
            print!("{}", pypgx::compare_genotypes(&a, &b, verbose));
        }
        Command::CountAlleles { results } => {
            let r = load(&results);
            for (allele, count) in pypgx::count_alleles(&r) {
                println!("{allele}\t{count}");
            }
        }
        Command::PrintMetadata { input } => {
            let archive = load(&input);
            for (k, v) in &archive.metadata {
                println!("{k}={v}");
            }
        }
        Command::PrintData { input } => {
            print_sample_table(&load(&input));
        }
    }
}
