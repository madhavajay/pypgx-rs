//! PyPGx CLI (port of `pypgx/__main__.py` + `pypgx/cli/*`).
//!
//! Implements the data-inspection and pure analytical commands that need no
//! external tools: the predict→genotype→phenotype→combine→report pipeline.
//! External-tool / CNV / plot subcommands are deferred (see TODO.md and the
//! `pypgx::external` module).

use clap::{Parser, Subcommand};
use std::process::Command as ShellCommand;

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
    /// Run the end-to-end NGS pipeline over target genes on a VCF (phasing via
    /// beagle-rs + 1KGP panels from the bundle), writing per-gene archives and a
    /// combined results.tsv. Requires the `beagle` feature for unphased input.
    RunNgsPipeline {
        /// Input VCF, bgzipped + tabix-indexed (e.g. sample.vcf.gz).
        #[arg(long)]
        vcf: String,
        /// Genome assembly: GRCh37 or GRCh38.
        #[arg(long, default_value = "GRCh38")]
        assembly: String,
        /// Output directory (one subdir per gene + a results.tsv summary).
        #[arg(long)]
        output: String,
        /// pypgx-bundle path (1KGP panels + CNV models). Defaults to $PYPGX_BUNDLE.
        #[arg(long)]
        bundle: Option<String>,
        /// Comma-separated gene list; default = all target genes.
        #[arg(long)]
        genes: Option<String>,
        /// Sequencing platform metadata.
        #[arg(long, default_value = "WGS")]
        platform: String,
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
        Command::RunNgsPipeline {
            vcf,
            assembly,
            output,
            bundle,
            genes,
            platform,
        } => {
            run_ngs_cli(&vcf, &assembly, &output, bundle, genes, &platform);
        }
    }
}

/// Whether the VCF's contigs are `chr`-prefixed (peek `tabix -l`).
fn vcf_has_chr_prefix(vcf: &str) -> bool {
    ShellCommand::new("tabix")
        .arg("-l")
        .arg(vcf)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).lines().next().unwrap_or("").starts_with("chr"))
        .unwrap_or(false)
}

/// `pypgx run-ngs-pipeline` — slice each gene region from the VCF with `tabix`,
/// run the NGS pipeline (phasing against the bundle's 1KGP panel), and collect a
/// per-gene `results.tsv`. One bad gene never aborts the run.
fn run_ngs_cli(
    vcf: &str,
    assembly: &str,
    output: &str,
    bundle: Option<String>,
    genes: Option<String>,
    platform: &str,
) {
    // Resolve the bundle and export it so predict_cnv's default-model lookup and
    // panel resolution both find it.
    let bundle = bundle
        .or_else(|| std::env::var("PYPGX_BUNDLE").ok())
        .unwrap_or_else(|| {
            eprintln!("error: no bundle path (pass --bundle or set $PYPGX_BUNDLE)");
            std::process::exit(2);
        });
    std::env::set_var("PYPGX_BUNDLE", &bundle);

    let gene_list: Vec<String> = match genes {
        Some(csv) => csv.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect(),
        None => pypgx::list_genes("target"),
    };

    std::fs::create_dir_all(output).expect("create output dir");
    let tmp = std::env::temp_dir().join(format!("pypgx_slices_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).expect("create slice dir");
    let chr = vcf_has_chr_prefix(vcf);
    // Internal panics become per-gene errors; keep the console quiet.
    std::panic::set_hook(Box::new(|_| {}));

    let mut summary = String::from("Gene\tStatus\tGenotype\tPhenotype\n");
    let (mut ok, mut failed) = (0usize, 0usize);
    for gene in &gene_list {
        let region = match pypgx::core::get_region(gene, assembly) {
            Ok(r) => if chr { format!("chr{r}") } else { r },
            Err(e) => {
                summary.push_str(&format!("{gene}\tERR:region:{e}\t\t\n"));
                failed += 1;
                continue;
            }
        };
        let sliced = match ShellCommand::new("tabix").arg("-h").arg(vcf).arg(&region).output() {
            Ok(o) if o.status.success() => o.stdout,
            _ => {
                summary.push_str(&format!("{gene}\tERR:tabix\t\t\n"));
                failed += 1;
                continue;
            }
        };
        let vf = pypgx::fuc::VcfFrame::from_string(&String::from_utf8_lossy(&sliced));
        let panel = format!("{bundle}/1kgp/{assembly}/{gene}.vcf.gz");
        let panel_opt = std::path::Path::new(&panel).exists().then_some(panel.as_str());
        let geneout = format!("{output}/{gene}");

        let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pypgx::run_ngs_pipeline(
                gene, &geneout, Some(&vf), None, None, platform, assembly, panel_opt, true, None,
                false, None, None,
            )
        }));
        match run {
            Ok(Ok(())) => {
                let (mut g, mut p) = (String::new(), String::new());
                if let Ok(r) = pypgx::Archive::from_file(&format!("{geneout}/results.zip")) {
                    let st = r.as_sample_table();
                    if let Some(gi) = st.columns.iter().position(|c| c == "Genotype") {
                        g = st.rows.first().map(|r| r[gi].clone()).unwrap_or_default();
                    }
                    if let Some(pi) = st.columns.iter().position(|c| c == "Phenotype") {
                        p = st.rows.first().map(|r| r[pi].clone()).unwrap_or_default();
                    }
                }
                summary.push_str(&format!("{gene}\tok\t{g}\t{p}\n"));
                ok += 1;
            }
            Ok(Err(e)) => {
                let msg: String = e.to_string().split_whitespace().collect::<Vec<_>>().join(" ");
                summary.push_str(&format!("{gene}\tERR:{}\t\t\n", msg.chars().take(80).collect::<String>()));
                failed += 1;
            }
            Err(_) => {
                summary.push_str(&format!("{gene}\tERR:panic\t\t\n"));
                failed += 1;
            }
        }
    }
    let _ = std::panic::take_hook();
    let summary_path = format!("{output}/results.tsv");
    std::fs::write(&summary_path, &summary).expect("write results.tsv");
    eprintln!("pypgx-rs: {ok} ok, {failed} failed over {} genes -> {summary_path}", gene_list.len());
}
