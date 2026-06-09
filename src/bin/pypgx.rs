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
        /// Parallel worker threads (0 = auto / all cores).
        #[arg(long, default_value_t = 0)]
        jobs: usize,
    },
    /// Run the chip/array pipeline over target genes on a VCF (phasing against
    /// the bundle's 1KGP panel). Requires the `beagle` feature for unphased input.
    RunChipPipeline {
        #[arg(long)]
        vcf: String,
        #[arg(long, default_value = "GRCh38")]
        assembly: String,
        #[arg(long)]
        output: String,
        #[arg(long)]
        bundle: Option<String>,
        #[arg(long)]
        genes: Option<String>,
        /// Impute ungenotyped markers during phasing.
        #[arg(long)]
        impute: bool,
        /// Parallel worker threads (0 = auto / all cores).
        #[arg(long, default_value_t = 0)]
        jobs: usize,
    },
    /// Run the long-read pipeline over target genes on a VCF (read-backed
    /// phasing; no bundle/panel needed).
    RunLongReadPipeline {
        #[arg(long)]
        vcf: String,
        #[arg(long, default_value = "GRCh38")]
        assembly: String,
        #[arg(long)]
        output: String,
        #[arg(long)]
        genes: Option<String>,
        /// Parallel worker threads (0 = auto / all cores).
        #[arg(long, default_value_t = 0)]
        jobs: usize,
    },
    /// Subset an archive to (or, with --exclude, away from) the given samples.
    FilterSamples {
        /// Input archive (.zip).
        input: String,
        /// Comma-separated sample names.
        #[arg(long)]
        samples: String,
        /// Exclude the listed samples instead of keeping only them.
        #[arg(long)]
        exclude: bool,
        /// Output archive (.zip); prints a summary if omitted.
        #[arg(long)]
        output: Option<String>,
    },
    /// Plot a copy-number profile (PNG) from a CovFrame[CopyNumber] archive.
    #[cfg(feature = "plots")]
    PlotBamCopyNumber {
        input: String,
        #[arg(long)]
        output: String,
    },
    /// Plot a read-depth profile (PNG) from a CovFrame[ReadDepth] archive.
    #[cfg(feature = "plots")]
    PlotBamReadDepth {
        input: String,
        #[arg(long)]
        output: String,
    },
    /// Plot allele fraction (PNG) from a VcfFrame[Imported] archive.
    #[cfg(feature = "plots")]
    PlotVcfAlleleFraction {
        input: String,
        #[arg(long)]
        output: String,
    },
    /// Plot VCF read depth (PNG) from a VcfFrame[Imported] archive.
    #[cfg(feature = "plots")]
    PlotVcfReadDepth {
        input: String,
        #[arg(long)]
        output: String,
    },
    /// Plot copy-number vs allele-fraction (PNG) from CopyNumber + Imported.
    #[cfg(feature = "plots")]
    PlotCnAf {
        copy_number: String,
        imported: String,
        #[arg(long)]
        output: String,
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
    // Turn any panic — CLI glue or a deep library `.unwrap()` — into a clean
    // `error: ...` message and exit code 1, never a Rust backtrace.
    std::panic::set_hook(Box::new(|_| {}));
    if let Err(payload) =
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| dispatch(cli.command)))
    {
        let msg = payload
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "unknown error".to_string());
        // A closed downstream pipe (`pypgx ... | head`) is normal, not an error.
        if msg.contains("Broken pipe") {
            std::process::exit(0);
        }
        eprintln!("error: {msg}");
        std::process::exit(1);
    }
}

/// Dispatch one parsed subcommand. Any panic here is caught by `main` and
/// reported as a clean error, so library `.unwrap()`/`.expect()` and explicit
/// `panic!`s surface as `error: ...` rather than aborting with a backtrace.
fn dispatch(command: Command) {
    match command {
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
            jobs,
        } => {
            run_pipeline_cli(
                PipelineKind::Ngs,
                &vcf,
                &assembly,
                &output,
                bundle,
                genes,
                &platform,
                jobs,
            );
        }
        Command::RunChipPipeline {
            vcf,
            assembly,
            output,
            bundle,
            genes,
            impute,
            jobs,
        } => {
            run_pipeline_cli(
                PipelineKind::Chip { impute },
                &vcf,
                &assembly,
                &output,
                bundle,
                genes,
                "Chip",
                jobs,
            );
        }
        Command::RunLongReadPipeline {
            vcf,
            assembly,
            output,
            genes,
            jobs,
        } => {
            run_pipeline_cli(
                PipelineKind::LongRead,
                &vcf,
                &assembly,
                &output,
                None,
                genes,
                "LongRead",
                jobs,
            );
        }
        Command::FilterSamples {
            input,
            samples,
            exclude,
            output,
        } => {
            let archive = load(&input);
            let names: Vec<String> = samples.split(',').map(|s| s.trim().to_string()).collect();
            let result = pypgx::filter_samples(&archive, &names, exclude);
            save_or_print(&result, output);
        }
        #[cfg(feature = "plots")]
        Command::PlotBamCopyNumber { input, output } => {
            std::fs::create_dir_all(&output).expect("create output dir");
            pypgx::plot_bam_copy_number(&load(&input), Some(&output), None).expect("plot");
        }
        #[cfg(feature = "plots")]
        Command::PlotBamReadDepth { input, output } => {
            std::fs::create_dir_all(&output).expect("create output dir");
            pypgx::plot_bam_read_depth(&load(&input), Some(&output), None).expect("plot");
        }
        #[cfg(feature = "plots")]
        Command::PlotVcfAlleleFraction { input, output } => {
            std::fs::create_dir_all(&output).expect("create output dir");
            pypgx::plot_vcf_allele_fraction(&load(&input), Some(&output), None).expect("plot");
        }
        #[cfg(feature = "plots")]
        Command::PlotVcfReadDepth { input, output } => {
            std::fs::create_dir_all(&output).expect("create output dir");
            let archive = load(&input);
            let gene = archive.get("Gene").expect("Gene metadata");
            let assembly = archive.get("Assembly").expect("Assembly metadata");
            pypgx::plot_vcf_read_depth(gene, &archive.as_vcf(), assembly, Some(&output), None)
                .expect("plot");
        }
        #[cfg(feature = "plots")]
        Command::PlotCnAf {
            copy_number,
            imported,
            output,
        } => {
            std::fs::create_dir_all(&output).expect("create output dir");
            pypgx::plot_cn_af(&load(&copy_number), &load(&imported), Some(&output), None)
                .expect("plot");
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
        .map(|o| {
            String::from_utf8_lossy(&o.stdout)
                .lines()
                .next()
                .unwrap_or("")
                .starts_with("chr")
        })
        .unwrap_or(false)
}

/// Best-effort memory budget in bytes: the smaller of this process's cgroup
/// memory limit (v2 `memory.max`, else v1 `memory.limit_in_bytes`) and total
/// system RAM (`/proc/meminfo`). `None` if nothing usable is found (e.g. non-
/// Linux), in which case the caller leaves the job count uncapped.
fn memory_budget_bytes() -> Option<u64> {
    let parse = |s: &str| s.trim().parse::<u64>().ok();
    // cgroup v2: a number, or "max" for unlimited.
    let cgroup_v2 = std::fs::read_to_string("/sys/fs/cgroup/memory.max")
        .ok()
        .and_then(|s| parse(&s));
    // cgroup v1: a number; "unlimited" is a huge sentinel near u64::MAX.
    let cgroup_v1 = std::fs::read_to_string("/sys/fs/cgroup/memory/memory.limit_in_bytes")
        .ok()
        .and_then(|s| parse(&s))
        .filter(|&v| v < (1u64 << 62));
    let cgroup = cgroup_v2.or(cgroup_v1);
    // Total system RAM (MemTotal is in kB).
    let sys = std::fs::read_to_string("/proc/meminfo").ok().and_then(|s| {
        s.lines()
            .find(|l| l.starts_with("MemTotal:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|kb| kb.parse::<u64>().ok())
            .map(|kb| kb * 1024)
    });
    match (cgroup, sys) {
        (Some(c), Some(s)) => Some(c.min(s)),
        (c, s) => c.or(s),
    }
}

/// Which per-gene pipeline a `run-*-pipeline` subcommand drives.
#[derive(Clone, Copy)]
enum PipelineKind {
    Ngs,
    Chip { impute: bool },
    LongRead,
}

/// Shared driver for the `run-*-pipeline` subcommands: slice each gene region
/// from the VCF with `tabix`, run the chosen pipeline (NGS/chip phase against the
/// bundle's 1KGP panel; long-read needs none), and collect a per-gene
/// `results.tsv`. One bad gene never aborts the run — main's quiet panic hook +
/// the per-gene catch_unwind turn any failure into an `ERR:` row.
#[allow(clippy::too_many_arguments)]
fn run_pipeline_cli(
    kind: PipelineKind,
    vcf: &str,
    assembly: &str,
    output: &str,
    bundle: Option<String>,
    genes: Option<String>,
    platform: &str,
    jobs: usize,
) {
    // NGS/chip need the bundle (panel + default CNV model); long-read does not.
    let needs_bundle = !matches!(kind, PipelineKind::LongRead);
    let bundle = bundle.or_else(|| std::env::var("PYPGX_BUNDLE").ok());
    if needs_bundle {
        match &bundle {
            // Set before spawning workers (env is process-global) so predict_cnv's
            // default-model lookup and the panel path both resolve in every thread.
            Some(b) => std::env::set_var("PYPGX_BUNDLE", b),
            None => {
                eprintln!("error: no bundle path (pass --bundle or set $PYPGX_BUNDLE)");
                std::process::exit(2);
            }
        }
    }

    let mut gene_list: Vec<String> = match genes {
        Some(csv) => csv
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect(),
        None => pypgx::list_genes("target"),
    };

    std::fs::create_dir_all(output).expect("create output dir");
    let chr = vcf_has_chr_prefix(vcf);

    // LPT scheduling: phasing cost tracks panel size, so start the biggest-panel
    // genes (e.g. DPYD) first — they overlap with the many small genes and don't
    // become a serial tail at the end.
    if let Some(b) = bundle.as_deref() {
        gene_list.sort_by_key(|g| {
            std::cmp::Reverse(
                std::fs::metadata(format!("{b}/1kgp/{assembly}/{g}.vcf.gz"))
                    .map(|m| m.len())
                    .unwrap_or(0),
            )
        });
    }

    let n = gene_list.len();
    let mut jobs = match jobs {
        0 => std::thread::available_parallelism()
            .map(|j| j.get())
            .unwrap_or(1),
        j => j,
    }
    .clamp(1, n.max(1));

    // Memory guard: each concurrent gene runs a beagle-rs subprocess (~0.5 GB
    // peak on the biggest panel). In a memory-limited container (e.g. Docker
    // Desktop), all-cores phasing can OOM, so cap jobs to the memory budget —
    // the smaller of the cgroup limit and total RAM. Tunable via
    // $PYPGX_JOB_MEM_MB (per-job budget; default 1024 MB).
    let per_job_mb: u64 = std::env::var("PYPGX_JOB_MEM_MB")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&v| v > 0)
        .unwrap_or(1024);
    if let Some(budget) = memory_budget_bytes() {
        let mem_jobs = (budget / (per_job_mb * 1024 * 1024)).max(1) as usize;
        if mem_jobs < jobs {
            eprintln!(
                "pypgx-rs: capping jobs {jobs} -> {mem_jobs} to fit ~{} MB memory \
                 (~{per_job_mb} MB/job; set PYPGX_JOB_MEM_MB to override)",
                budget / 1024 / 1024
            );
            jobs = mem_jobs;
        }
    }

    // Each gene is independent (own tabix slice → pipeline → output dir), so fan
    // them out across `jobs` worker threads pulling from a shared index. beagle-rs
    // runs single-threaded per gene (nthreads=1), so N genes saturate N cores.
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    let genes = Arc::new(gene_list);
    let rows: Arc<Vec<Mutex<GeneCliResult>>> = Arc::new(
        (0..n)
            .map(|_| Mutex::new(GeneCliResult::default()))
            .collect(),
    );
    let next = Arc::new(AtomicUsize::new(0));

    let handles: Vec<_> = (0..jobs)
        .map(|_| {
            let (genes, rows, next) = (genes.clone(), rows.clone(), next.clone());
            let (vcf, assembly, output, platform) = (
                vcf.to_string(),
                assembly.to_string(),
                output.to_string(),
                platform.to_string(),
            );
            let bundle = bundle.clone();
            std::thread::spawn(move || loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= genes.len() {
                    break;
                }
                let row = process_gene(
                    kind,
                    &vcf,
                    &assembly,
                    &output,
                    bundle.as_deref(),
                    chr,
                    &platform,
                    &genes[i],
                );
                *rows[i].lock().unwrap() = row;
            })
        })
        .collect();
    for h in handles {
        let _ = h.join();
    }

    // Stable, execution-order-independent output: sort the rows by gene name.
    let mut out: Vec<(String, GeneCliResult)> = (0..n)
        .map(|i| {
            let row = std::mem::take(&mut *rows[i].lock().unwrap());
            (genes[i].clone(), row)
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));

    let mut summary = format!("Gene\tStatus\t{}\n", RESULT_COLUMNS.join("\t"));
    let mut details = String::from(
        "Gene\tStatus\tSample\tGenotype\tPhenotype\tHaplotype\tCandidateRank\tAllele\tVariant\tAlleleFraction\tAlternativePhase\tCNV\n",
    );
    let (mut ok, mut failed) = (0usize, 0usize);
    for (_, result) in &out {
        summary.push_str(&result.summary);
        details.push_str(&result.details);
        if result.ok {
            ok += 1;
        } else {
            failed += 1;
        }
    }
    let summary_path = format!("{output}/results.tsv");
    std::fs::write(&summary_path, &summary).expect("write results.tsv");
    let details_path = format!("{output}/results.long.tsv");
    std::fs::write(&details_path, &details).expect("write results.long.tsv");
    eprintln!(
        "pypgx-rs: {ok} ok, {failed} failed over {n} genes ({jobs} jobs) -> {summary_path}, {details_path}"
    );
}

#[derive(Default)]
struct GeneCliResult {
    summary: String,
    details: String,
    ok: bool,
}

/// Run one gene's pipeline (tabix slice → chosen pipeline → read back the call)
/// and return its `results.tsv` row plus whether it succeeded. Pure w.r.t. shared
/// state: writes only to `{output}/{gene}/`, so it is safe to call concurrently.
#[allow(clippy::too_many_arguments)]
fn process_gene(
    kind: PipelineKind,
    vcf: &str,
    assembly: &str,
    output: &str,
    bundle: Option<&str>,
    chr: bool,
    platform: &str,
    gene: &str,
) -> GeneCliResult {
    let region = match pypgx::core::get_region(gene, assembly) {
        Ok(r) => {
            if chr {
                format!("chr{r}")
            } else {
                r
            }
        }
        Err(e) => return error_result(gene, &format!("ERR:region:{e}")),
    };
    let sliced = match ShellCommand::new("tabix")
        .arg("-h")
        .arg(vcf)
        .arg(&region)
        .output()
    {
        Ok(o) if o.status.success() => o.stdout,
        _ => return error_result(gene, "ERR:tabix"),
    };
    let vf = pypgx::fuc::VcfFrame::from_string(&String::from_utf8_lossy(&sliced));
    let panel = bundle.map(|b| format!("{b}/1kgp/{assembly}/{gene}.vcf.gz"));
    let panel_opt = panel
        .as_deref()
        .filter(|p| std::path::Path::new(p).exists());
    let geneout = format!("{output}/{gene}");

    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| match kind {
        PipelineKind::Ngs => pypgx::run_ngs_pipeline(
            gene,
            &geneout,
            Some(&vf),
            None,
            None,
            platform,
            assembly,
            panel_opt,
            true,
            None,
            false,
            None,
            None,
        ),
        PipelineKind::Chip { impute } => pypgx::run_chip_pipeline(
            gene, &geneout, &vf, assembly, panel_opt, impute, true, None, false,
        ),
        PipelineKind::LongRead => {
            pypgx::run_long_read_pipeline(gene, &geneout, &vf, assembly, true, None, false)
        }
    }));
    match run {
        Ok(Ok(())) => {
            let mut values = vec![String::new(); RESULT_COLUMNS.len()];
            let mut details = String::new();
            if let Ok(r) = pypgx::Archive::from_file(&format!("{geneout}/results.zip")) {
                let st = r.as_sample_table();
                if let Some(row) = st.rows.first() {
                    for (i, col) in RESULT_COLUMNS.iter().enumerate() {
                        if let Some(ci) = st.columns.iter().position(|c| c == col) {
                            values[i] = row[ci].clone();
                        }
                    }
                }
                details = long_result_rows(gene, "ok", &r);
            }
            GeneCliResult {
                summary: format!(
                    "{gene}\tok\t{}\n",
                    values
                        .iter()
                        .map(|v| tsv_cell(v))
                        .collect::<Vec<_>>()
                        .join("\t")
                ),
                details,
                ok: true,
            }
        }
        Ok(Err(e)) => {
            let msg: String = e
                .to_string()
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ");
            error_result(
                gene,
                &format!("ERR:{}", msg.chars().take(80).collect::<String>()),
            )
        }
        Err(_) => error_result(gene, "ERR:panic"),
    }
}

const RESULT_COLUMNS: [&str; 7] = [
    "Genotype",
    "Phenotype",
    "Haplotype1",
    "Haplotype2",
    "AlternativePhase",
    "VariantData",
    "CNV",
];

fn error_result(gene: &str, status: &str) -> GeneCliResult {
    GeneCliResult {
        summary: format!(
            "{gene}\t{}\t{}\n",
            tsv_cell(status),
            vec![String::new(); RESULT_COLUMNS.len()].join("\t")
        ),
        details: format!("{gene}\t{}\t\t\t\t\t\t\t\t\t\t\n", tsv_cell(status)),
        ok: false,
    }
}

fn long_result_rows(gene: &str, status: &str, archive: &pypgx::Archive) -> String {
    let st = archive.as_sample_table();
    let col = |name: &str| st.columns.iter().position(|c| c == name);
    let genotype_c = col("Genotype");
    let phenotype_c = col("Phenotype");
    let haplotype1_c = col("Haplotype1");
    let haplotype2_c = col("Haplotype2");
    let alt_phase_c = col("AlternativePhase");
    let variant_data_c = col("VariantData");
    let cnv_c = col("CNV");

    let mut out = String::new();
    for (ri, sample) in st.index.iter().enumerate() {
        let row = &st.rows[ri];
        let genotype = genotype_c.map(|i| row[i].as_str()).unwrap_or("");
        let phenotype = phenotype_c.map(|i| row[i].as_str()).unwrap_or("");
        let alt_phase = alt_phase_c.map(|i| row[i].as_str()).unwrap_or("");
        let cnv = cnv_c.map(|i| row[i].as_str()).unwrap_or("");
        let variant_data = variant_data_c.map(|i| row[i].as_str()).unwrap_or("");
        let variants = parse_variant_data_entries(variant_data);

        for (label, ci) in [
            ("Haplotype1", haplotype1_c),
            ("Haplotype2", haplotype2_c),
            ("AlternativePhase", alt_phase_c),
        ] {
            let Some(ci) = ci else { continue };
            for (rank, allele) in split_semis(&row[ci]).iter().enumerate() {
                let entries = variants
                    .iter()
                    .find(|(a, _)| a == allele)
                    .map(|(_, entries)| entries.as_slice())
                    .unwrap_or(&[]);
                if entries.is_empty() {
                    push_long_row(
                        &mut out,
                        gene,
                        status,
                        sample,
                        genotype,
                        phenotype,
                        label,
                        rank + 1,
                        allele,
                        "",
                        "",
                        alt_phase,
                        cnv,
                    );
                } else {
                    for (variant, fraction) in entries {
                        push_long_row(
                            &mut out,
                            gene,
                            status,
                            sample,
                            genotype,
                            phenotype,
                            label,
                            rank + 1,
                            allele,
                            variant,
                            fraction,
                            alt_phase,
                            cnv,
                        );
                    }
                }
            }
        }
    }
    out
}

fn split_semis(value: &str) -> Vec<String> {
    value
        .trim_matches(';')
        .split(';')
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn parse_variant_data_entries(value: &str) -> Vec<(String, Vec<(String, String)>)> {
    let mut out = Vec::new();
    for entry in split_semis(value) {
        let fields: Vec<&str> = entry.splitn(3, ':').collect();
        if fields.len() < 2 {
            continue;
        }
        let allele = fields[0].to_string();
        let entries = if fields[1] == "default" {
            vec![("default".to_string(), String::new())]
        } else if fields.len() == 3 {
            let variants: Vec<&str> = fields[1].split(',').collect();
            let fractions: Vec<&str> = fields[2].split(',').collect();
            variants
                .iter()
                .enumerate()
                .map(|(i, variant)| {
                    (
                        (*variant).to_string(),
                        fractions.get(i).copied().unwrap_or("").to_string(),
                    )
                })
                .collect()
        } else {
            vec![(fields[1].to_string(), String::new())]
        };
        out.push((allele, entries));
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn push_long_row(
    out: &mut String,
    gene: &str,
    status: &str,
    sample: &str,
    genotype: &str,
    phenotype: &str,
    haplotype: &str,
    rank: usize,
    allele: &str,
    variant: &str,
    fraction: &str,
    alt_phase: &str,
    cnv: &str,
) {
    out.push_str(
        &[
            tsv_cell(gene),
            tsv_cell(status),
            tsv_cell(sample),
            tsv_cell(genotype),
            tsv_cell(phenotype),
            tsv_cell(haplotype),
            rank.to_string(),
            tsv_cell(allele),
            tsv_cell(variant),
            tsv_cell(fraction),
            tsv_cell(alt_phase),
            tsv_cell(cnv),
        ]
        .join("\t"),
    );
    out.push('\n');
}

fn tsv_cell(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}
