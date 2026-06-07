//! Run pypgx-rs's NGS pipeline over every target gene on a whole-genome VCF and
//! time it. pypgx-rs has no file/region reader, so we shell out to `tabix` to
//! slice each gene region (chr-prefixed) — that slice time is included in the
//! per-gene total, for a fair comparison with PyPGx (which slices internally).
//!
//! Args: <vcf.gz> <out_dir> <assembly> <bundle_dir>
//! Env:  BEAGLE_RS_BIN=<beagle-rs binary>. Build with `--features beagle`.
use std::process::Command;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (vcf, out_dir, assembly, bundle) = (&a[1], &a[2], &a[3], &a[4]);
    let genes = pypgx::list_genes("target");
    // pypgx-rs panics (rather than returning typed errors) on some real-data
    // edge cases; isolate per gene so one panic can't abort the whole run.
    std::panic::set_hook(Box::new(|_| {}));
    std::fs::create_dir_all(out_dir).unwrap();
    let tmp = std::env::temp_dir().join("pypgxrs_slices");
    std::fs::create_dir_all(&tmp).unwrap();

    let mut rows: Vec<(String, f64, String, String, String)> = Vec::new();
    let total = Instant::now();
    for gene in &genes {
        let t = Instant::now();
        let (mut status, mut g, mut p) = ("ok".to_string(), String::new(), String::new());

        let region = match pypgx::core::get_region(gene, assembly) {
            Ok(r) => format!("chr{r}"), // input VCF is chr-prefixed (GRCh38)
            Err(e) => {
                push(&mut rows, gene, t, format!("ERR:region:{e}"), g, p);
                continue;
            }
        };
        // Slice the gene region with tabix (random access via the .tbi).
        let sliced = match Command::new("tabix").arg("-h").arg(vcf).arg(&region).output() {
            Ok(o) if o.status.success() => o.stdout,
            _ => {
                push(&mut rows, gene, t, "ERR:tabix".into(), g, p);
                continue;
            }
        };
        let slice_path = tmp.join(format!("{gene}.vcf"));
        std::fs::write(&slice_path, &sliced).unwrap();
        let vf = pypgx::fuc::VcfFrame::from_string(&String::from_utf8_lossy(&sliced));

        let panel = format!("{bundle}/1kgp/{assembly}/{gene}.vcf.gz");
        let panel_opt = std::path::Path::new(&panel).exists().then_some(panel.as_str());
        let geneout = format!("{out_dir}/{gene}");

        let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pypgx::run_ngs_pipeline(
                gene, &geneout, Some(&vf), None, None, "WGS", assembly, panel_opt, true, None,
                false, None, None,
            )
        }));
        match run {
            Ok(Ok(())) => {
                if let Ok(r) = pypgx::Archive::from_file(&format!("{geneout}/results.zip")) {
                    let st = r.as_sample_table();
                    let gi = st.columns.iter().position(|c| c == "Genotype").unwrap();
                    let pi = st.columns.iter().position(|c| c == "Phenotype").unwrap();
                    g = st.rows.first().map(|r| r[gi].clone()).unwrap_or_default();
                    p = st.rows.first().map(|r| r[pi].clone()).unwrap_or_default();
                }
            }
            Ok(Err(e)) => {
                // Collapse all whitespace (subprocess stderr can be multi-line).
                let msg: String = e.to_string().split_whitespace().collect::<Vec<_>>().join(" ");
                status = format!("ERR:{}", msg.chars().take(70).collect::<String>());
            }
            Err(_) => status = "ERR:panic".to_string(),
        }
        let dt = t.elapsed().as_secs_f64();
        println!("{gene}\t{dt:.3}\t{status}\t{g}\t{p}");
        rows.push((gene.clone(), dt, status, g, p));
    }
    let tot = total.elapsed().as_secs_f64();

    let mut s = String::from("Gene\tSeconds\tStatus\tGenotype\tPhenotype\n");
    for (gene, dt, st, g, p) in &rows {
        s.push_str(&format!("{gene}\t{dt:.3}\t{st}\t{g}\t{p}\n"));
    }
    s.push_str(&format!("#TOTAL\t{tot:.3}\t{} genes\n", genes.len()));
    std::fs::write(format!("{out_dir}/timing.tsv"), s).unwrap();
    let ok = rows.iter().filter(|r| r.2 == "ok").count();
    eprintln!("pypgx-rs: {ok}/{} genes ok in {tot:.1}s total", genes.len());
}

fn push(
    rows: &mut Vec<(String, f64, String, String, String)>,
    gene: &str,
    t: Instant,
    status: String,
    g: String,
    p: String,
) {
    let dt = t.elapsed().as_secs_f64();
    println!("{gene}\t{dt:.3}\t{status}\t{g}\t{p}");
    rows.push((gene.to_string(), dt, status, g, p));
}
