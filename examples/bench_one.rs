//! Run pypgx-rs's NGS pipeline for one gene on a (pre-sliced, plain-text) VCF
//! and print `gene\tseconds\tstatus\tgenotype\tphenotype`. Build with the
//! `beagle` feature so the phasing step is live.
//!
//! Args: <gene> <sliced_vcf> <assembly> <out_dir> [panel.vcf.gz]
//! (pypgx-rs has no file/region reader, so the caller pre-slices with tabix.)
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (gene, vcf_path, assembly, out) = (&a[1], &a[2], &a[3], &a[4]);
    let panel = a.get(5).filter(|p| !p.is_empty()).map(|s| s.as_str());

    let text = std::fs::read_to_string(vcf_path).expect("read sliced vcf");
    let vf = pypgx::fuc::VcfFrame::from_string(&text);

    let t = Instant::now();
    let res = pypgx::run_ngs_pipeline(
        gene, out, Some(&vf), None, None, "WGS", assembly, panel, true, None, false, None, None,
    );
    let dt = t.elapsed().as_secs_f64();

    match res {
        Ok(()) => {
            let r = pypgx::Archive::from_file(&format!("{out}/results.zip")).expect("results");
            let st = r.as_sample_table();
            let gi = st.columns.iter().position(|c| c == "Genotype").unwrap();
            let pi = st.columns.iter().position(|c| c == "Phenotype").unwrap();
            let g = st.rows.first().map(|row| row[gi].clone()).unwrap_or_default();
            let p = st.rows.first().map(|row| row[pi].clone()).unwrap_or_default();
            println!("{gene}\t{dt:.3}\tok\t{g}\t{p}");
        }
        Err(e) => {
            let msg = e.to_string().replace('\t', " ");
            let msg: String = msg.chars().take(80).collect();
            println!("{gene}\t{dt:.3}\tERR:{msg}\t\t");
        }
    }
}
