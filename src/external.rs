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

/// `estimate_phase_beagle` — statistical phasing of a `VcfFrame[Imported]`
/// archive into `VcfFrame[Phased]` by shelling out to the **beagle-rs** binary
/// (a drop-in for the Beagle jar): `beagle-rs gt=input.vcf chrom=<region>
/// [ref=<panel>] out=<prefix> impute=<bool> em=<bool>`, with PyPGx's EM-skip
/// retry (`em=true` then `em=false`). The bgzf `output.vcf.gz` is read back into
/// a `VcfFrame[Phased]`. Binary discovery: `$BEAGLE_RS_BIN`, else `beagle-rs` on
/// PATH. We deliberately do not pass `nthreads`, matching PyPGx's Java Beagle
/// invocation; Beagle phase orientation can be thread-count-sensitive.
///
/// ⚠️ **Not yet byte-parity with PyPGx.** PyPGx bundles Beagle **22Jul22.46e**;
/// beagle-rs targets **27Feb25.75f** — phasing output can differ across Beagle
/// versions. And `panel=None` here means *pure phasing with no reference panel*;
/// PyPGx instead loads the 1KGP panel from `pypgx-bundle` (absent here). Treat
/// this as a working integration pending version + panel reconciliation.
#[cfg(feature = "beagle")]
pub fn estimate_phase_beagle(
    imported_variants: &Archive,
    panel: Option<&str>,
    impute: bool,
) -> Result<Archive, PgxError> {
    use crate::sdk::ArchiveData;
    use std::io::Read;

    let io = |e: std::io::Error| PgxError::External(e.to_string());
    imported_variants.check_type(&["VcfFrame[Imported]"])?;
    let gene = imported_variants
        .get("Gene")
        .ok_or_else(|| PgxError::External("missing Gene metadata".into()))?;
    let assembly = imported_variants
        .get("Assembly")
        .ok_or_else(|| PgxError::External("missing Assembly metadata".into()))?;
    let region = crate::core::get_region(gene, assembly)?;

    let mut metadata = imported_variants.copy_metadata();
    if let Some(e) = metadata.iter_mut().find(|(k, _)| k == "SemanticType") {
        e.1 = "VcfFrame[Phased]".to_string();
    }
    metadata.push(("Program".to_string(), "Beagle".to_string()));

    let vf = imported_variants.as_vcf();
    if vf.rows.is_empty() {
        return Ok(Archive::new(metadata, ArchiveData::Vcf(vf.clone())));
    }

    // Panel-driven chr-prefix (mirrors PyPGx): if the reference panel is
    // `chr`-prefixed, add `chr` to the (de-chr'd) input + region before phasing,
    // then strip it back off the output.
    let add_chr = match panel {
        Some(p) => panel_has_chr_prefix(p).unwrap_or(false),
        None => false,
    };
    let vf_run = if add_chr {
        std::borrow::Cow::Owned(vf.update_chr_prefix("add"))
    } else {
        std::borrow::Cow::Borrowed(vf)
    };
    let region_run = if add_chr && !region.starts_with("chr") {
        format!("chr{region}")
    } else {
        region.clone()
    };

    // Overlap pre-check (mirrors PyPGx): Beagle errors when 0 or 1 input markers
    // overlap the reference panel, so handle those cases without running it.
    if let Some(p) = panel {
        let panel_vars: std::collections::HashSet<String> =
            panel_variants(p).map_err(io)?.into_iter().collect();
        let common: std::collections::HashSet<String> = vf_run
            .to_variants()
            .into_iter()
            .filter(|v| panel_vars.contains(v))
            .collect();
        if common.len() <= 1 {
            let phased = if let Some(variant) = common.iter().next() {
                // One overlapping marker: pseudo-phase just that variant's row.
                let pos = variant.split('-').nth(1).unwrap_or("");
                let rows: Vec<Vec<String>> = vf_run
                    .rows
                    .iter()
                    .filter(|r| r[1] == pos)
                    .cloned()
                    .collect();
                crate::fuc::VcfFrame::new(Vec::new(), vf_run.columns.clone(), rows)
                    .pseudophase()
                    .strip("GT")
            } else {
                // No overlapping markers: skip statistical phasing entirely.
                crate::fuc::VcfFrame::new(Vec::new(), vf_run.columns.clone(), Vec::new())
            };
            let phased = if add_chr {
                phased.update_chr_prefix("remove")
            } else {
                phased
            };
            return Ok(Archive::new(metadata, ArchiveData::Vcf(phased)));
        }
    }

    // Unique per call (pid + atomic seq) so concurrent phasing of multiple genes
    // in the same process never collides on the temp input/output files.
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let seq = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("pypgx_beagle_{}_{seq}", std::process::id()));
    std::fs::create_dir_all(&dir).map_err(io)?;
    let input = dir.join("input.vcf");
    let out_prefix = dir.join("output");

    // Imported frames carry no meta, so prepend a ##fileformat header.
    let mut vcf = String::from("##fileformat=VCFv4.2\n");
    for m in &vf_run.meta {
        vcf.push_str(m);
        vcf.push('\n');
    }
    vcf.push('#');
    vcf.push_str(&vf_run.columns.join("\t"));
    vcf.push('\n');
    for r in &vf_run.rows {
        vcf.push_str(&r.join("\t"));
        vcf.push('\n');
    }
    std::fs::write(&input, vcf).map_err(io)?;

    let bin = std::env::var("BEAGLE_RS_BIN").unwrap_or_else(|_| "beagle-rs".to_string());
    let run = |em: bool| -> std::io::Result<std::process::Output> {
        let mut cmd = beagle_command(&bin, &input, &region_run, &out_prefix, impute, em, panel);
        cmd.output()
    };

    // EM-skip retry, mirroring PyPGx (em=true, then em=false on failure).
    let mut output = run(true).map_err(io)?;
    if !output.status.success() {
        output = run(false).map_err(io)?;
    }
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        let _ = std::fs::remove_dir_all(&dir);
        // Markers too distant land in separate windows (Beagle: "Window has
        // only one position"); PyPGx skips phasing in that case rather than crash.
        if stderr.contains("Window has only one position") {
            let empty = crate::fuc::VcfFrame::new(Vec::new(), vf_run.columns.clone(), Vec::new());
            let empty = if add_chr {
                empty.update_chr_prefix("remove")
            } else {
                empty
            };
            return Ok(Archive::new(metadata, ArchiveData::Vcf(empty)));
        }
        return Err(PgxError::External(format!("beagle-rs failed: {stderr}")));
    }

    let gz = dir.join("output.vcf.gz");
    let file = std::fs::File::open(&gz)
        .map_err(|e| PgxError::External(format!("missing beagle output {}: {e}", gz.display())))?;
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(file)
        .read_to_string(&mut text)
        .map_err(io)?;
    let _ = std::fs::remove_dir_all(&dir);

    let mut phased = crate::fuc::VcfFrame::from_string(&text);
    if add_chr {
        phased = phased.update_chr_prefix("remove");
    }
    Ok(Archive::new(metadata, ArchiveData::Vcf(phased)))
}

#[cfg(feature = "beagle")]
fn beagle_command(
    bin: &str,
    input: &std::path::Path,
    region: &str,
    out_prefix: &std::path::Path,
    impute: bool,
    em: bool,
    panel: Option<&str>,
) -> std::process::Command {
    let mut cmd = std::process::Command::new(bin);
    cmd.arg(format!("gt={}", input.display()))
        .arg(format!("chrom={region}"))
        .arg(format!("out={}", out_prefix.display()))
        .arg(format!("impute={impute}"))
        .arg(format!("em={em}"));
    if let Some(p) = panel {
        cmd.arg(format!("ref={p}"));
    }
    cmd
}

/// Does the (bgzf) reference panel use `chr`-prefixed contigs? Peeks the first
/// data record. Mirrors PyPGx's `pyvcf.has_chr_prefix(panel)`.
#[cfg(feature = "beagle")]
fn panel_has_chr_prefix(panel: &str) -> std::io::Result<bool> {
    use std::io::BufRead;
    let file = std::fs::File::open(panel)?;
    let reader = std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file));
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('#') {
            return Ok(line.starts_with("chr"));
        }
    }
    Ok(false)
}

/// All `CHROM-POS-REF-ALT` variants in a (bgzf) reference panel, for the
/// input∩panel overlap check. Streams the panel and keeps only the 5 needed
/// fields: parsing it into a full `VcfFrame` would materialize every 1KGP
/// sample-genotype column (~2500 of them) and can reach multiple GB for
/// big-panel genes (e.g. DPYD) — this stays at a few MB.
#[cfg(feature = "beagle")]
fn panel_variants(panel: &str) -> std::io::Result<Vec<String>> {
    use std::io::BufRead;
    let file = std::fs::File::open(panel)?;
    let reader = std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file));
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.split('\t');
        let chrom = f.next().unwrap_or("");
        let pos = f.next().unwrap_or("");
        let _id = f.next();
        let rf = f.next().unwrap_or("");
        let alt = f.next().unwrap_or("");
        // Match VcfFrame::to_variants: split multiallelic ALT on commas.
        for a in alt.split(',') {
            out.push(format!("{chrom}-{pos}-{rf}-{a}"));
        }
    }
    Ok(out)
}

#[cfg(not(feature = "beagle"))]
deferred!(
    estimate_phase_beagle(imported_variants: &Archive, panel: Option<&str>, impute: bool) -> Archive,
    "the beagle-rs binary (enable feature `beagle`; binary via $BEAGLE_RS_BIN)",
    "`estimate_phase_beagle` — statistical phasing of a `VcfFrame[Imported]` \
     into `VcfFrame[Phased]`. Enable the `beagle` feature to shell out to the \
     beagle-rs binary (`gt=`/`chrom=`/`ref=`/`out=`/`impute=`/`em=`)."
);

// ---- Genuinely-not-yet-ported externals -----------------------------------
//
// Everything else PyPGx exposes via external tools is implemented natively and
// exported from its own module: depth/CNV (`import_read_depth`,
// `compute_copy_number`, `predict_cnv`, `test_cnv_caller`, the depth fns) live
// in `api`; the three pipelines in `pipeline`; the five plots in `plot`. Only
// these three have no native implementation yet — each needs a tool that has no
// usable Rust port — so they remain `NotPorted` stubs naming the dependency.

deferred!(
    slice_bam(input: &str, output: &str, assembly: &str, genes: Option<&[String]>, exclude: bool) -> (),
    "samtools (via pysam/pybam)",
    "`slice_bam` — subset a BAM to PyPGx's gene regions."
);
deferred!(
    create_input_vcf(vcf: &str, fasta: &str, bams: &[String], assembly: &str) -> Archive,
    "bcftools call/mpileup (not yet in bcftools-rs)",
    "`create_input_vcf` — call SNVs/indels from BAMs into a `VcfFrame[Imported]`."
);
deferred!(
    train_cnv_caller(copy_number: &Archive, cnv_calls: &Archive, confusion_matrix: Option<&str>) -> Archive,
    "scikit-learn SVM training (inherently non-parity)",
    "`train_cnv_caller` — train a `Model[CNV]` SVM classifier. No Rust trainer \
     reproduces PyPGx's libsvm-fitted models; the parity route is \
     `tools/convert_cnv_model.py` + native `api::predict_cnv`."
);

#[cfg(all(test, feature = "beagle"))]
mod tests {
    use super::beagle_command;

    #[test]
    fn beagle_args_match_pypgx_thread_default_for_g6pd_parity() {
        let cmd = beagle_command(
            "beagle-rs",
            std::path::Path::new("/tmp/input.vcf"),
            "chrX:154528389-154550018",
            std::path::Path::new("/tmp/output"),
            false,
            true,
            Some("/bundle/1kgp/GRCh38/G6PD.vcf.gz"),
        );
        let args: Vec<String> = cmd
            .get_args()
            .map(|arg| arg.to_string_lossy().into_owned())
            .collect();

        assert!(args.contains(&"gt=/tmp/input.vcf".to_string()));
        assert!(args.contains(&"chrom=chrX:154528389-154550018".to_string()));
        assert!(args.contains(&"out=/tmp/output".to_string()));
        assert!(args.contains(&"impute=false".to_string()));
        assert!(args.contains(&"em=true".to_string()));
        assert!(args.contains(&"ref=/bundle/1kgp/GRCh38/G6PD.vcf.gz".to_string()));
        assert!(
            args.iter().all(|arg| !arg.starts_with("nthreads=")),
            "PyPGx does not pass nthreads; forcing it changed G6PD HG01621/HG01808/HG02614 phase orientation: {args:?}"
        );
    }
}
