//! Port of `pypgx.api.plot` — the five diagnostic plots. Rendering uses the
//! pure-Rust `ruviz` (tiny-skia + cosmic-text) instead of matplotlib, so parity
//! is **visual, not pixel-exact**. Each function writes one PNG per sample and
//! returns the written paths.
//!
//! Gated behind the `plots` feature so the analytical core stays lean. The thin
//! exon-annotation track (`_plot_exons`) and the `fitted` copy-number overlay
//! (`_process_copy_number`) are not reproduced (ruviz lacks gridspec height
//! ratios / rectangle primitives); the main data panel is faithful.
#![cfg(feature = "plots")]

use ruviz::prelude::*;

use crate::core;
use crate::fuc::{parse_region, VcfFrame};
use crate::sdk::Archive;

type PlotResult = std::result::Result<Vec<String>, Box<dyn std::error::Error>>;

fn out_path(path: Option<&str>, sample: &str) -> String {
    match path {
        Some(p) => format!("{p}/{sample}.png"),
        None => format!("{sample}.png"),
    }
}

fn region_bounds(gene: &str, assembly: &str) -> (f64, f64) {
    let region = core::get_region(gene, assembly).expect("region");
    let (_, start, end) = parse_region(&region);
    (start.unwrap_or(0) as f64, end.unwrap_or(0) as f64)
}

/// Resolve the requested samples (or all), as owned strings.
fn resolve(all: Vec<String>, samples: Option<&[String]>) -> Vec<String> {
    samples.map(|s| s.to_vec()).unwrap_or(all)
}

/// Shared per-position line plot for a CovFrame column (read depth / copy number).
fn cov_line_plots(
    archive: &Archive,
    ylabel: &str,
    path: Option<&str>,
    samples: Option<&[String]>,
) -> PlotResult {
    let cf = archive.as_cov();
    let gene = archive.get("Gene").expect("Gene");
    let assembly = archive.get("Assembly").expect("Assembly");
    let (start, end) = region_bounds(gene, assembly);
    let chosen = resolve(cf.samples(), samples);
    let mut written = Vec::new();
    for sample in &chosen {
        let col = cf.columns.iter().position(|c| c == sample).expect("sample");
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for r in &cf.rows {
            let pos: f64 = r[1].parse().unwrap_or(f64::NAN);
            if pos < start || pos > end {
                continue;
            }
            if let Ok(v) = r[col].parse::<f64>() {
                xs.push(pos);
                ys.push(v);
            }
        }
        let output = out_path(path, sample);
        Plot::new()
            .line(&xs, &ys)
            .title(&format!("{sample} — {gene}"))
            .xlabel("Position (bp)")
            .ylabel(ylabel)
            .save(&output)?;
        written.push(output);
    }
    Ok(written)
}

/// `plot_bam_read_depth` — per-position read depth from a `CovFrame[ReadDepth]`.
pub fn plot_bam_read_depth(read_depth: &Archive, path: Option<&str>, samples: Option<&[String]>) -> PlotResult {
    read_depth.check_type(&["CovFrame[ReadDepth]"])?;
    cov_line_plots(read_depth, "Read depth", path, samples)
}

/// `plot_bam_copy_number` — per-position copy number from `CovFrame[CopyNumber]`
/// (the `fitted` overlay is not reproduced).
pub fn plot_bam_copy_number(copy_number: &Archive, path: Option<&str>, samples: Option<&[String]>) -> PlotResult {
    copy_number.check_type(&["CovFrame[CopyNumber]"])?;
    cov_line_plots(copy_number, "Copy number", path, samples)
}

/// Extract a FORMAT subfield (e.g. `DP`) per variant for one sample.
fn vcf_field(vf: &VcfFrame, sample_col: usize, key: &str) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for r in &vf.rows {
        let pos: f64 = match r[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let idx = r[8].split(':').position(|k| k == key);
        if let Some(i) = idx {
            if let Some(field) = r[sample_col].split(':').nth(i) {
                if let Ok(v) = field.parse::<f64>() {
                    out.push((pos, v));
                }
            }
        }
    }
    out
}

/// `plot_vcf_read_depth` — per-variant read depth (`DP`) from a VcfFrame.
pub fn plot_vcf_read_depth(
    gene: &str,
    vcf: &VcfFrame,
    assembly: &str,
    path: Option<&str>,
    samples: Option<&[String]>,
) -> PlotResult {
    let (start, end) = region_bounds(gene, assembly);
    let vf = vcf.slice(&core::get_region(gene, assembly)?);
    let chosen = resolve(vf.samples(), samples);
    let mut written = Vec::new();
    for sample in &chosen {
        let col = 9 + vf.samples().iter().position(|s| s == sample).expect("sample");
        let pts: Vec<(f64, f64)> = vcf_field(&vf, col, "DP")
            .into_iter()
            .filter(|(p, _)| *p >= start && *p <= end)
            .collect();
        let (xs, ys): (Vec<f64>, Vec<f64>) = pts.into_iter().unzip();
        let output = out_path(path, sample);
        Plot::new()
            .scatter(&xs, &ys)
            .title(&format!("{sample} — {gene}"))
            .xlabel("Position (bp)")
            .ylabel("Read depth")
            .save(&output)?;
        written.push(output);
    }
    Ok(written)
}

/// REF and ALT allele fractions per variant for one sample, from `AD`.
fn allele_fractions(vf: &VcfFrame, sample_col: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut xs, mut refs, mut alts) = (Vec::new(), Vec::new(), Vec::new());
    for r in &vf.rows {
        let pos: f64 = match r[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let Some(ad_i) = r[8].split(':').position(|k| k == "AD") else {
            continue;
        };
        let Some(ad) = r[sample_col].split(':').nth(ad_i) else {
            continue;
        };
        let depths: Vec<f64> = ad.split(',').filter_map(|x| x.parse().ok()).collect();
        let total: f64 = depths.iter().sum();
        if depths.len() < 2 || total == 0.0 {
            continue;
        }
        xs.push(pos);
        refs.push(depths[0] / total);
        alts.push(depths[1] / total);
    }
    (xs, refs, alts)
}

/// `plot_vcf_allele_fraction` — REF/ALT allele fractions from an Imported/
/// Consolidated VcfFrame archive.
pub fn plot_vcf_allele_fraction(imported_variants: &Archive, path: Option<&str>, samples: Option<&[String]>) -> PlotResult {
    imported_variants.check_type(&["VcfFrame[Imported]", "VcfFrame[Consolidated]"])?;
    let vf = imported_variants.as_vcf();
    let gene = imported_variants.get("Gene").expect("Gene");
    let chosen = resolve(vf.samples(), samples);
    let mut written = Vec::new();
    for sample in &chosen {
        let col = 9 + vf.samples().iter().position(|s| s == sample).expect("sample");
        let (xs, refs, alts) = allele_fractions(vf, col);
        let output = out_path(path, sample);
        Plot::new()
            .scatter(&xs, &refs)
            .label("REF")
            .scatter(&xs, &alts)
            .label("ALT")
            .legend(Position::TopRight)
            .title(&format!("{sample} — {gene}"))
            .xlabel("Position (bp)")
            .ylabel("Allele fraction")
            .save(&output)?;
        written.push(output);
    }
    Ok(written)
}

/// `plot_cn_af` — copy number (line) above allele fraction (scatter), per sample.
pub fn plot_cn_af(copy_number: &Archive, imported_variants: &Archive, path: Option<&str>, samples: Option<&[String]>) -> PlotResult {
    copy_number.check_type(&["CovFrame[CopyNumber]"])?;
    imported_variants.check_type(&["VcfFrame[Imported]", "VcfFrame[Consolidated]"])?;
    let cf = copy_number.as_cov();
    let vf = imported_variants.as_vcf();
    let gene = copy_number.get("Gene").expect("Gene");
    let chosen = resolve(cf.samples(), samples);
    let mut written = Vec::new();
    for sample in &chosen {
        // Copy-number panel.
        let cn_col = cf.columns.iter().position(|c| c == sample).expect("sample");
        let (mut cx, mut cy) = (Vec::new(), Vec::new());
        for r in &cf.rows {
            if let (Ok(p), Ok(v)) = (r[1].parse::<f64>(), r[cn_col].parse::<f64>()) {
                cx.push(p);
                cy.push(v);
            }
        }
        let cn_plot = Plot::new()
            .line(&cx, &cy)
            .title(&format!("{sample} — {gene} copy number"))
            .xlabel("Position (bp)")
            .ylabel("Copy number");
        // Allele-fraction panel.
        let af_col = 9 + vf.samples().iter().position(|s| s == sample).unwrap_or(0);
        let (ax, aref, aalt) = allele_fractions(vf, af_col);
        let af_plot = Plot::new()
            .scatter(&ax, &aref)
            .label("REF")
            .scatter(&ax, &aalt)
            .label("ALT")
            .legend(Position::TopRight)
            .xlabel("Position (bp)")
            .ylabel("Allele fraction");

        let output = out_path(path, sample);
        subplots(2, 1, 1000, 800)?
            .subplot_at(0, cn_plot.into())?
            .subplot_at(1, af_plot.into())?
            .save(&output)?;
        written.push(output);
    }
    Ok(written)
}
