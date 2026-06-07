//! Smoke + error-handling tests for the `pypgx` CLI binary. These exercise the
//! clap dispatch (previously untested) and, crucially, confirm that bad input
//! produces a clean `error: ...` message with exit code 1 — never a Rust panic
//! backtrace (the dispatch is wrapped in a catch_unwind in `main`).

use std::process::Command;

fn pypgx() -> Command {
    Command::new(env!("CARGO_BIN_EXE_pypgx"))
}

fn run(args: &[&str]) -> (Option<i32>, String, String) {
    let out = pypgx().args(args).output().expect("spawn pypgx");
    (
        out.status.code(),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn list_genes_succeeds() {
    let (code, stdout, _) = run(&["list-genes", "--mode", "target"]);
    assert_eq!(code, Some(0));
    assert_eq!(stdout.lines().count(), 88, "expected 88 target genes");
}

#[test]
fn unknown_gene_errors_cleanly_not_panic() {
    let (code, _, stderr) = run(&["list-alleles", "NOT_A_GENE"]);
    assert_eq!(code, Some(1), "stderr: {stderr}");
    assert!(stderr.starts_with("error:"), "stderr: {stderr}");
    assert!(!stderr.contains("panicked"), "must not panic: {stderr}");
}

#[test]
fn combine_results_no_input_errors_cleanly() {
    let (code, _, stderr) = run(&["combine-results"]);
    assert_eq!(code, Some(1), "stderr: {stderr}");
    assert!(stderr.contains("No input data"), "stderr: {stderr}");
    assert!(!stderr.contains("panicked"), "must not panic: {stderr}");
}

#[test]
fn missing_archive_errors_cleanly() {
    let (code, _, stderr) = run(&["predict-alleles", "/no/such/file.zip"]);
    assert_eq!(code, Some(1), "stderr: {stderr}");
    assert!(stderr.starts_with("error:"), "stderr: {stderr}");
    assert!(!stderr.contains("panicked"), "must not panic: {stderr}");
}

#[test]
fn call_phenotypes_on_missing_file_errors_cleanly() {
    let (code, _, stderr) = run(&["call-phenotypes", "/no/such/genotypes.zip"]);
    assert_eq!(code, Some(1), "stderr: {stderr}");
    assert!(!stderr.contains("panicked"), "must not panic: {stderr}");
}
