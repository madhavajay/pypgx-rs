#!/usr/bin/env python
"""Run the original PyPGx NGS pipeline over every target gene and time it.

Usage:
    PYPGX_BUNDLE=/home/linux/pypgx-bundle \
      python tools/bench_pypgx.py <input.vcf.gz> <out_dir> <assembly> [gene[,gene...]]

For each gene: run_ngs_pipeline(variants=<vcf>) — slices the gene region,
statistically phases against the 1KGP panel (Beagle), predicts star alleles,
genotypes, phenotypes. No depth_of_coverage -> the SV/CNV arm is skipped.

Writes <out_dir>/<gene>/ per gene and <out_dir>/timing.tsv (gene, seconds,
status, genotype, phenotype). Prints the total wall-clock.
"""
import sys
import os
import time
import warnings
import io
import contextlib

warnings.filterwarnings("ignore")
import pypgx
from pypgx import sdk


def main():
    vcf, out_dir, assembly = sys.argv[1], sys.argv[2], sys.argv[3]
    genes = (
        sys.argv[4].split(",")
        if len(sys.argv) > 4
        else pypgx.list_genes(mode="target")
    )
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    t0 = time.perf_counter()
    for gene in genes:
        gout = os.path.join(out_dir, gene)
        t = time.perf_counter()
        status, geno, pheno = "ok", "", ""
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                pypgx.run_ngs_pipeline(
                    gene, gout, variants=vcf, assembly=assembly,
                    do_not_plot_copy_number=True,
                    do_not_plot_allele_fraction=True, force=True,
                )
            res = sdk.Archive.from_file(os.path.join(gout, "results.zip")).data
            geno = ";".join(map(str, res["Genotype"].tolist()))
            pheno = ";".join(map(str, res["Phenotype"].tolist()))
        except Exception as e:  # noqa: BLE001
            status = f"ERR:{type(e).__name__}:{str(e)[:80].replace(chr(9),' ')}"
        dt = time.perf_counter() - t
        rows.append((gene, f"{dt:.3f}", status, geno, pheno))
        print(f"{gene:12} {dt:7.3f}s  {status}  {geno} / {pheno}", flush=True)
    total = time.perf_counter() - t0

    with open(os.path.join(out_dir, "timing.tsv"), "w") as f:
        f.write("Gene\tSeconds\tStatus\tGenotype\tPhenotype\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
        f.write(f"#TOTAL\t{total:.3f}\t{len(genes)} genes\n")
    ok = sum(1 for r in rows if r[2] == "ok")
    print(f"\nPyPGx: {ok}/{len(genes)} genes ok in {total:.1f}s total")


if __name__ == "__main__":
    main()
