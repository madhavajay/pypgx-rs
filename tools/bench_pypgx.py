#!/usr/bin/env python
"""Run the original PyPGx NGS pipeline over every target gene and time it.

Usage:
    PYPGX_BUNDLE=/home/linux/pypgx-bundle \
      python tools/bench_pypgx.py <input.vcf.gz> <out_dir> <assembly> [gene[,gene...]]

For each gene: run_ngs_pipeline(variants=<vcf>) — slices the gene region,
statistically phases against the 1KGP panel (Beagle), predicts star alleles,
genotypes, phenotypes. No depth_of_coverage -> the SV/CNV arm is skipped.

Writes <out_dir>/<gene>/ per gene and <out_dir>/timing.tsv. The timing table
keeps the quick genotype/phenotype columns and the detailed PyPGx result
columns so genes like COMT do not lose haplotype candidate details.
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

RESULT_COLUMNS = [
    "Genotype",
    "Phenotype",
    "Haplotype1",
    "Haplotype2",
    "AlternativePhase",
    "VariantData",
    "CNV",
]

LONG_COLUMNS = [
    "Gene",
    "Status",
    "Sample",
    "Genotype",
    "Phenotype",
    "Haplotype",
    "CandidateRank",
    "Allele",
    "Variant",
    "AlleleFraction",
    "AlternativePhase",
    "CNV",
]


def tsv_cell(value):
    return str(value).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def split_semis(value):
    value = "" if value != value else str(value)
    return [x for x in value.strip(";").split(";") if x]


def parse_variant_data(value):
    result = {}
    for entry in split_semis(value):
        fields = entry.split(":", 2)
        if len(fields) < 2:
            continue
        allele = fields[0]
        if fields[1] == "default":
            result[allele] = [("default", "")]
        elif len(fields) == 3:
            variants = fields[1].split(",")
            fractions = fields[2].split(",")
            result[allele] = [
                (variant, fractions[i] if i < len(fractions) else "")
                for i, variant in enumerate(variants)
            ]
        else:
            result[allele] = [(fields[1], "")]
    return result


def long_rows(gene, status, res):
    rows = []
    for sample, row in res.iterrows():
        genotype = "" if "Genotype" not in row else row["Genotype"]
        phenotype = "" if "Phenotype" not in row else row["Phenotype"]
        alt_phase = "" if "AlternativePhase" not in row else row["AlternativePhase"]
        cnv = "" if "CNV" not in row else row["CNV"]
        variant_data = "" if "VariantData" not in row else row["VariantData"]
        variants = parse_variant_data(variant_data)
        for haplotype in ["Haplotype1", "Haplotype2", "AlternativePhase"]:
            if haplotype not in row:
                continue
            for rank, allele in enumerate(split_semis(row[haplotype]), start=1):
                entries = variants.get(allele, [("", "")])
                for variant, fraction in entries:
                    rows.append([
                        gene,
                        status,
                        sample,
                        genotype,
                        phenotype,
                        haplotype,
                        rank,
                        allele,
                        variant,
                        fraction,
                        alt_phase,
                        cnv,
                    ])
    return rows


def main():
    vcf, out_dir, assembly = sys.argv[1], sys.argv[2], sys.argv[3]
    genes = (
        sys.argv[4].split(",")
        if len(sys.argv) > 4
        else pypgx.list_genes(mode="target")
    )
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    detail_rows = []
    t0 = time.perf_counter()
    for gene in genes:
        gout = os.path.join(out_dir, gene)
        t = time.perf_counter()
        status = "ok"
        values = [""] * len(RESULT_COLUMNS)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                pypgx.run_ngs_pipeline(
                    gene, gout, variants=vcf, assembly=assembly,
                    do_not_plot_copy_number=True,
                    do_not_plot_allele_fraction=True, force=True,
                )
            res = sdk.Archive.from_file(os.path.join(gout, "results.zip")).data
            for i, col in enumerate(RESULT_COLUMNS):
                if col in res:
                    values[i] = ";".join("" if x != x else str(x) for x in res[col].tolist())
            detail_rows.extend(long_rows(gene, status, res))
        except Exception as e:  # noqa: BLE001
            status = f"ERR:{type(e).__name__}:{str(e)[:80].replace(chr(9),' ')}"
            detail_rows.append([gene, status, "", "", "", "", "", "", "", "", "", ""])
        dt = time.perf_counter() - t
        rows.append((gene, f"{dt:.3f}", status, *values))
        print(f"{gene:12} {dt:7.3f}s  {status}  {values[0]} / {values[1]}", flush=True)
    total = time.perf_counter() - t0

    with open(os.path.join(out_dir, "timing.tsv"), "w") as f:
        f.write("Gene\tSeconds\tStatus\t" + "\t".join(RESULT_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(tsv_cell(x) for x in r) + "\n")
        f.write(f"#TOTAL\t{total:.3f}\t{len(genes)} genes\n")
    with open(os.path.join(out_dir, "timing.long.tsv"), "w") as f:
        f.write("\t".join(LONG_COLUMNS) + "\n")
        for r in detail_rows:
            f.write("\t".join(tsv_cell(x) for x in r) + "\n")
    ok = sum(1 for r in rows if r[2] == "ok")
    print(f"\nPyPGx: {ok}/{len(genes)} genes ok in {total:.1f}s total")


if __name__ == "__main__":
    main()
