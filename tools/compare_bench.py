#!/usr/bin/env python
"""Compare the PyPGx vs pypgx-rs per-gene benchmark results.

Usage: python tools/compare_bench.py <pypgx_timing.tsv> <pypgxrs_timing.tsv> <out.tsv>

Writes a merged per-gene table (gene, seconds each, calls each, agreement) and
prints a summary: totals/speedup, and counts of match / differ / one-side-failed.
"""
import sys

RESULT_COLUMNS = [
    "Genotype",
    "Phenotype",
    "Haplotype1",
    "Haplotype2",
    "AlternativePhase",
    "VariantData",
    "CNV",
]


def load(path):
    d, total = {}, None
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#TOTAL":
                total = float(parts[1])
                continue
            gene = parts[cols.get("Gene", 0)]
            sec_i = cols.get("Seconds")
            status_i = cols.get("Status", 2)
            sec = float(parts[sec_i]) if sec_i is not None else 0.0
            status = parts[status_i] if len(parts) > status_i else ""
            values = {
                col: parts[i] if (i := cols.get(col, -1)) >= 0 and len(parts) > i else ""
                for col in RESULT_COLUMNS
            }
            d[gene] = (sec, status, values)
    return d, total


def main():
    pp, pp_total = load(sys.argv[1])
    rs, rs_total = load(sys.argv[2])
    genes = sorted(set(pp) | set(rs))

    cats = {"match": 0, "differ": 0, "rs_failed": 0, "pypgx_failed": 0, "both_failed": 0}
    rows = []
    for g in genes:
        empty = {col: "" for col in RESULT_COLUMNS}
        ps, pst, pv = pp.get(g, (0, "MISSING", empty))
        rss, rst, rv = rs.get(g, (0, "MISSING", empty))
        pg, pph = pv["Genotype"], pv["Phenotype"]
        rg, rph = rv["Genotype"], rv["Phenotype"]
        p_ok, r_ok = pst == "ok", rst == "ok"
        if p_ok and r_ok:
            agree = "match" if (pg, pph) == (rg, rph) else "DIFFER"
            cats["match" if agree == "match" else "differ"] += 1
        elif p_ok and not r_ok:
            agree = "rs-FAILED"
            cats["rs_failed"] += 1
        elif r_ok and not p_ok:
            agree = "pypgx-FAILED"
            cats["pypgx_failed"] += 1
        else:
            agree = "both-failed"
            cats["both_failed"] += 1
        row = [
            g,
            f"{ps:.3f}",
            f"{rss:.3f}",
            pst,
            rst,
            f"{pg}|{pph}",
            f"{rg}|{rph}",
        ]
        row.extend(pv[col] for col in RESULT_COLUMNS[2:])
        row.extend(rv[col] for col in RESULT_COLUMNS[2:])
        row.append(agree)
        rows.append(row)

    with open(sys.argv[3], "w") as f:
        detail_cols = (
            [f"pypgx_{col}" for col in RESULT_COLUMNS[2:]]
            + [f"pypgxrs_{col}" for col in RESULT_COLUMNS[2:]]
        )
        f.write(
            "\t".join([
                "Gene",
                "pypgx_s",
                "pypgxrs_s",
                "pypgx_status",
                "pypgxrs_status",
                "pypgx_call",
                "pypgxrs_call",
                *detail_cols,
                "agreement",
            ])
            + "\n"
        )
        for r in rows:
            f.write("\t".join(r) + "\n")

    if pp_total is not None:
        print(f"PyPGx total:    {pp_total:.1f}s  ({len(pp)} genes)")
    else:
        print(f"PyPGx total:    N/A  ({len(pp)} genes)")
    if rs_total is not None:
        print(f"pypgx-rs total: {rs_total:.1f}s  ({len(rs)} genes)")
    else:
        print(f"pypgx-rs total: N/A  ({len(rs)} genes)")
    if pp_total is not None and rs_total:
        print(f"speedup:        {pp_total / rs_total:.1f}x faster (pypgx-rs)")
    print(f"\nagreement: {cats}")
    both_called = cats["match"] + cats["differ"]
    if both_called:
        print(f"of {both_called} genes both called: {cats['match']} identical, {cats['differ']} differ")
    if cats["differ"]:
        print("\nDIFFERing genes (likely Beagle-version phasing differences):")
        for r in rows:
            if r[-1] == "DIFFER":
                print(f"  {r[0]:10} pypgx={r[5]}   pypgx-rs={r[6]}")
    if cats["rs_failed"]:
        print("\npypgx-rs FAILED where PyPGx succeeded (port robustness gaps):")
        for r in rows:
            if r[-1] == "rs-FAILED":
                print(f"  {r[0]:10} pypgx={r[5]}   pypgx-rs={r[4]}")


if __name__ == "__main__":
    main()
