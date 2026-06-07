#!/usr/bin/env python
"""Compare the PyPGx vs pypgx-rs per-gene benchmark results.

Usage: python tools/compare_bench.py <pypgx_timing.tsv> <pypgxrs_timing.tsv> <out.tsv>

Writes a merged per-gene table (gene, seconds each, calls each, agreement) and
prints a summary: totals/speedup, and counts of match / differ / one-side-failed.
"""
import sys


def load(path):
    d, total = {}, None
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#TOTAL":
                total = float(parts[1])
                continue
            gene, sec, status = parts[0], float(parts[1]), parts[2]
            geno = parts[3] if len(parts) > 3 else ""
            pheno = parts[4] if len(parts) > 4 else ""
            d[gene] = (sec, status, geno, pheno)
    return d, total


def main():
    pp, pp_total = load(sys.argv[1])
    rs, rs_total = load(sys.argv[2])
    genes = sorted(set(pp) | set(rs))

    cats = {"match": 0, "differ": 0, "rs_failed": 0, "pypgx_failed": 0, "both_failed": 0}
    rows = []
    for g in genes:
        ps, pst, pg, pph = pp.get(g, (0, "MISSING", "", ""))
        rss, rst, rg, rph = rs.get(g, (0, "MISSING", "", ""))
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
        rows.append((g, f"{ps:.3f}", f"{rss:.3f}", pst, rst, f"{pg}|{pph}", f"{rg}|{rph}", agree))

    with open(sys.argv[3], "w") as f:
        f.write("Gene\tpypgx_s\tpypgxrs_s\tpypgx_status\tpypgxrs_status\tpypgx_call\tpypgxrs_call\tagreement\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"PyPGx total:    {pp_total:.1f}s  ({len(pp)} genes)")
    print(f"pypgx-rs total: {rs_total:.1f}s  ({len(rs)} genes)")
    if rs_total:
        print(f"speedup:        {pp_total / rs_total:.1f}x faster (pypgx-rs)")
    print(f"\nagreement: {cats}")
    both_called = cats["match"] + cats["differ"]
    if both_called:
        print(f"of {both_called} genes both called: {cats['match']} identical, {cats['differ']} differ")
    if cats["differ"]:
        print("\nDIFFERing genes (likely Beagle-version phasing differences):")
        for r in rows:
            if r[7] == "DIFFER":
                print(f"  {r[0]:10} pypgx={r[5]}   pypgx-rs={r[6]}")
    if cats["rs_failed"]:
        print("\npypgx-rs FAILED where PyPGx succeeded (port robustness gaps):")
        for r in rows:
            if r[7] == "rs-FAILED":
                print(f"  {r[0]:10} pypgx={r[5]}   pypgx-rs={r[4]}")


if __name__ == "__main__":
    main()
