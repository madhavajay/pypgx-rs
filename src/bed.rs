//! Minimal port of the `fuc.pybed.BedFrame` surface PyPGx touches: a BED-style
//! table (`Chromosome`, `Start`, `End`, optional `Name`) with the
//! `update_chr_prefix` and `merge` operations used by `create_regions_bed`.
//!
//! Coordinates and ordering are reproduced to match `pybed`/`pyranges` exactly
//! (verified against `tests/fixtures/regions_bed.json`): chromosomes sort
//! numerically first, then the remainder lexicographically (so `M` precedes
//! `X`); within a chromosome the input (gene-table) row order is preserved.

/// `pyranges` chromosome sort key: numeric contigs ascending, then the rest
/// lexicographically. Stable sorting on this key keeps the original row order
/// within each chromosome.
pub fn chrom_sort_key(chrom: &str) -> (u8, i64, String) {
    match chrom.parse::<i64>() {
        Ok(n) => (0, n, String::new()),
        Err(_) => (1, 0, chrom.to_string()),
    }
}

/// A BED table mirroring `pybed.BedFrame.gr.df`: named columns plus string-typed
/// rows (matching pandas' `astype(str)` rendering used in the ground truth).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BedFrame {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl BedFrame {
    /// Build from `(chrom, start, end, name)` tuples, stable-sorted by the
    /// `pyranges` chromosome order (input order preserved within a chromosome).
    pub fn from_regions(mut data: Vec<(String, i64, i64, String)>) -> Self {
        data.sort_by(|a, b| chrom_sort_key(&a.0).cmp(&chrom_sort_key(&b.0)));
        let rows = data
            .into_iter()
            .map(|(c, s, e, n)| vec![c, s.to_string(), e.to_string(), n])
            .collect();
        BedFrame {
            columns: ["Chromosome", "Start", "End", "Name"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            rows,
        }
    }

    /// `BedFrame.update_chr_prefix(mode='add')` — prepend `chr` to each contig.
    pub fn add_chr_prefix(&mut self) {
        for r in &mut self.rows {
            r[0] = format!("chr{}", r[0]);
        }
    }

    /// `BedFrame.merge()` — merge overlapping/bookended intervals (dropping the
    /// `Name` column), matching `pyranges` merge with the default `slack=0`
    /// (merge when `next.start <= prev.end`).
    pub fn merge(&self) -> BedFrame {
        let mut items: Vec<(String, i64, i64)> = self
            .rows
            .iter()
            .map(|r| (r[0].clone(), r[1].parse().unwrap(), r[2].parse().unwrap()))
            .collect();
        items.sort_by(|a, b| {
            chrom_sort_key(&a.0)
                .cmp(&chrom_sort_key(&b.0))
                .then(a.1.cmp(&b.1))
        });
        let mut out: Vec<(String, i64, i64)> = Vec::new();
        for (c, s, e) in items {
            if let Some(last) = out.last_mut() {
                if last.0 == c && s <= last.2 {
                    if e > last.2 {
                        last.2 = e;
                    }
                    continue;
                }
            }
            out.push((c, s, e));
        }
        BedFrame {
            columns: ["Chromosome", "Start", "End"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            rows: out
                .into_iter()
                .map(|(c, s, e)| vec![c, s.to_string(), e.to_string()])
                .collect(),
        }
    }

    /// `BedFrame.to_regions()` — `chrom:start-end` strings, one per interval.
    pub fn to_regions(&self) -> Vec<String> {
        self.rows
            .iter()
            .map(|r| format!("{}:{}-{}", r[0], r[1], r[2]))
            .collect()
    }

    /// Render as BED text (tab-separated, one interval per line, no header).
    pub fn to_bed_string(&self) -> String {
        let mut out = String::new();
        for r in &self.rows {
            out.push_str(&r.join("\t"));
            out.push('\n');
        }
        out
    }
}
