//! Port of the small slice of the `fuc` library that PyPGx depends on:
//! `common.parse_variant`, `common.sort_variants`, and the `pyvcf.VcfFrame`
//! methods used by `core`/`api` (`from_string`, `from_dict`, `sort`, `samples`,
//! `to_variants`, `get_af`).

/// Contig ordering used by `fuc` for genomic sorting (`pyvcf.CONTIGS`).
pub const CONTIGS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "M", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
    "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM",
];

/// VCF column headers (`pyvcf.HEADERS`).
pub const HEADERS: &[&str] = &[
    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
];

/// Parsed genomic variant: `(chrom, pos, ref, alt)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Variant {
    pub chrom: String,
    pub pos: i64,
    pub r#ref: String,
    pub alt: String,
}

/// `common.parse_variant`: split on any of `-`, `:`, `>` and take the first
/// four fields, defaulting ref/alt to empty string when absent.
pub fn parse_variant(variant: &str) -> Variant {
    let fields: Vec<&str> = variant.split(['-', ':', '>']).collect();
    let chrom = fields[0].to_string();
    let pos = fields[1].parse::<i64>().expect("variant position");
    let r#ref = fields.get(2).copied().unwrap_or("").to_string();
    let alt = fields.get(3).copied().unwrap_or("").to_string();
    Variant {
        chrom,
        pos,
        r#ref,
        alt,
    }
}

/// `common.parse_region`: parse a `chrom:start-end` region string into
/// `(chrom, start, end)`. `start`/`end` are `None` when absent (matching fuc's
/// NaN); PyPGx's gene-table regions always carry all three. Returned verbatim
/// (no 0-based conversion), matching `pybed`'s use of these values.
pub fn parse_region(region: &str) -> (String, Option<i64>, Option<i64>) {
    let (chrom, rest) = match region.split_once(':') {
        Some((c, r)) => (c.to_string(), Some(r)),
        None => (region.to_string(), None),
    };
    let (start, end) = match rest {
        Some(r) => {
            let mut it = r.split('-');
            let s = it.next().and_then(|x| x.parse::<i64>().ok());
            let e = it.next().and_then(|x| x.parse::<i64>().ok());
            (s, e)
        }
        None => (None, None),
    };
    (chrom, start, end)
}

/// Add/remove the `chr` prefix on a `chrom:start-end` region string
/// (`common.update_chr_prefix` applied to a region).
fn region_chr_prefix(region: &str, add: bool) -> String {
    if add {
        if region.starts_with("chr") {
            region.to_string()
        } else {
            format!("chr{region}")
        }
    } else {
        region.replacen("chr", "", 1)
    }
}

/// Ploidy of a genotype string = number of alleles in its `GT` field.
fn gt_ploidy(g: &str) -> usize {
    g.split(':').next().unwrap_or("").split(['/', '|']).count()
}

/// `pyvcf.gt_unphase`: unphase a genotype call. Unphased calls pass through;
/// phased calls with a missing allele just swap `|`→`/`; otherwise alleles are
/// sorted ascending and joined with `/`.
pub fn gt_unphase(g: &str) -> String {
    let l: Vec<&str> = g.split(':').collect();
    let gt = l[0];
    if !gt.contains('|') {
        return g.to_string();
    }
    if gt.contains('.') {
        return g.replace('|', "/");
    }
    let mut alleles: Vec<i64> = gt.split('|').map(|a| a.parse().unwrap()).collect();
    alleles.sort_unstable();
    let new_gt = alleles
        .iter()
        .map(|b| b.to_string())
        .collect::<Vec<_>>()
        .join("/");
    std::iter::once(new_gt)
        .chain(l[1..].iter().map(|s| s.to_string()))
        .collect::<Vec<_>>()
        .join(":")
}

/// `pyvcf.gt_het`: true when the genotype call is heterozygous.
pub fn gt_het(g: &str) -> bool {
    let gt = g.split(':').next().unwrap_or("");
    let parts: Vec<&str> = if gt.contains('/') {
        gt.split('/').collect()
    } else if gt.contains('|') {
        gt.split('|').collect()
    } else {
        return false;
    };
    parts[0] != parts[1]
}

/// `pyvcf.gt_pseudophase`: turn an unphased call into a pseudo-phased one by
/// replacing `/` with `|` in the GT field (allele order preserved).
pub fn gt_pseudophase(g: &str) -> String {
    let mut l: Vec<String> = g.split(':').map(|s| s.to_string()).collect();
    l[0] = l[0].replace('/', "|");
    l.join(":")
}

/// True when every sample call in a VCF row is phased (`|` in GT).
pub fn row_phased(samples: &[String]) -> bool {
    samples
        .iter()
        .all(|g| g.split(':').next().unwrap_or("").contains('|'))
}

/// `pyvcf.gt_diploidize`: promote a haploid genotype to diploid (`0`→`0/0`-style
/// `0/`, `.`→`./`); diploid+ calls pass through.
pub fn gt_diploidize(g: &str) -> String {
    if gt_ploidy(g) != 1 {
        return g.to_string();
    }
    let gt = g.split(':').next().unwrap_or("");
    if gt == "." {
        format!("./{g}")
    } else {
        format!("0/{g}")
    }
}

/// Sort key matching `fuc`: contig index (or `len(CONTIGS)` for unknown
/// contigs) then position, reference, alternate.
fn variant_key(v: &str) -> (usize, i64, String, String) {
    let p = parse_variant(v);
    let chrom = CONTIGS
        .iter()
        .position(|c| *c == p.chrom)
        .unwrap_or(CONTIGS.len());
    (chrom, p.pos, p.r#ref, p.alt)
}

/// `common.sort_variants`: return the variants sorted by genomic coordinate.
pub fn sort_variants<I, S>(variants: I) -> Vec<String>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let mut v: Vec<String> = variants.into_iter().map(Into::into).collect();
    v.sort_by_key(|x| variant_key(x));
    v
}

/// Format a float the way Python's `str(float)` does (e.g. `0.5` -> "0.5",
/// `1.0` -> "1.0"). Used for allele-fraction output in `predict_alleles`.
pub fn python_float_str(x: f64) -> String {
    if x.is_nan() {
        "nan".to_string()
    } else if x.is_finite() && x.fract() == 0.0 {
        format!("{x:.1}")
    } else {
        format!("{x}")
    }
}

/// Minimal `pyvcf.VcfFrame`: metadata lines plus a tabular body whose first
/// nine columns are [`HEADERS`] and whose remaining columns are samples.
#[derive(Clone, Debug)]
pub struct VcfFrame {
    pub meta: Vec<String>,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl VcfFrame {
    /// Construct directly from parts (used by `build_definition_table`).
    pub fn new(meta: Vec<String>, columns: Vec<String>, rows: Vec<Vec<String>>) -> Self {
        VcfFrame {
            meta,
            columns,
            rows,
        }
    }

    /// `VcfFrame.from_string`: parse VCF text. `##` lines are metadata, the
    /// `#CHROM` line names the columns (first renamed to `CHROM`), and data
    /// rows are tab-separated with pandas-style quoting (a fully quoted field
    /// is unquoted). Blank lines are skipped (`skip_blank_lines=True`).
    pub fn from_string(s: &str) -> Self {
        let mut meta = Vec::new();
        let mut columns: Vec<String> = Vec::new();
        let mut data_lines: Vec<&str> = Vec::new();
        let mut in_body = false;
        for line in s.split('\n') {
            if !in_body {
                if line.starts_with("##") {
                    meta.push(line.trim_end_matches(['\r']).trim().to_string());
                    continue;
                } else if line.starts_with("#CHROM") {
                    columns = line
                        .trim_end_matches(['\r'])
                        .trim()
                        .split('\t')
                        .map(|s| s.to_string())
                        .collect();
                    columns[0] = "CHROM".to_string();
                    in_body = true;
                    continue;
                }
            }
            if in_body {
                let l = line.trim_end_matches(['\r']);
                if !l.is_empty() {
                    data_lines.push(l);
                }
            }
        }
        let joined = data_lines.join("\n");
        let mut rows = Vec::new();
        if !joined.is_empty() {
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .flexible(true)
                .from_reader(joined.as_bytes());
            for rec in rdr.records() {
                let rec = rec.expect("VCF record");
                rows.push(rec.iter().map(|f| f.to_string()).collect());
            }
        }
        VcfFrame {
            meta,
            columns,
            rows,
        }
    }

    /// Column index by name.
    fn col(&self, name: &str) -> usize {
        self.columns
            .iter()
            .position(|c| c == name)
            .expect("vcf column")
    }

    /// `VcfFrame.samples`: column names after the nine fixed VCF headers.
    pub fn samples(&self) -> Vec<String> {
        self.columns[9..].to_vec()
    }

    /// `VcfFrame.subset(samples, exclude)`: keep the nine fixed columns plus the
    /// selected sample columns. With `exclude=false` samples appear in the
    /// requested order; with `exclude=true` the remaining samples keep their
    /// original order.
    pub fn subset(&self, samples: &[String], exclude: bool) -> VcfFrame {
        let all = self.samples();
        let kept: Vec<String> = if exclude {
            all.iter().filter(|s| !samples.contains(s)).cloned().collect()
        } else {
            samples
                .iter()
                .filter(|s| all.contains(s))
                .cloned()
                .collect()
        };
        let mut columns: Vec<String> = self.columns[..9].to_vec();
        columns.extend(kept.iter().cloned());
        let idx: Vec<usize> = kept.iter().map(|s| self.col(s)).collect();
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let mut nr: Vec<String> = r[..9].to_vec();
                for &i in &idx {
                    nr.push(r[i].clone());
                }
                nr
            })
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns,
            rows,
        }
    }

    fn get(&self, row: &[String], name: &str) -> String {
        row[self.col(name)].clone()
    }

    /// `VcfFrame.has_chr_prefix` — true when every contig is `chr`-prefixed.
    pub fn has_chr_prefix(&self) -> bool {
        !self.rows.is_empty() && self.rows.iter().all(|r| r[0].starts_with("chr"))
    }

    /// `VcfFrame.slice(region)` — keep variants within `chrom:start-end`,
    /// matching the frame's own contig-prefix convention.
    pub fn slice(&self, region: &str) -> VcfFrame {
        let region = region_chr_prefix(region, self.has_chr_prefix());
        let (chrom, start, end) = parse_region(&region);
        let rows = self
            .rows
            .iter()
            .filter(|r| {
                if r[0] != chrom {
                    return false;
                }
                let pos: i64 = r[1].parse().unwrap_or(i64::MIN);
                start.map(|s| pos >= s).unwrap_or(true) && end.map(|e| pos <= e).unwrap_or(true)
            })
            .cloned()
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `VcfFrame.update_chr_prefix(mode)` — add/remove the `chr` contig prefix.
    pub fn update_chr_prefix(&self, mode: &str) -> VcfFrame {
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let mut nr = r.clone();
                nr[0] = match mode {
                    "remove" => nr[0].replace("chr", ""),
                    "add" => {
                        if nr[0].contains("chr") {
                            nr[0].clone()
                        } else {
                            format!("chr{}", nr[0])
                        }
                    }
                    _ => panic!("Incorrect mode: {mode}"),
                };
                nr
            })
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `VcfFrame.duplicated(subset, keep='first')` — per-row duplicate mask over
    /// the given columns.
    pub fn duplicated(&self, subset: &[&str]) -> Vec<bool> {
        let cols: Vec<usize> = subset.iter().map(|c| self.col(c)).collect();
        let mut seen = std::collections::HashSet::new();
        self.rows
            .iter()
            .map(|r| {
                let key: Vec<String> = cols.iter().map(|&c| r[c].clone()).collect();
                !seen.insert(key)
            })
            .collect()
    }

    /// `VcfFrame.drop_duplicates(subset, keep='first')`.
    pub fn drop_duplicates(&self, subset: &[&str]) -> VcfFrame {
        let mask = self.duplicated(subset);
        let rows = self
            .rows
            .iter()
            .zip(mask)
            .filter(|(_, dup)| !*dup)
            .map(|(r, _)| r.clone())
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `VcfFrame.strip(format)` — reduce each genotype to the listed FORMAT
    /// fields (missing → `.`), blank ID/QUAL/FILTER/INFO, and set FORMAT.
    /// Metadata is dropped (`metadata=False`).
    pub fn strip(&self, format: &str) -> VcfFrame {
        let new_keys: Vec<&str> = format.split(':').collect();
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let old_keys: Vec<&str> = r[8].split(':').collect();
                let indices: Vec<Option<usize>> = new_keys
                    .iter()
                    .map(|k| old_keys.iter().position(|o| o == k))
                    .collect();
                let mut nr = r.clone();
                for &c in &[2usize, 5, 6, 7] {
                    nr[c] = ".".to_string();
                }
                for (s, cell) in nr.iter_mut().enumerate().skip(9) {
                    let old_fields: Vec<&str> = r[s].split(':').collect();
                    *cell = indices
                        .iter()
                        .map(|idx| match idx {
                            Some(i) if *i < old_fields.len() => old_fields[*i],
                            _ => ".",
                        })
                        .collect::<Vec<_>>()
                        .join(":");
                }
                nr[8] = format.to_string();
                nr
            })
            .collect();
        VcfFrame {
            meta: Vec::new(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `VcfFrame.add_af(decimals=3)` — append an `AF` field per genotype,
    /// computed from `AD` (`x/total` to 3 decimals, `.` when AD is missing/0).
    pub fn add_af(&self) -> VcfFrame {
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let ad_idx = r[8].split(':').position(|k| k == "AD");
                let mut nr = r.clone();
                for (s, cell) in nr.iter_mut().enumerate().skip(9) {
                    let af = match ad_idx {
                        None => ".".to_string(),
                        Some(i) => {
                            let fields: Vec<&str> = r[s].split(':').collect();
                            let ad = fields.get(i).copied().unwrap_or(".");
                            if ad == "." {
                                ".".to_string()
                            } else {
                                let depths: Vec<i64> =
                                    ad.split(',').map(|x| x.parse().unwrap()).collect();
                                let total: i64 = depths.iter().sum();
                                if total == 0 {
                                    ".".to_string()
                                } else {
                                    depths
                                        .iter()
                                        .map(|x| format!("{:.3}", *x as f64 / total as f64))
                                        .collect::<Vec<_>>()
                                        .join(",")
                                }
                            }
                        }
                    };
                    *cell = format!("{}:{}", r[s], af);
                }
                nr[8] = format!("{}:AF", r[8]);
                nr
            })
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `VcfFrame.unphase()` — convert phased genotypes (`0|1`) to unphased,
    /// sorted (`0/1`); see [`gt_unphase`].
    pub fn unphase(&self) -> VcfFrame {
        self.map_genotypes(gt_unphase)
    }

    /// `VcfFrame.diploidize()` — promote haploid genotypes to diploid; see
    /// [`gt_diploidize`].
    pub fn diploidize(&self) -> VcfFrame {
        self.map_genotypes(gt_diploidize)
    }

    /// `VcfFrame.pseudophase()` — mark every call phased (`/`→`|`) without
    /// reordering; see [`gt_pseudophase`]. Used for haploid MT-RNR1 variants.
    pub fn pseudophase(&self) -> VcfFrame {
        self.map_genotypes(|g| gt_pseudophase(g))
    }

    /// `VcfFrame.phased` — true when every genotype call is phased (`|`).
    pub fn phased(&self) -> bool {
        if self.rows.is_empty() {
            return false;
        }
        self.rows.iter().all(|r| {
            r[9..]
                .iter()
                .all(|g| g.split(':').next().unwrap_or("").contains('|'))
        })
    }

    /// `VcfFrame.fetch(variant)` — the row whose `CHROM-POS-REF-ALT` equals
    /// `variant` (whole ALT, not split), or `None`.
    pub fn fetch(&self, variant: &str) -> Option<Vec<String>> {
        self.rows
            .iter()
            .find(|r| format!("{}-{}-{}-{}", r[0], r[1], r[3], r[4]) == variant)
            .cloned()
    }

    /// `VcfFrame.filter_vcf(other, opposite)` — keep rows whose
    /// `(CHROM, POS, REF, ALT)` key appears in `other` (or, with `opposite`, the
    /// rows that do not).
    pub fn filter_vcf(&self, other: &VcfFrame, opposite: bool) -> VcfFrame {
        let keys: std::collections::HashSet<(String, String, String, String)> = other
            .rows
            .iter()
            .map(|r| (r[0].clone(), r[1].clone(), r[3].clone(), r[4].clone()))
            .collect();
        let rows = self
            .rows
            .iter()
            .filter(|r| {
                let present = keys.contains(&(r[0].clone(), r[1].clone(), r[3].clone(), r[4].clone()));
                present != opposite
            })
            .cloned()
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// Apply a per-genotype transform to every sample cell.
    fn map_genotypes(&self, f: fn(&str) -> String) -> VcfFrame {
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let mut nr = r.clone();
                for cell in nr.iter_mut().skip(9) {
                    *cell = f(cell);
                }
                nr
            })
            .collect();
        VcfFrame {
            meta: self.meta.clone(),
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `VcfFrame.to_variants`: every `CHROM-POS-REF-ALT` (multiallelic ALT is
    /// split on commas), in row order, duplicates preserved.
    pub fn to_variants(&self) -> Vec<String> {
        let mut out = Vec::new();
        for r in &self.rows {
            let chrom = self.get(r, "CHROM");
            let pos = self.get(r, "POS");
            let rf = self.get(r, "REF");
            let alt = self.get(r, "ALT");
            for a in alt.split(',') {
                out.push(format!("{chrom}-{pos}-{rf}-{a}"));
            }
        }
        out
    }

    /// `VcfFrame.get_af`: allele fraction for a (sample, variant) pair, or
    /// `None` (NaN) when the variant is absent or has no `AF` FORMAT field.
    pub fn get_af(&self, sample: &str, variant: &str) -> Option<f64> {
        let v = parse_variant(variant);
        let row = self.rows.iter().find(|r| {
            self.get(r, "CHROM") == v.chrom
                && self.get(r, "POS") == v.pos.to_string()
                && self.get(r, "REF") == v.r#ref
                && self.get(r, "ALT").split(',').any(|a| a == v.alt)
        })?;
        let format = self.get(row, "FORMAT");
        let af_idx = format.split(':').position(|k| k == "AF")?;
        let alt_field = self.get(row, "ALT");
        let j = alt_field.split(',').position(|a| a == v.alt)?;
        let field_full = row[self.col(sample)].clone();
        let field = field_full.split(':').nth(af_idx)?;
        if field == "." {
            return None;
        }
        // AF lists ref + each alt; skip the reference value with `j + 1`.
        field
            .split(',')
            .nth(j + 1)
            .and_then(|x| x.parse::<f64>().ok())
    }

    /// `VcfFrame.sort`: order rows by contig index then position (stable).
    pub fn sort(mut self) -> Self {
        let chrom_c = self.col("CHROM");
        let pos_c = self.col("POS");
        self.rows.sort_by(|a, b| {
            let ka = CONTIGS
                .iter()
                .position(|c| *c == a[chrom_c])
                .unwrap_or(CONTIGS.len());
            let kb = CONTIGS
                .iter()
                .position(|c| *c == b[chrom_c])
                .unwrap_or(CONTIGS.len());
            let pa: i64 = a[pos_c].parse().unwrap_or(i64::MAX);
            let pb: i64 = b[pos_c].parse().unwrap_or(i64::MAX);
            (ka, pa).cmp(&(kb, pb))
        });
        self
    }
}

/// `VcfFrame.to_string`: metadata lines, the `#CHROM` header (first column
/// re-prefixed with `#`), then tab-separated rows. Implemented via `Display`
/// so `vcf.to_string()` works as in PyPGx.
impl std::fmt::Display for VcfFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for m in &self.meta {
            writeln!(f, "{m}")?;
        }
        writeln!(f, "#{}", self.columns.join("\t"))?;
        for r in &self.rows {
            writeln!(f, "{}", r.join("\t"))?;
        }
        Ok(())
    }
}

/// `pycov.CovFrame`: a per-position coverage table — `Chromosome`, `Position`,
/// then one read-depth column per sample. Only the slice pypgx touches is
/// ported (`slice`/`subset`/`update_chr_prefix`); `from_bam` lives behind the
/// `bam` feature.
#[derive(Clone, Debug)]
pub struct CovFrame {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl CovFrame {
    /// Parse the archive's `data.tsv` (`Chromosome`, `Position`, samples…).
    pub fn from_string(s: &str) -> Self {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_reader(s.as_bytes());
        let columns = rdr
            .headers()
            .expect("cov header")
            .iter()
            .map(|x| x.to_string())
            .collect();
        let rows = rdr
            .records()
            .map(|rec| rec.expect("cov row").iter().map(|f| f.to_string()).collect())
            .collect();
        CovFrame { columns, rows }
    }

    /// Sample (depth) column names — everything after `Chromosome`, `Position`.
    pub fn samples(&self) -> Vec<String> {
        self.columns[2..].to_vec()
    }

    /// True when every contig is `chr`-prefixed.
    pub fn has_chr_prefix(&self) -> bool {
        !self.rows.is_empty() && self.rows.iter().all(|r| r[0].starts_with("chr"))
    }

    /// `CovFrame.update_chr_prefix(mode)` — add/remove the `chr` contig prefix.
    pub fn update_chr_prefix(&self, mode: &str) -> CovFrame {
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let mut nr = r.clone();
                nr[0] = match mode {
                    "remove" => nr[0].replace("chr", ""),
                    "add" => {
                        if nr[0].contains("chr") {
                            nr[0].clone()
                        } else {
                            format!("chr{}", nr[0])
                        }
                    }
                    _ => panic!("Incorrect mode: {mode}"),
                };
                nr
            })
            .collect();
        CovFrame {
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `CovFrame.slice(region)` — keep positions within `chrom:start-end`.
    pub fn slice(&self, region: &str) -> CovFrame {
        let region = region_chr_prefix(region, self.has_chr_prefix());
        let (chrom, start, end) = parse_region(&region);
        let rows = self
            .rows
            .iter()
            .filter(|r| {
                if r[0] != chrom {
                    return false;
                }
                let pos: i64 = r[1].parse().unwrap_or(i64::MIN);
                start.map(|s| pos >= s).unwrap_or(true) && end.map(|e| pos <= e).unwrap_or(true)
            })
            .cloned()
            .collect();
        CovFrame {
            columns: self.columns.clone(),
            rows,
        }
    }

    /// `CovFrame.subset(samples, exclude)` — keep `Chromosome`, `Position` plus
    /// the selected sample columns (requested order on include; original order
    /// on exclude).
    pub fn subset(&self, samples: &[String], exclude: bool) -> CovFrame {
        let all = self.samples();
        let kept: Vec<String> = if exclude {
            all.iter().filter(|s| !samples.contains(s)).cloned().collect()
        } else {
            samples
                .iter()
                .filter(|s| all.contains(s))
                .cloned()
                .collect()
        };
        let mut columns: Vec<String> = self.columns[..2].to_vec();
        columns.extend(kept.iter().cloned());
        let idx: Vec<usize> = kept
            .iter()
            .map(|s| self.columns.iter().position(|c| c == s).unwrap())
            .collect();
        let rows = self
            .rows
            .iter()
            .map(|r| {
                let mut nr: Vec<String> = r[..2].to_vec();
                for &i in &idx {
                    nr.push(r[i].clone());
                }
                nr
            })
            .collect();
        CovFrame { columns, rows }
    }
}

/// Serialize a `CovFrame` to TSV (`Chromosome\tPosition\t…`, one row per line).
impl std::fmt::Display for CovFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}", self.columns.join("\t"))?;
        for r in &self.rows {
            writeln!(f, "{}", r.join("\t"))?;
        }
        Ok(())
    }
}

#[cfg(feature = "bam")]
impl CovFrame {
    /// `pycov.CovFrame.from_bam` — per-position read depth across `region` for
    /// each BAM (one depth column per BAM; sample name = file stem). Uses
    /// samtools-rs's `native::depth` (verified against C samtools).
    ///
    /// NOTE: byte-parity with PyPGx is **unverified** in this environment —
    /// PyPGx computes depth via pysam, which may differ from `samtools depth`,
    /// and there are no BAM fixtures here to diff against.
    pub fn from_bam(bams: &[String], region: &str, zero: bool) -> std::io::Result<CovFrame> {
        use std::collections::{BTreeSet, HashMap};
        let mut sample_names = Vec::new();
        let mut per_sample: Vec<HashMap<(String, usize), u32>> = Vec::new();
        let mut positions: BTreeSet<(usize, String, usize)> = BTreeSet::new();
        for bam in bams {
            let name = std::path::Path::new(bam)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("sample")
                .to_string();
            sample_names.push(name);
            let depths = samtools_rs::native::depth(bam, region, zero, None)?;
            let mut m = HashMap::new();
            for d in depths {
                let contig_rank = CONTIGS
                    .iter()
                    .position(|c| *c == d.reference_name)
                    .unwrap_or(CONTIGS.len());
                positions.insert((contig_rank, d.reference_name.clone(), d.position));
                m.insert((d.reference_name, d.position), d.depth);
            }
            per_sample.push(m);
        }
        let mut columns = vec!["Chromosome".to_string(), "Position".to_string()];
        columns.extend(sample_names.iter().cloned());
        let rows = positions
            .into_iter()
            .map(|(_, chrom, pos)| {
                let mut r = vec![chrom.clone(), pos.to_string()];
                for m in &per_sample {
                    r.push(m.get(&(chrom.clone(), pos)).copied().unwrap_or(0).to_string());
                }
                r
            })
            .collect();
        Ok(CovFrame { columns, rows })
    }
}
