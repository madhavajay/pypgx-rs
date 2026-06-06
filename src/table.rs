//! A minimal pandas-like data frame with faithful NaN semantics.
//!
//! PyPGx loads its reference tables with `pandas.read_csv`. To reproduce its
//! behavior byte-for-byte we mirror two pandas details exactly:
//!
//! 1. **NA tokens.** With the default `na_filter=True`, pandas converts a fixed
//!    set of tokens (including `''`, `'NA'`, `'N/A'`, `'None'`, `'NaN'`, ...) to
//!    `NaN`. The recommendation table is loaded with `na_filter=False`, so no
//!    conversion happens there.
//! 2. **Row/column order** is preserved as in the source CSV.

/// The pandas default `STR_NA_VALUES` set (as of the reference environment).
/// Cells equal to one of these become [`Cell::Null`] when `na_filter` is true.
pub const NA_VALUES: &[&str] = &[
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN",
    "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null",
];

/// A single cell: either a string value or NaN (missing).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Cell {
    Str(String),
    Null,
}

impl Cell {
    /// `Some(&str)` for a value, `None` for NaN. Mirrors a non-NaN check.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Cell::Str(s) => Some(s.as_str()),
            Cell::Null => None,
        }
    }

    /// True when the cell is NaN — equivalent to `pandas.isna(cell)`.
    pub fn is_null(&self) -> bool {
        matches!(self, Cell::Null)
    }

    /// Interpret a boolean-typed column. pandas infers a bool dtype from the
    /// CSV tokens `True`/`TRUE`/`true` (and the false spellings); match those.
    pub fn is_true(&self) -> bool {
        matches!(self, Cell::Str(s) if s == "True" || s == "TRUE" || s == "true")
    }
}

/// A column-named, row-ordered table.
#[derive(Clone, Debug)]
pub struct Frame {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<Cell>>,
}

impl Frame {
    /// Parse CSV text the way `pandas.read_csv` would.
    ///
    /// `na_filter`: when true, [`NA_VALUES`] tokens become [`Cell::Null`];
    /// when false (as for the recommendation table) every field is kept as-is.
    pub fn from_csv(text: &str, na_filter: bool) -> Frame {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(true)
            .flexible(false)
            .from_reader(text.as_bytes());
        let columns: Vec<String> = rdr
            .headers()
            .expect("CSV header")
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut rows = Vec::new();
        for rec in rdr.records() {
            let rec = rec.expect("CSV record");
            let row: Vec<Cell> = rec
                .iter()
                .map(|f| {
                    if na_filter && NA_VALUES.contains(&f) {
                        Cell::Null
                    } else {
                        Cell::Str(f.to_string())
                    }
                })
                .collect();
            rows.push(row);
        }
        Frame { columns, rows }
    }

    /// Index of a named column.
    pub fn col(&self, name: &str) -> usize {
        self.columns
            .iter()
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("no such column: {name}"))
    }

    /// Borrow a cell by row index and column name.
    pub fn at(&self, row: usize, name: &str) -> &Cell {
        &self.rows[row][self.col(name)]
    }

    /// All values of a column, in row order (NaN included as `Cell::Null`).
    pub fn column(&self, name: &str) -> Vec<&Cell> {
        let c = self.col(name);
        self.rows.iter().map(|r| &r[c]).collect()
    }

    /// Rows where `name == value` (string equality; NaN never matches).
    pub fn filter_eq(&self, name: &str, value: &str) -> Vec<&Vec<Cell>> {
        let c = self.col(name);
        self.rows
            .iter()
            .filter(|r| r[c].as_str() == Some(value))
            .collect()
    }

    /// `pandas.Series.unique()` over a column: distinct values in order of
    /// first appearance. NaN is included once if present (as `None`).
    pub fn unique(&self, name: &str) -> Vec<Option<String>> {
        let c = self.col(name);
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for r in &self.rows {
            let key = r[c].as_str().map(|s| s.to_string());
            if seen.insert(key.clone()) {
                out.push(key);
            }
        }
        out
    }

    /// `value_counts()[value]` — number of rows whose column equals `value`.
    pub fn value_count(&self, name: &str, value: &str) -> usize {
        let c = self.col(name);
        self.rows
            .iter()
            .filter(|r| r[c].as_str() == Some(value))
            .count()
    }
}
