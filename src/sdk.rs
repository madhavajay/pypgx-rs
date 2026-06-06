//! Port of `pypgx.sdk` — the `Archive` container, its semantic types, and the
//! exception hierarchy.

use std::io::{Read, Write};

use crate::cnv::CnvModel;
use crate::fuc::{CovFrame, VcfFrame};

/// PyPGx error hierarchy (`pypgx/sdk/utils.py`). Each Python exception maps to
/// one variant; faithful ports raise these where Python does.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PgxError {
    AlleleNotFound {
        gene: String,
        allele: String,
    },
    GeneNotFound(String),
    IncorrectMetadata(String),
    IncorrectSemanticType(String),
    NotTargetGene(String),
    PhenotypeNotFound(String),
    SemanticTypeNotFound(String),
    VariantNotFound(String),
    BundleNotFound(String),
    /// A function whose implementation is deferred because it requires an
    /// external program, the `pypgx-bundle`, sklearn, or matplotlib. The string
    /// names the missing dependency.
    NotPorted(String),
}

impl std::fmt::Display for PgxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PgxError::AlleleNotFound { gene, allele } => write!(f, "{gene}/{allele}"),
            PgxError::GeneNotFound(g) => write!(f, "{g}"),
            PgxError::IncorrectMetadata(m) => write!(f, "{m}"),
            PgxError::IncorrectSemanticType(m) => write!(f, "{m}"),
            PgxError::NotTargetGene(g) => write!(f, "{g}"),
            PgxError::PhenotypeNotFound(p) => write!(f, "{p}"),
            PgxError::SemanticTypeNotFound(s) => write!(f, "{s}"),
            PgxError::VariantNotFound(v) => write!(f, "{v}"),
            PgxError::BundleNotFound(b) => write!(f, "{b}"),
            PgxError::NotPorted(dep) => {
                write!(f, "not yet ported (requires {dep})")
            }
        }
    }
}

impl std::error::Error for PgxError {}

/// A `SampleTable[*]` payload: a sample-indexed table (pandas `DataFrame` with
/// a string index). This is what `predict_alleles` returns.
#[derive(Clone, Debug)]
pub struct SampleTable {
    pub index: Vec<String>,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl SampleTable {
    /// `DataFrame.loc[label]` — the row for a given index label.
    pub fn loc(&self, label: &str) -> &Vec<String> {
        let i = self
            .index
            .iter()
            .position(|x| x == label)
            .unwrap_or_else(|| panic!("no such index: {label}"));
        &self.rows[i]
    }
}

/// The data payload of an [`Archive`], discriminated by semantic type.
#[derive(Clone, Debug)]
pub enum ArchiveData {
    Vcf(VcfFrame),
    SampleTable(SampleTable),
    Cov(CovFrame),
    /// `Model[CNV]` — a fitted RBF OvR-SVM. PyPGx ships these as pickled sklearn
    /// objects (`data.sav`); the Rust-native form stores extracted params as
    /// `data.json` (convert with `tools/convert_cnv_model.py`).
    Model(CnvModel),
    Unsupported,
}

/// Port of `pypgx.Archive`: metadata (ordered key/value lines) plus a typed
/// data payload, serialized as a ZIP container.
#[derive(Clone, Debug)]
pub struct Archive {
    pub metadata: Vec<(String, String)>,
    pub data: ArchiveData,
}

impl Archive {
    pub fn new(metadata: Vec<(String, String)>, data: ArchiveData) -> Self {
        Archive { metadata, data }
    }

    /// `Archive.type` — the `SemanticType` metadata value.
    pub fn semantic_type(&self) -> &str {
        self.get("SemanticType").expect("SemanticType in metadata")
    }

    /// Look up a metadata value by key.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.metadata
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// `Archive.copy_metadata` — a deep copy of the metadata.
    pub fn copy_metadata(&self) -> Vec<(String, String)> {
        self.metadata.clone()
    }

    /// `Archive.check_type` — error unless the archive has one of the given
    /// semantic types.
    pub fn check_type(&self, semantic_types: &[&str]) -> Result<(), PgxError> {
        if !semantic_types.contains(&self.semantic_type()) {
            return Err(PgxError::IncorrectSemanticType(format!(
                "Expected {}, but instead found {}",
                semantic_types.join("/"),
                self.semantic_type()
            )));
        }
        Ok(())
    }

    /// `Archive.from_file` — read a ZIP archive. The first entry's top-level
    /// directory is the parent; `<parent>/metadata.txt` holds `key=value`
    /// lines and the payload file depends on the semantic type.
    pub fn from_file(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file = std::fs::File::open(path)?;
        let mut zip = zip::ZipArchive::new(file)?;

        let first = zip.by_index(0)?.name().to_string();
        let parent = first.split('/').next().unwrap_or("").to_string();

        // metadata.txt
        let mut meta_text = String::new();
        zip.by_name(&format!("{parent}/metadata.txt"))?
            .read_to_string(&mut meta_text)?;
        let mut metadata = Vec::new();
        for line in meta_text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            // Python: line.split('=') then take fields[0], fields[1].
            if let Some((k, v)) = line.split_once('=') {
                metadata.push((k.to_string(), v.to_string()));
            }
        }
        let semantic = metadata
            .iter()
            .find(|(k, _)| k == "SemanticType")
            .map(|(_, v)| v.clone())
            .ok_or_else(|| PgxError::SemanticTypeNotFound("<missing>".into()))?;

        let data = if semantic.contains("VcfFrame") {
            let mut text = String::new();
            zip.by_name(&format!("{parent}/data.vcf"))?
                .read_to_string(&mut text)?;
            ArchiveData::Vcf(VcfFrame::from_string(&text))
        } else if semantic.contains("CovFrame") {
            let mut text = String::new();
            zip.by_name(&format!("{parent}/data.tsv"))?
                .read_to_string(&mut text)?;
            ArchiveData::Cov(CovFrame::from_string(&text))
        } else if semantic.contains("SampleTable") {
            let mut text = String::new();
            zip.by_name(&format!("{parent}/data.tsv"))?
                .read_to_string(&mut text)?;
            ArchiveData::SampleTable(parse_sample_table(&text))
        } else if semantic.contains("Model") {
            let mut text = String::new();
            zip.by_name(&format!("{parent}/data.json"))?
                .read_to_string(&mut text)?;
            ArchiveData::Model(serde_json::from_str(&text)?)
        } else {
            ArchiveData::Unsupported
        };

        Ok(Archive { metadata, data })
    }

    /// `Archive.to_file` — write the archive as a ZIP. The payload directory is
    /// named after the output file's stem (PyPGx uses a temp-dir name; the dir
    /// name is irrelevant to readers, which take it from the first entry).
    pub fn to_file(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let parent = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("archive")
            .to_string();

        let file = std::fs::File::create(path)?;
        let mut zip = zip::ZipWriter::new(file);
        let opts: zip::write::FileOptions<()> =
            zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Deflated);

        // metadata.txt
        zip.start_file(format!("{parent}/metadata.txt"), opts)?;
        for (k, v) in &self.metadata {
            zip.write_all(format!("{k}={v}\n").as_bytes())?;
        }

        match &self.data {
            ArchiveData::Vcf(vcf) => {
                let mut vcf = vcf.clone();
                vcf.meta.push("##fileformat=VCFv4.2".to_string());
                zip.start_file(format!("{parent}/data.vcf"), opts)?;
                zip.write_all(vcf.to_string().as_bytes())?;
            }
            ArchiveData::SampleTable(table) => {
                zip.start_file(format!("{parent}/data.tsv"), opts)?;
                zip.write_all(serialize_sample_table(table).as_bytes())?;
            }
            ArchiveData::Cov(cov) => {
                zip.start_file(format!("{parent}/data.tsv"), opts)?;
                zip.write_all(cov.to_string().as_bytes())?;
            }
            ArchiveData::Model(model) => {
                zip.start_file(format!("{parent}/data.json"), opts)?;
                zip.write_all(serde_json::to_string(model)?.as_bytes())?;
            }
            ArchiveData::Unsupported => {
                return Err(Box::new(PgxError::SemanticTypeNotFound(
                    self.semantic_type().to_string(),
                )));
            }
        }

        zip.finish()?;
        Ok(())
    }

    /// Borrow the VcfFrame payload (panics if the archive is not a VcfFrame).
    pub fn as_vcf(&self) -> &VcfFrame {
        match &self.data {
            ArchiveData::Vcf(v) => v,
            _ => panic!("archive is not a VcfFrame"),
        }
    }

    /// Borrow the SampleTable payload (panics otherwise).
    pub fn as_sample_table(&self) -> &SampleTable {
        match &self.data {
            ArchiveData::SampleTable(t) => t,
            _ => panic!("archive is not a SampleTable"),
        }
    }

    /// Borrow the CovFrame payload (panics otherwise).
    pub fn as_cov(&self) -> &CovFrame {
        match &self.data {
            ArchiveData::Cov(c) => c,
            _ => panic!("archive is not a CovFrame"),
        }
    }

    /// Borrow the Model[CNV] payload (panics otherwise).
    pub fn as_model(&self) -> &CnvModel {
        match &self.data {
            ArchiveData::Model(m) => m,
            _ => panic!("archive is not a Model"),
        }
    }
}

/// Parse a `SampleTable` TSV (pandas `to_csv(sep='\t')` with the index in the
/// first, unnamed column).
fn parse_sample_table(text: &str) -> SampleTable {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_reader(text.as_bytes());
    let header: Vec<String> = rdr
        .headers()
        .expect("tsv header")
        .iter()
        .map(|s| s.to_string())
        .collect();
    // First column is the (unnamed) index.
    let columns = header[1..].to_vec();
    let mut index = Vec::new();
    let mut rows = Vec::new();
    for rec in rdr.records() {
        let rec = rec.expect("tsv record");
        let fields: Vec<String> = rec.iter().map(|s| s.to_string()).collect();
        index.push(fields[0].clone());
        rows.push(fields[1..].to_vec());
    }
    SampleTable {
        index,
        columns,
        rows,
    }
}

/// Serialize a `SampleTable` to TSV (pandas `to_csv(sep='\t')`: an empty header
/// cell for the index column, then one row per sample).
fn serialize_sample_table(table: &SampleTable) -> String {
    let mut out = String::new();
    out.push('\t');
    out.push_str(&table.columns.join("\t"));
    out.push('\n');
    for (i, idx) in table.index.iter().enumerate() {
        out.push_str(idx);
        for cell in &table.rows[i] {
            out.push('\t');
            out.push_str(cell);
        }
        out.push('\n');
    }
    out
}
