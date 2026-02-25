use std::path::{Path, PathBuf};

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::stage_error::StageError;

#[derive(Debug, Clone)]
pub struct MtxPaths {
    pub matrix: PathBuf,
    pub features: PathBuf,
    pub barcodes: PathBuf,
}

pub fn detect_mtx_dir(dir: &Path) -> Result<MtxPaths, StageError> {
    if !dir.is_dir() {
        return Err(StageError::Validation(format!(
            "input path is not a directory: {}",
            dir.display()
        )));
    }

    let ds = kira_scio::discover(dir).map_err(|e| StageError::Validation(e.message))?;
    let features = ds.features.or(ds.genes).ok_or_else(|| {
        StageError::Validation("missing features.tsv(.gz) or genes.tsv(.gz)".to_string())
    })?;
    let barcodes = ds
        .barcodes
        .ok_or_else(|| StageError::Validation("missing barcodes.tsv(.gz)".to_string()))?;

    Ok(MtxPaths {
        matrix: ds.matrix,
        features,
        barcodes,
    })
}

pub fn resolve_cache_bin_path(dir: &Path) -> Result<PathBuf, StageError> {
    let prefix = detect_dataset_prefix(dir)?;
    let filename = kira_shared_sc_cache::resolve_shared_cache_filename(prefix.as_deref());
    Ok(dir.join(filename))
}

pub fn detect_dataset_prefix(dir: &Path) -> Result<Option<String>, StageError> {
    kira_scio::detect_prefix(dir).map_err(|e| StageError::Validation(e.to_string()))
}

pub fn read_features(path: &Path) -> Result<Vec<String>, StageError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| StageError::Format(e.message))?;
    Ok(md.gene_symbols)
}

pub fn read_barcodes(path: &Path) -> Result<Vec<String>, StageError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| StageError::Format(e.message))?;
    Ok(md.barcodes)
}

pub fn read_matrix_csc(
    path: &Path,
) -> Result<(usize, usize, Vec<usize>, Vec<usize>, Vec<f32>), StageError> {
    let m = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| StageError::Format(e.message))?;

    Ok((m.n_genes, m.n_cells, m.col_ptr, m.row_idx, m.values))
}

#[cfg(test)]
#[path = "../../tests/src_inline/io/mtx.rs"]
mod tests;
