use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::stage_error::StageError;

pub fn read_obs_var(path: &Path) -> Result<(Vec<String>, Vec<String>), StageError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| StageError::Format(e.message))?;

    Ok((md.barcodes, md.gene_symbols))
}

pub fn read_sparse_matrix(
    path: &Path,
) -> Result<(usize, usize, Vec<usize>, Vec<usize>, Vec<f32>), StageError> {
    let m = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| StageError::Format(e.message))?;

    Ok((m.n_genes, m.n_cells, m.col_ptr, m.row_idx, m.values))
}
