use crate::expr::normalized::{CscMatrix, CscNormalized, NormalizedExpr};
use crate::io::{cache_bin, dense_tsv, mtx};
use crate::model::cli::RunMode;
use crate::model::ctx::{Ctx, InputType};
use crate::stage_error::StageError;
use tracing::warn;

#[cfg(feature = "h5ad")]
use crate::io::h5ad;

const DEFAULT_SCALE: f32 = 1e4;

pub fn run_stage1(ctx: &mut Ctx) -> Result<(), StageError> {
    let input_type = ctx
        .input_type
        .ok_or_else(|| StageError::Validation("input_type missing".to_string()))?;
    let input_path = ctx
        .input_path
        .as_ref()
        .ok_or_else(|| StageError::Validation("input_path missing".to_string()))?;

    let scale = DEFAULT_SCALE;
    if matches!(ctx.run_mode, RunMode::Pipeline)
        && matches!(input_type, InputType::Mtx10x | InputType::DenseTsv)
    {
        let cache_candidate = if let Some(path) = ctx.cache_path_override.clone() {
            Some(path)
        } else if input_path.is_dir() {
            Some(mtx::resolve_cache_bin_path(input_path)?)
        } else {
            None
        };
        if let Some(cache_path) = cache_candidate {
            if cache_path.is_file() {
                return load_from_pipeline_cache(ctx, &cache_path, scale);
            }
            if matches!(input_type, InputType::Mtx10x) && input_path.is_dir() {
                // Preserve legacy parsing behavior for MTX discovery in pipeline mode,
                // but fail explicitly because shared cache is required in this mode.
                let paths = mtx::detect_mtx_dir(input_path)?;
                let _ = mtx::read_barcodes(&paths.barcodes)?;
                let _ = mtx::read_features(&paths.features)?;
                let _ = mtx::read_matrix_csc(&paths.matrix)?;
                return Err(StageError::Validation(format!(
                    "pipeline mode expects shared cache file '{}' (produced by kira-mitoqc)",
                    cache_path.display()
                )));
            }
            warn!(
                expected_cache = %cache_path.display(),
                "shared cache not found in pipeline mode, falling back to direct input"
            );
        }
    }

    let matrix = match input_type {
        InputType::Mtx10x => {
            let paths = mtx::detect_mtx_dir(input_path)?;
            let (n_rows, n_cols, indptr, indices, data) = mtx::read_matrix_csc(&paths.matrix)?;
            CscMatrix {
                n_vars: n_rows,
                n_obs: n_cols,
                indptr,
                indices,
                data,
            }
        }
        InputType::DenseTsv => dense_tsv::load_dense(input_path)?.matrix,
        InputType::H5ad => {
            #[cfg(feature = "h5ad")]
            {
                let (n_rows, n_cols, indptr, indices, data) = h5ad::read_sparse_matrix(input_path)?;
                CscMatrix {
                    n_vars: n_rows,
                    n_obs: n_cols,
                    indptr,
                    indices,
                    data,
                }
            }
            #[cfg(not(feature = "h5ad"))]
            {
                return Err(StageError::FeatureDisabled(
                    "h5ad support requires the 'h5ad' feature".to_string(),
                ));
            }
        }
    };

    if matrix.n_obs != ctx.n_obs {
        return Err(StageError::Validation(
            "matrix n_obs does not match Stage 0 obs_ids".to_string(),
        ));
    }
    if matrix.n_vars != ctx.n_vars {
        return Err(StageError::Validation(
            "matrix n_vars does not match Stage 0 gene symbols".to_string(),
        ));
    }

    let normalized = CscNormalized::new(matrix, scale)?;
    let zero_lib = normalized.zero_libsize_count();
    if zero_lib > 0 {
        warn!("{zero_lib} observations have zero libsize");
    }
    let nnz = normalized.nnz_per_obs();
    ctx.obs_libsize = normalized.libsize_values().to_vec();
    ctx.obs_nnz = nnz.clone();
    ctx.obs_expressed_genes = nnz;

    ctx.normalized = Some(Box::new(normalized));
    ctx.normalization = Some(serde_json::json!({ "scale": scale }));

    Ok(())
}

fn load_from_pipeline_cache(
    ctx: &mut Ctx,
    cache_path: &std::path::Path,
    scale: f32,
) -> Result<(), StageError> {
    let normalized =
        cache_bin::MmapCscNormalized::new(cache_bin::SharedCache::open(cache_path)?, scale)?;

    if normalized.n_obs() != ctx.n_obs {
        return Err(StageError::Validation(
            "cache n_obs does not match Stage 0 obs_ids".to_string(),
        ));
    }
    if normalized.genes().len() != ctx.n_vars {
        return Err(StageError::Validation(
            "cache n_vars does not match Stage 0 gene symbols".to_string(),
        ));
    }
    if normalized.barcodes() != ctx.obs_ids {
        return Err(StageError::Validation(
            "cache barcodes do not match Stage 0 obs_ids".to_string(),
        ));
    }
    if !genes_compatible(normalized.genes(), &ctx.gene_symbols) {
        return Err(StageError::Validation(
            "cache genes do not match Stage 0 gene symbols".to_string(),
        ));
    }

    let zero_lib = normalized.zero_libsize_count();
    if zero_lib > 0 {
        warn!("{zero_lib} observations have zero libsize");
    }
    let nnz = normalized.nnz_per_obs()?;
    ctx.obs_libsize = normalized.libsize_values().to_vec();
    ctx.obs_nnz = nnz.clone();
    ctx.obs_expressed_genes = nnz;
    ctx.normalized = Some(Box::new(normalized));
    ctx.normalization = Some(serde_json::json!({ "scale": scale }));
    Ok(())
}

fn genes_compatible(cache_genes: &[String], stage0_genes: &[String]) -> bool {
    if cache_genes.len() != stage0_genes.len() {
        return false;
    }
    for (left, right) in cache_genes.iter().zip(stage0_genes.iter()) {
        if canonicalize_symbol(left) != canonicalize_symbol(right) {
            return false;
        }
    }
    true
}

fn canonicalize_symbol(symbol: &str) -> String {
    symbol.trim().to_ascii_uppercase()
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage1_normalization.rs"]
mod tests;
