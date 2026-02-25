use std::collections::BTreeMap;
use std::path::Path;
use tracing::warn;

use crate::genesets::v1;
use crate::io::{cache_bin, dense_tsv, mtx};
use crate::model::cli::CliArgs;
use crate::model::cli::RunMode;
use crate::model::ctx::{Ctx, GeneSetResolved, InputType, MissingGenes, Mode};
use crate::stage_error::StageError;

pub fn run_stage0(ctx: &mut Ctx, cli: &CliArgs) -> Result<(), StageError> {
    let input = &cli.input;
    ctx.run_mode = cli.run_mode;

    if matches!(cli.run_mode, RunMode::Pipeline) {
        let cache_candidate = if let Some(path) = cli.cache_path.clone() {
            Some(path)
        } else if input.is_dir() {
            Some(mtx::resolve_cache_bin_path(input)?)
        } else {
            None
        };
        if let Some(cache_path) = cache_candidate {
            if cache_path.is_file() {
                let cache = cache_bin::SharedCache::open(&cache_path)?;
                return finish_stage0(ctx, cli, InputType::Mtx10x, cache.barcodes, cache.genes);
            }
            if input.is_dir() {
                // Preserve directory validation behavior, but fail with explicit pipeline guidance.
                let _ = mtx::detect_mtx_dir(input)?;
                return Err(StageError::Validation(format!(
                    "pipeline mode expects shared cache file '{}' (produced by kira-mitoqc)",
                    cache_path.display()
                )));
            }
        }
    }

    let (input_type, obs_ids, gene_symbols_raw) = if input.is_dir() {
        if dense_tsv::resolve_dense_input_path(input).is_ok() {
            let (obs_ids, gene_symbols) = dense_tsv::load_dense_metadata(input)?;
            (InputType::DenseTsv, obs_ids, gene_symbols)
        } else {
            let paths = mtx::detect_mtx_dir(input)?;
            let obs_ids = mtx::read_barcodes(&paths.barcodes)?;
            let gene_symbols = mtx::read_features(&paths.features)?;
            (InputType::Mtx10x, obs_ids, gene_symbols)
        }
    } else if input.is_file() {
        if is_h5ad(input) {
            #[cfg(feature = "h5ad")]
            {
                let (obs_ids, gene_symbols) = crate::io::h5ad::read_obs_var(input)?;
                (InputType::H5ad, obs_ids, gene_symbols)
            }
            #[cfg(not(feature = "h5ad"))]
            {
                return Err(StageError::FeatureDisabled(
                    "h5ad support requires the 'h5ad' feature".to_string(),
                ));
            }
        } else if dense_tsv::is_dense_tsv_path(input) {
            let (obs_ids, gene_symbols) = dense_tsv::load_dense_metadata(input)?;
            (InputType::DenseTsv, obs_ids, gene_symbols)
        } else {
            return Err(StageError::Validation(format!(
                "unsupported input file: {}",
                input.display()
            )));
        }
    } else {
        return Err(StageError::Validation(format!(
            "input path does not exist: {}",
            input.display()
        )));
    };

    finish_stage0(ctx, cli, input_type, obs_ids, gene_symbols_raw)
}

fn finish_stage0(
    ctx: &mut Ctx,
    cli: &CliArgs,
    input_type: InputType,
    obs_ids: Vec<String>,
    gene_symbols_raw: Vec<String>,
) -> Result<(), StageError> {
    let input = &cli.input;

    let gene_symbols = canonicalize_symbols(&gene_symbols_raw);
    let (gene_symbols, gene_index) = build_gene_index(&gene_symbols)?;

    let n_obs = obs_ids.len();
    let n_vars = gene_symbols.len();

    let mode = cli.mode.unwrap_or_else(|| auto_mode(input_type, n_obs));

    let (gene_sets, missing_genes, gene_panel_mask) = resolve_gene_sets(&gene_index, n_vars)?;

    ctx.input_type = Some(input_type);
    ctx.mode = Some(mode);
    ctx.input_path = Some(input.clone());
    ctx.obs_ids = obs_ids;
    ctx.gene_symbols = gene_symbols;
    ctx.gene_index = gene_index;
    ctx.n_obs = n_obs;
    ctx.n_vars = n_vars;
    ctx.gene_sets = gene_sets;
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.missing_genes = missing_genes;

    Ok(())
}

fn is_h5ad(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("h5ad"))
        .unwrap_or(false)
}

fn canonicalize_symbols(symbols: &[String]) -> Vec<String> {
    symbols.iter().map(|sym| canonicalize_symbol(sym)).collect()
}

fn canonicalize_symbol(symbol: &str) -> String {
    symbol.trim().to_ascii_uppercase()
}

fn build_gene_index(
    gene_symbols: &[String],
) -> Result<(Vec<String>, BTreeMap<String, usize>), StageError> {
    let mut index = BTreeMap::new();
    let symbols = gene_symbols.to_vec();

    for (idx, symbol) in symbols.iter().enumerate() {
        if symbol.is_empty() {
            warn!("empty gene symbol at index {}; retained for alignment", idx);
            continue;
        }
        if let Some(existing) = index.get(symbol) {
            warn!(
                "duplicate gene symbol {}; first occurrence index {} retained",
                symbol, existing
            );
            continue;
        }
        index.insert(symbol.clone(), idx);
    }

    apply_aliases(&mut index, &symbols)?;

    Ok((symbols, index))
}

fn apply_aliases(
    index: &mut BTreeMap<String, usize>,
    gene_symbols: &[String],
) -> Result<(), StageError> {
    for (left, right) in alias_pairs() {
        let left = canonicalize_symbol(left);
        let right = canonicalize_symbol(right);
        let left_idx = index.get(&left).copied();
        let right_idx = index.get(&right).copied();
        match (left_idx, right_idx) {
            (Some(l), None) => {
                index.insert(right, l);
            }
            (None, Some(r)) => {
                index.insert(left, r);
            }
            (Some(l), Some(r)) => {
                if l != r {
                    let min = l.min(r);
                    let max = l.max(r);
                    let symbol = gene_symbols
                        .get(max)
                        .map(|s| s.as_str())
                        .unwrap_or("<unknown>");
                    warn!(
                        "alias collision between {} and {}; index {} retained, {} at index {} remapped",
                        left, right, min, symbol, max
                    );
                    index.insert(left, min);
                    index.insert(right, min);
                }
            }
            (None, None) => {}
        }
    }
    Ok(())
}

fn alias_pairs() -> &'static [(&'static str, &'static str)] {
    &[("RB1CC1", "FIP200"), ("SQSTM1", "P62")]
}

fn auto_mode(input_type: InputType, n_obs: usize) -> Mode {
    match input_type {
        InputType::Mtx10x | InputType::DenseTsv if n_obs >= 1000 => Mode::Cell,
        _ => Mode::Sample,
    }
}

fn resolve_gene_sets(
    gene_index: &BTreeMap<String, usize>,
    n_vars: usize,
) -> Result<(Vec<GeneSetResolved>, Vec<MissingGenes>, Vec<u64>), StageError> {
    let defs = v1::gene_sets();
    if defs.len() > 64 {
        return Err(StageError::Validation(
            "gene set count exceeds 64 panels".to_string(),
        ));
    }

    let mut resolved = Vec::with_capacity(defs.len());
    let mut missing = Vec::with_capacity(defs.len());
    let mut panel_mask = vec![0u64; n_vars];

    for (panel_idx, def) in defs.iter().enumerate() {
        let mut indices = Vec::new();
        let mut missing_genes = Vec::new();
        for gene in def.genes {
            let canonical = canonicalize_symbol(gene);
            if let Some(&idx) = gene_index.get(&canonical) {
                indices.push(idx);
                panel_mask[idx] |= 1u64 << panel_idx;
            } else {
                missing_genes.push(canonical);
            }
        }
        missing_genes.sort();
        if !missing_genes.is_empty() {
            warn!(
                "panel {} missing genes: {}",
                def.name,
                missing_genes.join(", ")
            );
        }
        resolved.push(GeneSetResolved {
            name: def.name.to_string(),
            indices,
        });
        missing.push(MissingGenes {
            panel: def.name.to_string(),
            missing: missing_genes,
        });
    }

    Ok((resolved, missing, panel_mask))
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage0_input.rs"]
mod tests;
