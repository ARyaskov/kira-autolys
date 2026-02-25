use std::path::Path;
use std::sync::Arc;

use crate::expr::normalized::NormalizedExpr;
use crate::stage_error::StageError;

#[derive(Debug)]
pub struct SharedCache {
    inner: Arc<kira_shared_sc_cache::SharedCacheMmap>,
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
}

impl SharedCache {
    pub fn open(path: &Path) -> Result<Self, StageError> {
        let inner = kira_shared_sc_cache::mmap_shared_cache(path).map_err(map_err)?;
        Ok(Self {
            n_genes: inner.n_genes,
            n_cells: inner.n_cells,
            nnz: inner.nnz,
            genes: inner.genes.clone(),
            barcodes: inner.barcodes.clone(),
            inner: Arc::new(inner),
        })
    }

    pub fn col_ptr(&self) -> Result<&[u64], StageError> {
        Ok(self.inner.col_ptr())
    }

    pub fn row_idx(&self) -> Result<&[u32], StageError> {
        Ok(self.inner.row_idx())
    }

    pub fn values_u32(&self) -> Result<&[u32], StageError> {
        Ok(self.inner.values_u32())
    }
}

pub struct MmapCscNormalized {
    cache: Arc<SharedCache>,
    libsize: Vec<f32>,
    scale: f32,
}

impl MmapCscNormalized {
    pub fn new(cache: SharedCache, scale: f32) -> Result<Self, StageError> {
        let cache = Arc::new(cache);
        let col_ptr = cache.col_ptr()?;
        let values = cache.values_u32()?;

        let mut libsize = vec![0.0_f32; cache.n_cells];
        for obs in 0..cache.n_cells {
            let start = usize_from_u64(col_ptr[obs])?;
            let end = usize_from_u64(col_ptr[obs + 1])?;
            let mut sum = 0.0_f32;
            for &value in &values[start..end] {
                sum += value as f32;
            }
            libsize[obs] = sum;
        }

        Ok(Self {
            cache,
            libsize,
            scale,
        })
    }

    pub fn zero_libsize_count(&self) -> usize {
        self.libsize.iter().filter(|v| **v == 0.0).count()
    }

    pub fn genes(&self) -> &[String] {
        &self.cache.genes
    }

    pub fn barcodes(&self) -> &[String] {
        &self.cache.barcodes
    }

    pub fn libsize_values(&self) -> &[f32] {
        &self.libsize
    }

    pub fn nnz_per_obs(&self) -> Result<Vec<u32>, StageError> {
        let col_ptr = self.cache.col_ptr()?;
        let mut out = Vec::with_capacity(self.cache.n_cells);
        for obs in 0..self.cache.n_cells {
            let start = usize_from_u64(col_ptr[obs])?;
            let end = usize_from_u64(col_ptr[obs + 1])?;
            out.push((end - start) as u32);
        }
        Ok(out)
    }
}

impl NormalizedExpr for MmapCscNormalized {
    fn n_obs(&self) -> usize {
        self.cache.n_cells
    }

    fn for_each_in_obs(&self, obs_idx: usize, f: &mut dyn FnMut(usize, f32)) {
        let Ok(col_ptr) = self.cache.col_ptr() else {
            return;
        };
        let Ok(row_idx) = self.cache.row_idx() else {
            return;
        };
        let Ok(values_u32) = self.cache.values_u32() else {
            return;
        };

        let start = col_ptr[obs_idx] as usize;
        let end = col_ptr[obs_idx + 1] as usize;
        let libsize = self.libsize[obs_idx];
        if libsize == 0.0 {
            return;
        }

        let scale = self.scale;
        for i in start..end {
            let gene_idx = row_idx[i] as usize;
            let count = values_u32[i] as f32;
            let norm = ((count / libsize) * scale).ln_1p();
            f(gene_idx, norm);
        }
    }
}

fn usize_from_u64(value: u64) -> Result<usize, StageError> {
    usize::try_from(value).map_err(|_| StageError::Format("value does not fit usize".to_string()))
}

fn map_err(err: kira_shared_sc_cache::SharedCacheError) -> StageError {
    match err {
        kira_shared_sc_cache::SharedCacheError::Io { source, .. } => StageError::Io(source),
        kira_shared_sc_cache::SharedCacheError::Format { message, .. } => {
            StageError::Format(message)
        }
    }
}

pub fn crc64_ecma(bytes: &[u8]) -> u64 {
    kira_shared_sc_cache::crc64_ecma(bytes)
}

#[cfg(test)]
#[path = "../../tests/src_inline/io/cache_bin.rs"]
mod tests;
