use std::sync::Arc;

use crate::stage_error::StageError;

pub trait NormalizedExpr: Send + Sync {
    fn n_obs(&self) -> usize;
    fn for_each_in_obs(&self, obs_idx: usize, f: &mut dyn FnMut(usize, f32));
}

#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n_vars: usize,
    pub n_obs: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<f32>,
}

impl CscMatrix {
    pub fn validate(&self) -> Result<(), StageError> {
        if self.indptr.len() != self.n_obs + 1 {
            return Err(StageError::Validation(
                "CSC indptr length does not match n_obs".to_string(),
            ));
        }
        if self.indices.len() != self.data.len() {
            return Err(StageError::Validation(
                "CSC indices/data length mismatch".to_string(),
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct CscNormalized {
    matrix: Arc<CscMatrix>,
    libsize: Vec<f32>,
    scale: f32,
}

impl CscNormalized {
    pub fn new(matrix: CscMatrix, scale: f32) -> Result<Self, StageError> {
        matrix.validate()?;
        let libsize = compute_libsize(&matrix);
        Ok(Self {
            matrix: Arc::new(matrix),
            libsize,
            scale,
        })
    }

    pub fn zero_libsize_count(&self) -> usize {
        self.libsize.iter().filter(|v| **v == 0.0).count()
    }

    pub fn libsize_values(&self) -> &[f32] {
        &self.libsize
    }

    pub fn nnz_per_obs(&self) -> Vec<u32> {
        let mut out = Vec::with_capacity(self.matrix.n_obs);
        for obs_idx in 0..self.matrix.n_obs {
            let start = self.matrix.indptr[obs_idx];
            let end = self.matrix.indptr[obs_idx + 1];
            out.push((end - start) as u32);
        }
        out
    }
}

impl NormalizedExpr for CscNormalized {
    fn n_obs(&self) -> usize {
        self.matrix.n_obs
    }

    fn for_each_in_obs(&self, obs_idx: usize, f: &mut dyn FnMut(usize, f32)) {
        let start = self.matrix.indptr[obs_idx];
        let end = self.matrix.indptr[obs_idx + 1];
        let libsize = self.libsize[obs_idx];
        if libsize == 0.0 {
            return;
        }
        let scale = self.scale;
        for i in start..end {
            let gene_idx = self.matrix.indices[i];
            let count = self.matrix.data[i];
            let norm = ((count / libsize) * scale).ln_1p();
            f(gene_idx, norm);
        }
    }
}

pub fn compute_libsize(matrix: &CscMatrix) -> Vec<f32> {
    let mut libsize = vec![0.0; matrix.n_obs];
    for obs_idx in 0..matrix.n_obs {
        let start = matrix.indptr[obs_idx];
        let end = matrix.indptr[obs_idx + 1];
        let mut sum = 0.0_f32;
        for i in start..end {
            sum += matrix.data[i];
        }
        libsize[obs_idx] = sum;
    }
    libsize
}
