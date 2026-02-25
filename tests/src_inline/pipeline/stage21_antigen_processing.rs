
use super::*;
use crate::model::ctx::Ctx;
use crate::model::ctx::NormalizedExpr;

struct DenseRows {
    rows: Vec<Vec<(usize, f32)>>,
}

impl NormalizedExpr for DenseRows {
    fn n_obs(&self) -> usize {
        self.rows.len()
    }

    fn for_each_in_obs(&self, obs_idx: usize, f: &mut dyn FnMut(usize, f32)) {
        for &(idx, value) in &self.rows[obs_idx] {
            f(idx, value);
        }
    }
}

#[test]
fn test_antigen_processing_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.mhc, bits.protease, bits.accessory];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 3.0), (2, 4.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["a".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage21(&mut ctx).unwrap();
    let metrics = ctx.antigen_processing.as_ref().unwrap();

    assert!((metrics.mhc_expression_load[0] - 2.0).abs() < 1e-6);
    assert!((metrics.protease_load[0] - 3.0).abs() < 1e-6);
    assert!((metrics.loading_accessory_load[0] - 4.0).abs() < 1e-6);
    assert!((metrics.antigen_processing_index[0] - 24.0).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.mhc];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage21(&mut ctx).unwrap();
    let metrics = ctx.antigen_processing.as_ref().unwrap();
    assert_eq!(metrics.mhc_expression_load, vec![1.0, 2.0]);
}
