
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
fn test_lipid_buffering_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.storage, bits.utilization, bits.context];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 1.0), (2, 4.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["a".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage22(&mut ctx).unwrap();
    let metrics = ctx.lipid_buffering.as_ref().unwrap();

    let storage = 2.0;
    let utilization = 1.0;
    let context = 4.0;
    let index = (storage * context) / (utilization + EPS);

    assert!((metrics.storage_buffering_load[0] - storage).abs() < 1e-6);
    assert!((metrics.utilization_load[0] - utilization).abs() < 1e-6);
    assert!((metrics.lipophagy_context[0] - context).abs() < 1e-6);
    assert!((metrics.lipid_buffering_index[0] - index).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.storage];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage22(&mut ctx).unwrap();
    let metrics = ctx.lipid_buffering.as_ref().unwrap();
    assert_eq!(metrics.storage_buffering_load, vec![1.0, 2.0]);
}
