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
fn test_positioning_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.perinuclear, bits.peripheral];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 1.0)], vec![(0, 1.0), (1, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage16(&mut ctx).unwrap();
    let metrics = ctx.lysosome_positioning.as_ref().unwrap();

    let ratio0 = (2.0 + EPS) / (1.0 + EPS);
    assert!((metrics.positioning_ratio[0] - ratio0).abs() < 1e-6);
    assert!(metrics.positioning_bias[0] > 0.0);

    let ratio1 = (1.0 + EPS) / (2.0 + EPS);
    assert!((metrics.positioning_ratio[1] - ratio1).abs() < 1e-6);
    assert!(metrics.positioning_bias[1] < 0.0);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.perinuclear];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage16(&mut ctx).unwrap();
    let metrics = ctx.lysosome_positioning.as_ref().unwrap();

    assert_eq!(metrics.perinuclear_mean, vec![1.0, 2.0]);
}
