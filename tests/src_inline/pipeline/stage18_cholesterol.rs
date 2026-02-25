
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
fn test_cholesterol_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.export, bits.import, bits.efflux];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0), (1, 3.0), (2, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["a".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage18(&mut ctx).unwrap();
    let metrics = ctx.cholesterol.as_ref().unwrap();

    let export = 1.0;
    let import = 3.0;
    let efflux = 2.0;
    let trap = import / (export + efflux + EPS);

    assert!((metrics.export_capacity[0] - export).abs() < 1e-6);
    assert!((metrics.import_pressure[0] - import).abs() < 1e-6);
    assert!((metrics.efflux_capacity[0] - efflux).abs() < 1e-6);
    assert!((metrics.cholesterol_trap_index[0] - trap).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.export];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage18(&mut ctx).unwrap();
    let metrics = ctx.cholesterol.as_ref().unwrap();
    assert_eq!(metrics.export_capacity, vec![1.0, 2.0]);
}
