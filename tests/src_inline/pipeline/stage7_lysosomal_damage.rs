use super::*;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{Ctx, LysosomeMetrics};

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
fn test_damage_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.lmp, bits.stress, 0];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 4.0), (2, 1.0)], vec![(0, 1.0), (1, 3.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["obs1".to_string(), "obs2".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0, 0.0],
        vatp: vec![0.0, 0.0],
        prot: vec![2.0, 1.0],
        mem: vec![1.0, 2.0],
        global_load: vec![0.0, 0.0],
    });

    run_stage7(&mut ctx).unwrap();
    let metrics = ctx.lysosomal_damage.as_ref().unwrap();

    let lmp0 = 2.0;
    let stress0 = 4.0;
    let imbalance0 = 2.0 / (1.0 + 1e-6);
    let ldi0 = 0.4 * lmp0 + 0.3 * stress0 + 0.3 * imbalance0;

    assert!((metrics.lmp[0] - lmp0).abs() < 1e-6);
    assert!((metrics.stress[0] - stress0).abs() < 1e-6);
    assert!((metrics.cathepsin_membrane_imbalance[0] - imbalance0).abs() < 1e-6);
    assert!((metrics.ldi[0] - ldi0).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.lmp];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0, 0.0],
        vatp: vec![0.0, 0.0],
        prot: vec![0.0, 0.0],
        mem: vec![1.0, 1.0],
        global_load: vec![0.0, 0.0],
    });

    run_stage7(&mut ctx).unwrap();
    let metrics = ctx.lysosomal_damage.as_ref().unwrap();

    assert_eq!(metrics.lmp, vec![1.0, 2.0]);
}
