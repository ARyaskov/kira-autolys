
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
fn test_ros_aggregation() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.ros, bits.iron, bits.anti, bits.damage];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 4.0), (2, 8.0), (3, 1.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["obs0".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage23(&mut ctx).unwrap();
    let metrics = ctx.lysosomal_ros.as_ref().unwrap();
    assert!((metrics.ros_generation_load[0] - 2.0).abs() < 1e-6);
    assert!((metrics.iron_redox_context[0] - 4.0).abs() < 1e-6);
    assert!((metrics.antioxidant_capacity[0] - 8.0).abs() < 1e-6);
    assert!((metrics.damage_response_load[0] - 1.0).abs() < 1e-6);
}

#[test]
fn test_ros_stress_index_eps() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.ros];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["obs0".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage23(&mut ctx).unwrap();
    let metrics = ctx.lysosomal_ros.as_ref().unwrap();
    assert!(metrics.lysosomal_ros_stress_index[0].is_finite());
}
