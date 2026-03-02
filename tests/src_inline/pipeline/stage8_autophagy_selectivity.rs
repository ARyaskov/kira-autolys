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
fn test_selectivity_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.mito, bits.aggre, bits.er, bits.ferr, bits.lipo];
    let normalized = DenseRows {
        rows: vec![
            vec![(0, 2.0), (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)],
            vec![(0, 1.0), (1, 3.0)],
        ],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["obs1".to_string(), "obs2".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage8(&mut ctx).unwrap();
    let metrics = ctx.autophagy_selectivity.as_ref().unwrap();

    let total0 = 2.0_f32 + 1.0_f32 + 1.0_f32 + 1.0_f32 + 1.0_f32 + 1e-6_f32;
    let mito_frac0 = 2.0_f32 / total0;
    let aggre_frac0 = 1.0_f32 / total0;
    let er_frac0 = 1.0_f32 / total0;
    let ferr_frac0 = 1.0_f32 / total0;
    let lipo_frac0 = 1.0_f32 / total0;
    let entropy0 = -(mito_frac0 * (mito_frac0 + 1e-6_f32).ln()
        + aggre_frac0 * (aggre_frac0 + 1e-6_f32).ln()
        + er_frac0 * (er_frac0 + 1e-6_f32).ln()
        + ferr_frac0 * (ferr_frac0 + 1e-6_f32).ln()
        + lipo_frac0 * (lipo_frac0 + 1e-6_f32).ln());

    assert!((metrics.mitophagy[0] - 2.0).abs() < 1e-6);
    assert!((metrics.aggrephagy[0] - 1.0).abs() < 1e-6);
    assert!((metrics.erphagy[0] - 1.0).abs() < 1e-6);
    assert!((metrics.ferritinophagy[0] - 1.0).abs() < 1e-6);
    assert!((metrics.lipophagy[0] - 1.0).abs() < 1e-6);
    assert!((metrics.mito_frac[0] - mito_frac0).abs() < 1e-6);
    assert!((metrics.aggre_frac[0] - aggre_frac0).abs() < 1e-6);
    assert!((metrics.er_frac[0] - er_frac0).abs() < 1e-6);
    assert!((metrics.ferr_frac[0] - ferr_frac0).abs() < 1e-6);
    assert!((metrics.lipo_frac[0] - lipo_frac0).abs() < 1e-6);
    assert!((metrics.entropy[0] - entropy0).abs() < 1e-6);
}

#[test]
fn test_fraction_sum_and_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.mito];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage8(&mut ctx).unwrap();
    let metrics = ctx.autophagy_selectivity.as_ref().unwrap();

    let sum0 = metrics.mito_frac[0]
        + metrics.aggre_frac[0]
        + metrics.er_frac[0]
        + metrics.ferr_frac[0]
        + metrics.lipo_frac[0];
    assert!((sum0 - 1.0).abs() < 1e-4);
    assert_eq!(metrics.mitophagy, vec![1.0, 2.0]);
}
