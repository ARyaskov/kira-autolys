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
fn test_regulatory_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.tfeb, bits.mtor];
    let normalized = DenseRows {
        rows: vec![
            vec![(0, 2.0), (1, 1.0)],
            vec![(0, 4.0), (1, 3.0)],
            vec![(0, 0.0), (1, 5.0)],
        ],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 3;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string(), "c".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage4(&mut ctx).unwrap();
    let metrics = ctx.regulatory.as_ref().unwrap();

    let tfeb = &metrics.tfeb;
    let mtor = &metrics.mtor;
    assert_eq!(tfeb, &vec![2.0, 4.0, 0.0]);
    assert_eq!(mtor, &vec![1.0, 3.0, 5.0]);

    assert!((metrics.median_tfeb - 2.0).abs() < 1e-6);
    assert!((metrics.median_mtor - 3.0).abs() < 1e-6);

    assert!((metrics.tfeb_act[0] - 0.0).abs() < 1e-6);
    assert!((metrics.tfeb_act[1] - 2.0).abs() < 1e-6);
    assert!((metrics.tfeb_act[2] + 2.0).abs() < 1e-6);

    assert!((metrics.mtor_supp[0] - 2.0).abs() < 1e-6);
    assert!((metrics.mtor_supp[1] - 0.0).abs() < 1e-6);
    assert!((metrics.mtor_supp[2] - 0.0).abs() < 1e-6);

    assert!((metrics.tfeb_mtor_diff[0] - 1.0).abs() < 1e-6);
    assert!((metrics.tfeb_mtor_diff[1] - 1.0).abs() < 1e-6);
    assert!((metrics.tfeb_mtor_diff[2] + 5.0).abs() < 1e-6);

    let ratio0 = (0.0_f32 + 1e-6_f32) / ((1.0_f32 - 3.0_f32).abs() + 1e-6_f32);
    let ratio1 = (2.0_f32 + 1e-6_f32) / ((3.0_f32 - 3.0_f32).abs() + 1e-6_f32);
    let ratio2 = (-2.0_f32 + 1e-6_f32) / ((5.0_f32 - 3.0_f32).abs() + 1e-6_f32);
    assert!((metrics.tfeb_mtor_ratio[0] - ratio0).abs() < 1e-6);
    assert!((metrics.tfeb_mtor_ratio[1] - ratio1).abs() < 1e-6);
    assert!((metrics.tfeb_mtor_ratio[2] - ratio2).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.tfeb];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage4(&mut ctx).unwrap();
    let metrics = ctx.regulatory.as_ref().unwrap();

    assert_eq!(metrics.tfeb, vec![1.0, 2.0]);
}
