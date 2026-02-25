
use super::*;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{Ctx, LysosomeMetrics, RegulatoryMetrics};

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
fn test_coupling_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.ampk, bits.ragulator, 0];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 4.0), (2, 1.0)], vec![(0, 1.0), (1, 3.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["obs1".to_string(), "obs2".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.regulatory = Some(RegulatoryMetrics {
        tfeb: vec![1.0, 2.0],
        mtor: vec![2.0, 1.0],
        tfeb_act: vec![0.0; 2],
        mtor_supp: vec![0.0; 2],
        tfeb_mtor_diff: vec![0.0; 2],
        tfeb_mtor_ratio: vec![0.0; 2],
        median_tfeb: 0.0,
        median_mtor: 0.0,
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![1.0, 2.0],
        vatp: vec![0.0; 2],
        prot: vec![0.0; 2],
        mem: vec![0.0; 2],
        global_load: vec![0.0; 2],
    });

    run_stage9(&mut ctx).unwrap();
    let metrics = ctx.coupling.as_ref().unwrap();

    let ampk0: f32 = 2.0;
    let rag0: f32 = 4.0;
    let ratio0 = ampk0 / (2.0_f32 + 1e-6_f32);
    let align0 = 1.0_f32 * 1.0_f32;
    let mismatch0: f32 = rag0 - 2.0_f32;
    let lsi0 = 0.35_f32 * ratio0 + 0.35_f32 * align0 + 0.30_f32 * mismatch0.max(0.0_f32);

    assert!((metrics.ampk[0] - ampk0).abs() < 1e-6);
    assert!((metrics.ragulator[0] - rag0).abs() < 1e-6);
    assert!((metrics.ampk_mtor_ratio[0] - ratio0).abs() < 1e-6);
    assert!((metrics.tfeb_lyso_alignment[0] - align0).abs() < 1e-6);
    assert!((metrics.rag_mtor_mismatch[0] - mismatch0).abs() < 1e-6);
    assert!((metrics.locked_survival_index[0] - lsi0).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.ampk];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.regulatory = Some(RegulatoryMetrics {
        tfeb: vec![0.0; 2],
        mtor: vec![1.0; 2],
        tfeb_act: vec![0.0; 2],
        mtor_supp: vec![0.0; 2],
        tfeb_mtor_diff: vec![0.0; 2],
        tfeb_mtor_ratio: vec![0.0; 2],
        median_tfeb: 0.0,
        median_mtor: 0.0,
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0; 2],
        vatp: vec![0.0; 2],
        prot: vec![0.0; 2],
        mem: vec![0.0; 2],
        global_load: vec![0.0; 2],
    });

    run_stage9(&mut ctx).unwrap();
    let metrics = ctx.coupling.as_ref().unwrap();

    assert_eq!(metrics.ampk, vec![1.0, 2.0]);
}
