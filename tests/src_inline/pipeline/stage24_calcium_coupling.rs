
use super::*;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{AutophagyMetrics, CrossOrganelleMetrics, Ctx, LysosomeMetrics};

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

fn base_ctx(n_obs: usize) -> Ctx {
    let mut ctx = Ctx::default();
    ctx.n_obs = n_obs;
    ctx.obs_ids = (0..n_obs).map(|i| format!("obs{i}")).collect();
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![0.0; n_obs],
        initiation: vec![0.0; n_obs],
        elongation: vec![0.0; n_obs],
        degradation: vec![0.0; n_obs],
        cargo: vec![0.0; n_obs],
        stall: vec![0.0; n_obs],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0; n_obs],
        vatp: vec![0.0; n_obs],
        prot: vec![0.0; n_obs],
        mem: vec![0.0; n_obs],
        global_load: vec![0.0; n_obs],
    });
    ctx.cross_organelle = Some(CrossOrganelleMetrics {
        mitophagy_reliance: vec![0.0; n_obs],
        mito_lyso_imbalance: vec![0.0; n_obs],
        mito_consumption_risk: vec![0.0; n_obs],
        energy_recycling_dependency: vec![0.0; n_obs],
    });
    ctx
}

#[test]
fn test_coupling_aggregation() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.lcrc, bits.mcuc, bits.csap];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 4.0), (2, 6.0)]],
    };

    let mut ctx = base_ctx(1);
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage24(&mut ctx).unwrap();
    let metrics = ctx.calcium_coupling.as_ref().unwrap();
    assert!((metrics.lcrc[0] - 2.0).abs() < 1e-6);
    assert!((metrics.mcuc[0] - 4.0).abs() < 1e-6);
    assert!((metrics.csap[0] - 6.0).abs() < 1e-6);
}

#[test]
fn test_ccb_and_lmcci_finite() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.lcrc];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)]],
    };

    let mut ctx = base_ctx(1);
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage24(&mut ctx).unwrap();
    let metrics = ctx.calcium_coupling.as_ref().unwrap();
    assert!(metrics.ccb[0].is_finite());
    assert!(metrics.lmcci[0].is_finite());
}
