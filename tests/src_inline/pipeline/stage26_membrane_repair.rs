
use super::*;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{AutophagyMetrics, Ctx, LysosomalDamageMetrics, LysosomeMetrics};

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
    ctx
}

#[test]
fn test_repair_aggregation() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.escrt, bits.ca, bits.lysophagy, bits.stabilization];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 4.0), (2, 6.0), (3, 8.0)]],
    };

    let mut ctx = base_ctx(1);
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0],
        stress: vec![0.0],
        cathepsin_membrane_imbalance: vec![0.0],
        ldi: vec![1.0],
    });

    run_stage26(&mut ctx).unwrap();
    let metrics = ctx.membrane_repair.as_ref().unwrap();
    assert!((metrics.rma[0] - 3.0).abs() < 1e-6);
    assert!((metrics.lpe[0] - 6.0).abs() < 1e-6);
    assert!((metrics.msr[0] - 8.0).abs() < 1e-6);
    assert!((metrics.erc[0] - 17.0).abs() < 1e-6);
}

#[test]
fn test_rdr_and_lmrci_finite() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.escrt];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)]],
    };

    let mut ctx = base_ctx(1);
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0],
        stress: vec![0.0],
        cathepsin_membrane_imbalance: vec![0.0],
        ldi: vec![0.0],
    });

    run_stage26(&mut ctx).unwrap();
    let metrics = ctx.membrane_repair.as_ref().unwrap();
    assert!(metrics.rdr[0].is_finite());
    assert!(metrics.lmrci[0].is_finite());
}
