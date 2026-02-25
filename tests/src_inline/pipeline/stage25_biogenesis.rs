
use super::*;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{
    AcidificationMetrics, Ctx, LysosomalDamageMetrics, LysosomeMetrics, RegulatoryMetrics,
};

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
fn test_biogenesis_aggregation() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.biogenesis, bits.maturation, 0u64];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 4.0)], vec![(0, 1.0), (1, 3.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.regulatory = Some(RegulatoryMetrics {
        tfeb: vec![0.0; 2],
        mtor: vec![0.0; 2],
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
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0; 2],
        stress: vec![0.0; 2],
        cathepsin_membrane_imbalance: vec![0.0; 2],
        ldi: vec![0.5, 1.5],
    });

    run_stage25(&mut ctx).unwrap();
    let metrics = ctx.biogenesis.as_ref().unwrap();
    assert!((metrics.biogenesis_drive[0] - 2.0).abs() < 1e-6);
    assert!((metrics.maturation_yield[0] - 4.0).abs() < 1e-6);
}

#[test]
fn test_functional_capacity_with_lasi() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.biogenesis];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["a".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.regulatory = Some(RegulatoryMetrics {
        tfeb: vec![0.0],
        mtor: vec![0.0],
        tfeb_act: vec![0.0],
        mtor_supp: vec![0.0],
        tfeb_mtor_diff: vec![0.0],
        tfeb_mtor_ratio: vec![0.0],
        median_tfeb: 0.0,
        median_mtor: 0.0,
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![2.0],
        vatp: vec![0.0],
        prot: vec![0.0],
        mem: vec![0.0],
        global_load: vec![0.0],
    });
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0],
        stress: vec![0.0],
        cathepsin_membrane_imbalance: vec![0.0],
        ldi: vec![1.0],
    });
    ctx.acidification = Some(AcidificationMetrics {
        acap: vec![0.0],
        pap: vec![0.0],
        acl: vec![0.0],
        air: vec![0.0],
        lasi: vec![2.0],
    });

    run_stage25(&mut ctx).unwrap();
    let metrics = ctx.biogenesis.as_ref().unwrap();
    assert!(metrics.functional_capacity[0].is_finite());
    assert!(metrics.bpr[0].is_finite());
    assert!(metrics.lbpi[0].is_finite());
}
