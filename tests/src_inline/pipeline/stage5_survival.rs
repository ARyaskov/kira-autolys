
use super::*;
use crate::math::simd::EPS;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{AutophagyMetrics, Ctx, LysosomeMetrics, RegulatoryMetrics};

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

fn dummy_reg(n_obs: usize) -> RegulatoryMetrics {
    RegulatoryMetrics {
        tfeb: vec![0.0; n_obs],
        mtor: vec![0.0; n_obs],
        tfeb_act: vec![0.0; n_obs],
        mtor_supp: vec![0.0; n_obs],
        tfeb_mtor_diff: vec![0.0; n_obs],
        tfeb_mtor_ratio: vec![0.0; n_obs],
        median_tfeb: 0.0,
        median_mtor: 0.0,
    }
}

#[test]
fn test_survival_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.prolif, bits.apop_pro, bits.apop_anti];
    let normalized = DenseRows {
        rows: vec![
            vec![(0, 2.0), (1, 3.0), (2, 1.0)],
            vec![(0, 1.0), (1, 2.0), (2, 4.0)],
        ],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![1.0, 3.0],
        initiation: vec![0.0; 2],
        elongation: vec![0.0; 2],
        degradation: vec![0.0; 2],
        cargo: vec![0.0; 2],
        stall: vec![0.0; 2],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![2.0, 4.0],
        vatp: vec![0.0; 2],
        prot: vec![0.0; 2],
        mem: vec![0.0; 2],
        global_load: vec![0.0; 2],
    });
    ctx.regulatory = Some(dummy_reg(2));

    run_stage5(&mut ctx).unwrap();
    let metrics = ctx.survival.as_ref().unwrap();
    let stats = ctx.survival_stats.unwrap();

    assert_eq!(metrics.prolif, vec![2.0, 1.0]);
    assert_eq!(metrics.apop_ready, vec![2.0, -2.0]);

    assert!((stats.mean_afp - 2.0).abs() < 1e-6);
    assert!((stats.mean_lds - 3.0).abs() < 1e-6);

    let z_afp0 = (1.0 - stats.mean_afp) / (stats.std_afp + EPS);
    let z_lds0 = (2.0 - stats.mean_lds) / (stats.std_lds + EPS);
    let z_prolif0 = (2.0 - stats.mean_prolif) / (stats.std_prolif + EPS);
    let z_apop0 = (2.0 - stats.mean_apop) / (stats.std_apop + EPS);
    let ssm0 = 0.35 * z_afp0 + 0.35 * z_lds0 - 0.15 * z_prolif0 - 0.15 * z_apop0;
    assert!((metrics.ssm[0] - ssm0).abs() < 1e-6);
}

#[test]
fn test_flags_and_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.prolif];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 10.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["x".to_string(), "y".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![0.0, 10.0],
        initiation: vec![0.0; 2],
        elongation: vec![0.0; 2],
        degradation: vec![0.0; 2],
        cargo: vec![0.0; 2],
        stall: vec![0.0; 2],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0, 10.0],
        vatp: vec![0.0; 2],
        prot: vec![0.0; 2],
        mem: vec![0.0; 2],
        global_load: vec![0.0; 2],
    });
    ctx.regulatory = Some(dummy_reg(2));

    run_stage5(&mut ctx).unwrap();
    let metrics = ctx.survival.as_ref().unwrap();
    assert_eq!(metrics.prolif, vec![1.0, 10.0]);
    assert_eq!(metrics.flag_afp_high.len(), 2);
    assert_eq!(metrics.flag_lds_high.len(), 2);
    assert_eq!(metrics.flag_ssm_high.len(), 2);
}

#[cfg(feature = "parallel")]
#[test]
fn test_parallel_matches_sequential() {
    let n = 16;
    let values: Vec<f32> = (0..n).map(|v| v as f32 * 0.25).collect();
    let mut z_seq = vec![0.0_f32; n];
    let mut z_par = vec![0.0_f32; n];
    let mean = 1.5;
    let std = 0.5;

    compute_z_scores(&values, mean, std, &mut z_seq);
    compute_z_scores_parallel(&values, mean, std, &mut z_par);

    for (a, b) in z_seq.iter().zip(z_par.iter()) {
        assert!((a - b).abs() < 1e-6);
    }

    let thresholds = Thresholds::default();
    let mut ssm_seq = vec![0.0_f32; n];
    let mut ssm_par = vec![0.0_f32; n];
    let mut ssm_high_seq = vec![false; n];
    let mut afp_high_seq = vec![false; n];
    let mut lds_high_seq = vec![false; n];
    let mut prolif_low_seq = vec![false; n];
    let mut apop_low_seq = vec![false; n];
    let mut flags_seq = SurvivalFlags {
        ssm_high: &mut ssm_high_seq,
        afp_high: &mut afp_high_seq,
        lds_high: &mut lds_high_seq,
        prolif_low: &mut prolif_low_seq,
        apop_low: &mut apop_low_seq,
    };

    let mut ssm_high_par = vec![false; n];
    let mut afp_high_par = vec![false; n];
    let mut lds_high_par = vec![false; n];
    let mut prolif_low_par = vec![false; n];
    let mut apop_low_par = vec![false; n];
    let mut flags_par = SurvivalFlags {
        ssm_high: &mut ssm_high_par,
        afp_high: &mut afp_high_par,
        lds_high: &mut lds_high_par,
        prolif_low: &mut prolif_low_par,
        apop_low: &mut apop_low_par,
    };

    compute_ssm_and_flags(
        &z_seq,
        &z_seq,
        &z_seq,
        &z_seq,
        &thresholds,
        &mut ssm_seq,
        &mut flags_seq,
    );
    compute_ssm_and_flags_parallel(
        &z_par,
        &z_par,
        &z_par,
        &z_par,
        &thresholds,
        &mut ssm_par,
        &mut flags_par,
    );

    for (a, b) in ssm_seq.iter().zip(ssm_par.iter()) {
        assert!((a - b).abs() < 1e-6);
    }
}
