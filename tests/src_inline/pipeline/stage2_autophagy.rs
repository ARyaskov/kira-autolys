use super::*;
use crate::model::ctx::Ctx;

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
fn test_autophagy_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.init, bits.elong, bits.late, bits.cargo];
    let normalized = DenseRows {
        rows: vec![
            vec![(0, 2.0), (1, 4.0), (2, 1.0), (3, 3.0)],
            vec![(0, 1.0), (2, 2.0)],
        ],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["obs1".to_string(), "obs2".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage2(&mut ctx).unwrap();
    let metrics = ctx.autophagy.as_ref().unwrap();
    let extension = ctx.autolys_extension.as_ref().unwrap();

    assert_eq!(metrics.initiation.len(), 2);
    assert_eq!(metrics.elongation.len(), 2);
    assert_eq!(metrics.degradation.len(), 2);
    assert_eq!(metrics.cargo.len(), 2);
    assert_eq!(metrics.stall.len(), 2);

    let init0 = 2.0_f32;
    let elong0 = 4.0_f32;
    let late0 = 1.0_f32;
    let cargo0 = 3.0_f32;
    let axis_early0 = 0.5_f32 * (init0 + elong0);
    let axis_late0 = late0;
    let stall0: f32 = (axis_early0 - axis_late0).max(0.0_f32);
    let afp_raw0 = 0.5_f32 * axis_early0 + 0.3_f32 * axis_late0 + 0.2_f32 * cargo0;
    let afp0 = afp_raw0 - 0.25_f32 * stall0;

    assert!((metrics.initiation[0] - init0).abs() < 1e-6);
    assert!((metrics.elongation[0] - elong0).abs() < 1e-6);
    assert!((metrics.degradation[0] - late0).abs() < 1e-6);
    assert!((metrics.cargo[0] - cargo0).abs() < 1e-6);
    assert!((metrics.stall[0] - stall0).abs() < 1e-6);
    assert!((metrics.afp[0] - afp0).abs() < 1e-6);

    assert_eq!(metrics.initiation[1], 1.0);
    assert_eq!(metrics.elongation[1], 0.0);
    assert_eq!(metrics.degradation[1], 2.0);
    assert_eq!(metrics.cargo[1], 0.0);
    assert_eq!(metrics.stall[1], 0.0);
    assert_eq!(extension.ais.len(), 2);
    assert_eq!(extension.asm.len(), 2);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.init];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage2(&mut ctx).unwrap();
    let metrics = ctx.autophagy.as_ref().unwrap();

    assert_eq!(metrics.initiation, vec![1.0, 2.0]);
}

#[cfg(feature = "parallel")]
#[test]
fn test_parallel_matches_sequential() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.init, bits.elong, bits.late, bits.cargo];
    let normalized = DenseRows {
        rows: vec![
            vec![(0, 2.0), (1, 4.0), (2, 1.0), (3, 3.0)],
            vec![(0, 1.0), (2, 2.0)],
            vec![(1, 1.5), (3, 2.5)],
            vec![(0, 0.5), (1, 0.5), (2, 0.5), (3, 0.5)],
        ],
    };

    let mut seq = AutophagyMetrics {
        afp: vec![0.0; normalized.n_obs()],
        initiation: vec![0.0; normalized.n_obs()],
        elongation: vec![0.0; normalized.n_obs()],
        degradation: vec![0.0; normalized.n_obs()],
        cargo: vec![0.0; normalized.n_obs()],
        stall: vec![0.0; normalized.n_obs()],
    };
    stage2_compute_block(
        &normalized,
        &gene_panel_mask,
        bits.init,
        bits.elong,
        bits.late,
        bits.cargo,
        &mut seq,
        0..normalized.n_obs(),
    )
    .unwrap();

    let mut par = AutophagyMetrics {
        afp: vec![0.0; normalized.n_obs()],
        initiation: vec![0.0; normalized.n_obs()],
        elongation: vec![0.0; normalized.n_obs()],
        degradation: vec![0.0; normalized.n_obs()],
        cargo: vec![0.0; normalized.n_obs()],
        stall: vec![0.0; normalized.n_obs()],
    };
    stage2_compute_parallel(
        &normalized,
        &gene_panel_mask,
        bits.init,
        bits.elong,
        bits.late,
        bits.cargo,
        &mut par,
    )
    .unwrap();

    for (a, b) in seq.afp.iter().zip(par.afp.iter()) {
        assert!((a - b).abs() < 1e-6);
    }
}
