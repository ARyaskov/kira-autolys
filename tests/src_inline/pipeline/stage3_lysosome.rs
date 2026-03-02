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
fn test_lysosome_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.vatp, bits.prot, bits.mem, 0];
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

    run_stage3(&mut ctx).unwrap();
    let metrics = ctx.lysosome.as_ref().unwrap();

    let vatp0 = 2.0;
    let prot0 = 4.0;
    let mem0 = 1.0;
    let global0 = (2.0 + 4.0 + 1.0 + 3.0) / 4.0;
    let lds_raw0 = 0.5 * vatp0 + 0.3 * prot0 + 0.2 * mem0;
    let lds0 = lds_raw0 / (global0 + 1e-6);

    assert!((metrics.vatp[0] - vatp0).abs() < 1e-6);
    assert!((metrics.prot[0] - prot0).abs() < 1e-6);
    assert!((metrics.mem[0] - mem0).abs() < 1e-6);
    assert!((metrics.global_load[0] - global0).abs() < 1e-6);
    assert!((metrics.lds[0] - lds0).abs() < 1e-6);

    assert_eq!(metrics.vatp[1], 1.0);
    assert_eq!(metrics.prot[1], 0.0);
    assert_eq!(metrics.mem[1], 2.0);
    assert_eq!(metrics.global_load[1], (1.0 + 2.0) / 2.0);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.vatp];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage3(&mut ctx).unwrap();
    let metrics = ctx.lysosome.as_ref().unwrap();

    assert_eq!(metrics.vatp, vec![1.0, 2.0]);
}

#[cfg(feature = "parallel")]
#[test]
fn test_parallel_matches_sequential() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.vatp, bits.prot, bits.mem, 0];
    let normalized = DenseRows {
        rows: vec![
            vec![(0, 2.0), (1, 4.0), (2, 1.0), (3, 3.0)],
            vec![(0, 1.0), (2, 2.0)],
            vec![(1, 1.5), (2, 2.5)],
            vec![(0, 0.5), (1, 0.5), (2, 0.5), (3, 0.5)],
        ],
    };

    let mut seq = LysosomeMetrics {
        lds: vec![0.0; normalized.n_obs()],
        vatp: vec![0.0; normalized.n_obs()],
        prot: vec![0.0; normalized.n_obs()],
        mem: vec![0.0; normalized.n_obs()],
        global_load: vec![0.0; normalized.n_obs()],
    };
    stage3_compute_block(
        &normalized,
        &gene_panel_mask,
        bits.vatp,
        bits.prot,
        bits.mem,
        &mut seq,
        0..normalized.n_obs(),
    )
    .unwrap();

    let mut par = LysosomeMetrics {
        lds: vec![0.0; normalized.n_obs()],
        vatp: vec![0.0; normalized.n_obs()],
        prot: vec![0.0; normalized.n_obs()],
        mem: vec![0.0; normalized.n_obs()],
        global_load: vec![0.0; normalized.n_obs()],
    };
    stage3_compute_parallel(
        &normalized,
        &gene_panel_mask,
        bits.vatp,
        bits.prot,
        bits.mem,
        &mut par,
    )
    .unwrap();

    for (a, b) in seq.lds.iter().zip(par.lds.iter()) {
        assert!((a - b).abs() < 1e-6);
    }
}
