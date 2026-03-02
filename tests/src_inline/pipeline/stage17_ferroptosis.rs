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
fn test_ferroptosis_metrics() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![
        bits.ferritin,
        bits.defense,
        bits.lipid,
        bits.iron,
        bits.iron,
    ];
    let normalized = DenseRows {
        rows: vec![vec![(0, 2.0), (1, 1.0), (2, 3.0), (3, 4.0), (4, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["a".to_string()];
    ctx.gene_symbols = vec![
        "NCOA4".to_string(),
        "GPX4".to_string(),
        "ACSL4".to_string(),
        "TFRC".to_string(),
        "SLC40A1".to_string(),
    ];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage17(&mut ctx).unwrap();
    let metrics = ctx.ferroptosis.as_ref().unwrap();

    let ferritin = 2.0;
    let defense = 1.0;
    let lipid = 3.0;
    let tfrc = 4.0;
    let slc40a1 = 2.0;
    let iron_bias = tfrc / (slc40a1 + EPS);
    let pressure = (ferritin * iron_bias * lipid) / (defense + EPS);

    assert!((metrics.ferritinophagy_load[0] - ferritin).abs() < 1e-6);
    assert!((metrics.iron_import_bias[0] - iron_bias).abs() < 1e-6);
    assert!((metrics.ferroptotic_pressure_index[0] - pressure).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let bits = panel_bits().unwrap();
    let gene_panel_mask = vec![bits.ferritin];
    let normalized = DenseRows {
        rows: vec![vec![(0, 1.0)], vec![(0, 2.0)]],
    };

    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.gene_symbols = vec!["NCOA4".to_string()];
    ctx.gene_panel_mask = gene_panel_mask;
    ctx.normalized = Some(Box::new(normalized));

    run_stage17(&mut ctx).unwrap();
    let metrics = ctx.ferroptosis.as_ref().unwrap();
    assert_eq!(metrics.ferritinophagy_load, vec![1.0, 2.0]);
}
