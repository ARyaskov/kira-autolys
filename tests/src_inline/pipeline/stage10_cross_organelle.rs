use super::*;
use crate::model::ctx::{
    AutophagySelectivityMetrics, CouplingMetrics, Ctx, LysosomalDamageMetrics, LysosomeMetrics,
    MitoMetrics,
};

fn base_ctx() -> Ctx {
    let mut ctx = Ctx::default();
    ctx.n_obs = 2;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
    ctx.autophagy_selectivity = Some(AutophagySelectivityMetrics {
        mitophagy: vec![2.0, 1.0],
        aggrephagy: vec![0.0; 2],
        erphagy: vec![0.0; 2],
        ferritinophagy: vec![0.0; 2],
        lipophagy: vec![0.0; 2],
        mito_frac: vec![0.0; 2],
        aggre_frac: vec![0.0; 2],
        er_frac: vec![0.0; 2],
        ferr_frac: vec![0.0; 2],
        lipo_frac: vec![0.0; 2],
        entropy: vec![0.0; 2],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![3.0, 2.0],
        vatp: vec![0.0; 2],
        prot: vec![0.0; 2],
        mem: vec![0.0; 2],
        global_load: vec![0.0; 2],
    });
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0; 2],
        stress: vec![0.0; 2],
        cathepsin_membrane_imbalance: vec![0.0; 2],
        ldi: vec![1.0, 2.0],
    });
    ctx.coupling = Some(CouplingMetrics {
        ampk: vec![0.0; 2],
        ragulator: vec![0.0; 2],
        ampk_mtor_ratio: vec![0.0; 2],
        tfeb_lyso_alignment: vec![0.0; 2],
        rag_mtor_mismatch: vec![0.0; 2],
        locked_survival_index: vec![1.5, -1.0],
    });
    ctx
}

#[test]
fn test_mitophagy_reliance_and_erdi() {
    let mut ctx = base_ctx();
    run_stage10(&mut ctx).unwrap();
    let metrics = ctx.cross_organelle.as_ref().unwrap();

    let mr0 = 2.0 * 3.0;
    let erdi0 = 0.4 * mr0 + 0.3 * 1.5 + 0.3 * 1.0;
    assert!((metrics.mitophagy_reliance[0] - mr0).abs() < 1e-6);
    assert!((metrics.energy_recycling_dependency[0] - erdi0).abs() < 1e-6);
}

#[test]
fn test_degraded_mode() {
    let mut ctx = base_ctx();
    run_stage10(&mut ctx).unwrap();
    let metrics = ctx.cross_organelle.as_ref().unwrap();
    assert!(metrics.mito_lyso_imbalance[0].is_nan());
    assert!(metrics.mito_consumption_risk[0].is_nan());
}

#[test]
fn test_with_mito_metrics() {
    let mut ctx = base_ctx();
    ctx.mito = Some(MitoMetrics {
        mito_stress: vec![4.0, 2.0],
        mito_decay: vec![0.0; 2],
        mito_mass_proxy: vec![2.0, 1.0],
    });

    run_stage10(&mut ctx).unwrap();
    let metrics = ctx.cross_organelle.as_ref().unwrap();

    let mr0 = 2.0 * 3.0;
    let imbalance0 = 4.0 / (3.0 + 1e-6);
    let risk0 = mr0 / (2.0 + 1e-6);
    assert!((metrics.mito_lyso_imbalance[0] - imbalance0).abs() < 1e-6);
    assert!((metrics.mito_consumption_risk[0] - risk0).abs() < 1e-6);
}

#[test]
fn test_deterministic_order() {
    let mut ctx = base_ctx();
    run_stage10(&mut ctx).unwrap();
    let metrics = ctx.cross_organelle.as_ref().unwrap();
    assert_eq!(metrics.mitophagy_reliance.len(), 2);
    assert_eq!(metrics.mitophagy_reliance[0], 6.0);
}
