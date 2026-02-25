
use super::*;
use crate::model::ctx::{
    AutophagyMetrics, CouplingMetrics, CrossOrganelleMetrics, Ctx, LysosomalDamageMetrics,
    LysosomeMetrics, SampleInfo, SurvivalMetrics,
};

fn base_ctx() -> Ctx {
    let mut ctx = Ctx::default();
    ctx.n_obs = 3;
    ctx.obs_ids = vec!["a".to_string(), "b".to_string(), "c".to_string()];
    ctx.samples = vec![
        SampleInfo {
            id: "a".to_string(),
            sample_group: "g1".to_string(),
            treatment: None,
            timepoint: 0,
            timepoint_label: None,
        },
        SampleInfo {
            id: "b".to_string(),
            sample_group: "g1".to_string(),
            treatment: None,
            timepoint: 1,
            timepoint_label: None,
        },
        SampleInfo {
            id: "c".to_string(),
            sample_group: "g2".to_string(),
            treatment: None,
            timepoint: 0,
            timepoint_label: None,
        },
    ];
    ctx.survival = Some(SurvivalMetrics {
        prolif: vec![2.0, 1.0, 1.0],
        apop_ready: vec![0.0; 3],
        z_afp: vec![0.0; 3],
        z_lds: vec![0.0; 3],
        z_prolif: vec![0.0; 3],
        z_apop: vec![0.0; 3],
        ssm: vec![0.5, 0.2, 0.1],
        flag_ssm_high: vec![false; 3],
        flag_afp_high: vec![false; 3],
        flag_lds_high: vec![false; 3],
        flag_prolif_low: vec![false; 3],
        flag_apop_low: vec![false; 3],
    });
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![1.0, 2.0, 1.5],
        initiation: vec![0.0; 3],
        elongation: vec![0.0; 3],
        degradation: vec![0.0; 3],
        cargo: vec![0.0; 3],
        stall: vec![0.0; 3],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![1.0, 1.5, 1.0],
        vatp: vec![0.0; 3],
        prot: vec![0.0; 3],
        mem: vec![0.0; 3],
        global_load: vec![0.0; 3],
    });
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0; 3],
        stress: vec![0.0; 3],
        cathepsin_membrane_imbalance: vec![0.0; 3],
        ldi: vec![0.2, 0.6, 0.1],
    });
    ctx.cross_organelle = Some(CrossOrganelleMetrics {
        mitophagy_reliance: vec![0.0; 3],
        mito_lyso_imbalance: vec![0.0; 3],
        mito_consumption_risk: vec![0.0; 3],
        energy_recycling_dependency: vec![0.4, 0.8, 0.2],
    });
    ctx.coupling = Some(CouplingMetrics {
        ampk: vec![0.0; 3],
        ragulator: vec![0.0; 3],
        ampk_mtor_ratio: vec![0.0; 3],
        tfeb_lyso_alignment: vec![0.0; 3],
        rag_mtor_mismatch: vec![0.0; 3],
        locked_survival_index: vec![0.0; 3],
    });
    ctx
}

#[test]
fn test_grouping_and_deltas() {
    let mut ctx = base_ctx();
    run_stage11(&mut ctx).unwrap();
    let metrics = ctx.therapy_delta.as_ref().unwrap();

    assert_eq!(metrics.delta_ssm.len(), 1);
    assert_eq!(metrics.from_obs[0], 0);
    assert_eq!(metrics.to_obs[0], 1);
    assert_eq!(metrics.delta_ssm[0], 0.2 - 0.5);
    assert_eq!(metrics.delta_afp[0], 2.0 - 1.0);
    assert_eq!(metrics.delta_lds[0], 1.5 - 1.0);
    assert_eq!(metrics.delta_ldi[0], 0.6 - 0.2);
    assert_eq!(metrics.delta_erdi[0], 0.8 - 0.4);
    assert_eq!(metrics.delta_prolif[0], 1.0 - 2.0);
}

#[test]
fn test_asi_and_response() {
    let mut ctx = base_ctx();
    run_stage11(&mut ctx).unwrap();
    let metrics = ctx.therapy_delta.as_ref().unwrap();
    let asi = 0.4 * (0.2 - 0.5) + 0.3 * (2.0 - 1.0) + 0.3 * (1.5 - 1.0);
    assert!((metrics.asi[0] - asi).abs() < 1e-6);
    assert_eq!(
        metrics.response_class[0],
        TherapyResponseClass::CytotoxicResponse
    );
}

#[test]
fn test_deterministic_order() {
    let mut ctx = base_ctx();
    run_stage11(&mut ctx).unwrap();
    let metrics = ctx.therapy_delta.as_ref().unwrap();
    assert_eq!(metrics.sample_group[0], "g1");
    assert_eq!(metrics.timepoint0[0], 0);
    assert_eq!(metrics.timepoint1[0], 1);
}
