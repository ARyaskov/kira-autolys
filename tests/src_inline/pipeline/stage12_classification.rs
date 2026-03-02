use super::*;
use crate::model::ctx::{
    AutophagyMetrics, CouplingMetrics, CrossOrganelleMetrics, Ctx, LysosomalDamageMetrics,
    LysosomeMetrics, SurvivalMetrics,
};

fn base_ctx() -> Ctx {
    let mut ctx = Ctx::default();
    ctx.n_obs = 1;
    ctx.obs_ids = vec!["a".to_string()];
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![0.0],
        initiation: vec![0.0],
        elongation: vec![0.0],
        degradation: vec![0.0],
        cargo: vec![0.0],
        stall: vec![0.0],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0],
        vatp: vec![0.0],
        prot: vec![0.0],
        mem: vec![0.0],
        global_load: vec![0.0],
    });
    ctx.survival = Some(SurvivalMetrics {
        prolif: vec![0.0],
        apop_ready: vec![0.0],
        z_afp: vec![0.0],
        z_lds: vec![0.0],
        z_prolif: vec![0.0],
        z_apop: vec![0.0],
        ssm: vec![0.0],
        flag_ssm_high: vec![false],
        flag_afp_high: vec![false],
        flag_lds_high: vec![false],
        flag_prolif_low: vec![false],
        flag_apop_low: vec![false],
    });
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0],
        stress: vec![0.0],
        cathepsin_membrane_imbalance: vec![0.0],
        ldi: vec![0.0],
    });
    ctx.coupling = Some(CouplingMetrics {
        ampk: vec![0.0],
        ragulator: vec![0.0],
        ampk_mtor_ratio: vec![0.0],
        tfeb_lyso_alignment: vec![0.0],
        rag_mtor_mismatch: vec![0.0],
        locked_survival_index: vec![0.0],
    });
    ctx.cross_organelle = Some(CrossOrganelleMetrics {
        mitophagy_reliance: vec![0.0],
        mito_lyso_imbalance: vec![0.0],
        mito_consumption_risk: vec![0.0],
        energy_recycling_dependency: vec![0.0],
    });
    ctx
}

#[test]
fn test_rule_priority_energy_recycling() {
    let mut ctx = base_ctx();
    ctx.cross_organelle
        .as_mut()
        .unwrap()
        .energy_recycling_dependency[0] = 2.0;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::EnergyRecyclingDependent
    );
}

#[test]
fn test_lysosome_overloaded() {
    let mut ctx = base_ctx();
    ctx.lysosome.as_mut().unwrap().lds[0] = 1.2;
    ctx.lysosomal_damage.as_mut().unwrap().ldi[0] = 1.2;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::LysosomeOverloaded
    );
}

#[test]
fn test_autophagy_adaptive() {
    let mut ctx = base_ctx();
    ctx.survival.as_mut().unwrap().flag_ssm_high[0] = true;
    ctx.autophagy.as_mut().unwrap().afp[0] = 1.2;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::AutophagyAdaptive
    );
}

#[test]
fn test_stalled_autophagy() {
    let mut ctx = base_ctx();
    ctx.autophagy.as_mut().unwrap().stall[0] = 0.6;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::StalledAutophagy
    );
}

#[test]
fn test_lysosome_dependent() {
    let mut ctx = base_ctx();
    ctx.lysosome.as_mut().unwrap().lds[0] = 1.2;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::LysosomeDependent
    );
}

#[test]
fn test_non_lysosomal() {
    let mut ctx = base_ctx();
    ctx.lysosome.as_mut().unwrap().lds[0] = 0.1;
    ctx.autophagy.as_mut().unwrap().afp[0] = 0.1;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::NonLysosomal
    );
}

#[test]
fn test_unclassified() {
    let mut ctx = base_ctx();
    ctx.lysosome.as_mut().unwrap().lds[0] = 0.6;
    ctx.autophagy.as_mut().unwrap().afp[0] = 0.6;
    run_stage12(&mut ctx).unwrap();
    assert_eq!(
        ctx.lysosome_class.as_ref().unwrap().class[0],
        LysosomeDependencyClass::Unclassified
    );
}
