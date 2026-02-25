
use super::*;
use crate::model::ctx::{
    AutophagyMetrics, CouplingMetrics, CrossOrganelleMetrics, Ctx, LysosomalDamageMetrics,
    LysosomeMetrics, SurvivalMetrics,
};
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn test_percentiles() {
    let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    assert_eq!(percentile(&values, 0.5), 3.0);
    assert_eq!(percentile(&values, 0.9), 5.0);
}

#[test]
fn test_threshold_fraction() {
    let summary = summarize_metric(&[0.5, 1.0, 2.0], 1.0);
    assert!((summary.fraction_high - (2.0 / 3.0)).abs() < 1e-6);
}

#[test]
fn test_json_output() {
    let mut ctx = Ctx::default();
    ctx.survival = Some(SurvivalMetrics {
        prolif: vec![0.0; 3],
        apop_ready: vec![0.0; 3],
        z_afp: vec![0.0; 3],
        z_lds: vec![0.0; 3],
        z_prolif: vec![0.0; 3],
        z_apop: vec![0.0; 3],
        ssm: vec![0.5, 1.0, 1.5],
        flag_ssm_high: vec![false; 3],
        flag_afp_high: vec![false; 3],
        flag_lds_high: vec![false; 3],
        flag_prolif_low: vec![false; 3],
        flag_apop_low: vec![false; 3],
    });
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![0.0; 3],
        initiation: vec![0.0; 3],
        elongation: vec![0.0; 3],
        degradation: vec![0.0; 3],
        cargo: vec![0.0; 3],
        stall: vec![0.0; 3],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![0.0; 3],
        vatp: vec![0.0; 3],
        prot: vec![0.0; 3],
        mem: vec![0.0; 3],
        global_load: vec![0.0; 3],
    });
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0; 3],
        stress: vec![0.0; 3],
        cathepsin_membrane_imbalance: vec![0.0; 3],
        ldi: vec![0.0; 3],
    });
    ctx.coupling = Some(CouplingMetrics {
        ampk: vec![0.0; 3],
        ragulator: vec![0.0; 3],
        ampk_mtor_ratio: vec![0.0; 3],
        tfeb_lyso_alignment: vec![0.0; 3],
        rag_mtor_mismatch: vec![0.0; 3],
        locked_survival_index: vec![0.0; 3],
    });
    ctx.cross_organelle = Some(CrossOrganelleMetrics {
        mitophagy_reliance: vec![0.0; 3],
        mito_lyso_imbalance: vec![0.0; 3],
        mito_consumption_risk: vec![0.0; 3],
        energy_recycling_dependency: vec![0.0; 3],
    });

    let out_dir = temp_out_dir("metrics");
    write_cohort_metric_summary(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_metrics.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!(json.get("SSM").is_some());
    let _ = std::fs::remove_dir_all(out_dir);
}

fn temp_out_dir(label: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("kira_autolys_{label}_{nanos}"))
}
