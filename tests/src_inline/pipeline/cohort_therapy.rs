
use super::*;
use crate::model::ctx::{Ctx, TherapyDeltaMetrics};
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn test_therapy_summary() {
    let mut ctx = Ctx::default();
    ctx.therapy_delta = Some(TherapyDeltaMetrics {
        delta_ssm: vec![0.0; 3],
        delta_afp: vec![0.0; 3],
        delta_lds: vec![0.0; 3],
        delta_ldi: vec![0.0; 3],
        delta_erdi: vec![0.0; 3],
        delta_prolif: vec![0.0; 3],
        asi: vec![0.0; 3],
        response_class: vec![
            TherapyResponseClass::AdaptiveSurvival,
            TherapyResponseClass::AdaptiveSurvival,
            TherapyResponseClass::NoResponse,
        ],
        from_obs: vec![0; 3],
        to_obs: vec![1; 3],
        sample_group: vec!["g1".to_string(); 3],
        timepoint0: vec![0; 3],
        timepoint1: vec![1; 3],
    });

    let out_dir = temp_out_dir("therapy");
    write_cohort_therapy_summary(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_therapy_response.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!((json.get("ADAPTIVE_SURVIVAL").unwrap().as_f64().unwrap() - (2.0 / 3.0)).abs() < 1e-6);
    assert!((json.get("NO_RESPONSE").unwrap().as_f64().unwrap() - (1.0 / 3.0)).abs() < 1e-6);
    let _ = std::fs::remove_dir_all(out_dir);
}

fn temp_out_dir(label: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("kira_autolys_{label}_{nanos}"))
}
