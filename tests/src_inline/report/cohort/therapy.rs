
use super::*;
use crate::config::thresholds::Thresholds;
use crate::model::cohort_view::{CohortView, SampleInfo};
use crate::model::ctx::{DrugVulnerabilityTag, LysosomeDependencyClass};

#[test]
fn test_therapy_summary() {
    let view = CohortView {
        obs_ids: vec!["a".to_string()],
        sample_info: vec![SampleInfo {
            id: "a".to_string(),
            sample_group: "g".to_string(),
            timepoint: 0,
        }],
        thresholds: Thresholds::default(),
        ssm: vec![0.0],
        afp: vec![0.0],
        lds: vec![0.0],
        ldi: vec![0.0],
        lsi: vec![0.0],
        erdi: vec![0.0],
        classes: vec![LysosomeDependencyClass::Unclassified],
        vulnerabilities: vec![vec![DrugVulnerabilityTag::LowVulnerability]],
        therapy_response: Some(vec![TherapyResponseClass::AdaptiveSurvival]),
    };

    let out_dir = std::env::temp_dir().join("kira_autolys_cohort_therapy");
    write_therapy(&view, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_therapy_response.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!((json.get("ADAPTIVE_SURVIVAL").unwrap().as_f64().unwrap() - 1.0).abs() < 1e-6);
    let _ = std::fs::remove_dir_all(out_dir);
}
