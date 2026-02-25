
use super::*;
use crate::config::thresholds::Thresholds;
use crate::model::cohort_view::{CohortView, SampleInfo};
use crate::model::ctx::{DrugVulnerabilityTag, LysosomeDependencyClass, TherapyResponseClass};

#[test]
fn test_metric_summary() {
    let view = CohortView {
        obs_ids: vec!["a".to_string(), "b".to_string(), "c".to_string()],
        sample_info: vec![
            SampleInfo {
                id: "a".to_string(),
                sample_group: "g".to_string(),
                timepoint: 0,
            },
            SampleInfo {
                id: "b".to_string(),
                sample_group: "g".to_string(),
                timepoint: 1,
            },
            SampleInfo {
                id: "c".to_string(),
                sample_group: "g".to_string(),
                timepoint: 2,
            },
        ],
        thresholds: Thresholds::default(),
        ssm: vec![0.5, 1.0, 1.5],
        afp: vec![0.0; 3],
        lds: vec![0.0; 3],
        ldi: vec![0.0; 3],
        lsi: vec![0.0; 3],
        erdi: vec![0.0; 3],
        classes: vec![LysosomeDependencyClass::Unclassified; 3],
        vulnerabilities: vec![vec![DrugVulnerabilityTag::LowVulnerability]; 3],
        therapy_response: Some(vec![TherapyResponseClass::NoResponse]),
    };

    let out_dir = std::env::temp_dir().join("kira_autolys_metric_summary");
    let _ = std::fs::create_dir_all(&out_dir);
    write_metrics(&view, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_metrics.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!(json.get("SSM").is_some());
    let _ = std::fs::remove_dir_all(out_dir);
}
