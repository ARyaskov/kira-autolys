use super::*;
use crate::config::thresholds::Thresholds;
use crate::model::cohort_view::{CohortView, SampleInfo};
use crate::model::ctx::{DrugVulnerabilityTag, TherapyResponseClass};

fn base_view() -> CohortView {
    CohortView {
        obs_ids: vec!["a".to_string(), "b".to_string()],
        sample_info: vec![
            SampleInfo {
                id: "a".to_string(),
                sample_group: "g1".to_string(),
                timepoint: 0,
            },
            SampleInfo {
                id: "b".to_string(),
                sample_group: "g2".to_string(),
                timepoint: 1,
            },
        ],
        thresholds: Thresholds::default(),
        ssm: vec![0.0, 0.0],
        afp: vec![0.0, 0.0],
        lds: vec![0.0, 0.0],
        ldi: vec![0.0, 0.0],
        lsi: vec![0.0, 0.0],
        erdi: vec![0.0, 0.0],
        classes: vec![
            LysosomeDependencyClass::AutophagyAdaptive,
            LysosomeDependencyClass::NonLysosomal,
        ],
        vulnerabilities: vec![vec![DrugVulnerabilityTag::LowVulnerability]; 2],
        therapy_response: Some(vec![TherapyResponseClass::NoResponse]),
    }
}

#[test]
fn test_by_group() {
    let view = base_view();
    let out_dir = std::env::temp_dir().join("kira_autolys_cohort_classes");
    let _ = std::fs::create_dir_all(&out_dir);
    write_classes(&view, &out_dir, true).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_classes.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!(json.get("by_sample_group").unwrap().get("g1").is_some());
    let _ = std::fs::remove_dir_all(out_dir);
}
