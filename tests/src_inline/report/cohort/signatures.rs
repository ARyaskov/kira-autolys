
use super::*;
use crate::config::thresholds::Thresholds;
use crate::model::cohort_view::{CohortView, SampleInfo};
use crate::model::ctx::{DrugVulnerabilityTag, TherapyResponseClass};

#[test]
fn test_signatures() {
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
        ssm: vec![1.2, 0.2, 1.2],
        afp: vec![0.0; 3],
        lds: vec![1.2, 0.2, 1.2],
        ldi: vec![1.2, 0.2, 1.2],
        lsi: vec![0.0; 3],
        erdi: vec![1.2, 0.2, 1.2],
        classes: vec![
            LysosomeDependencyClass::AutophagyAdaptive,
            LysosomeDependencyClass::Unclassified,
            LysosomeDependencyClass::AutophagyAdaptive,
        ],
        vulnerabilities: vec![vec![DrugVulnerabilityTag::LowVulnerability]; 3],
        therapy_response: Some(vec![TherapyResponseClass::NoResponse]),
    };

    let out_dir = std::env::temp_dir().join("kira_autolys_signatures");
    write_signatures(&view, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_signatures.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!((json.get("SSM_high").unwrap().as_f64().unwrap() - (2.0 / 3.0)).abs() < 1e-6);
    let _ = std::fs::remove_dir_all(out_dir);
}
