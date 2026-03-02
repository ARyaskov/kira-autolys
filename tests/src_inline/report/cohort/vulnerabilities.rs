use super::*;
use crate::config::thresholds::Thresholds;
use crate::model::cohort_view::{CohortView, SampleInfo};
use crate::model::ctx::{LysosomeDependencyClass, TherapyResponseClass};

#[test]
fn test_vulnerability_fractions() {
    let view = CohortView {
        obs_ids: vec!["a".to_string(), "b".to_string()],
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
        ],
        thresholds: Thresholds::default(),
        ssm: vec![0.0, 0.0],
        afp: vec![0.0, 0.0],
        lds: vec![0.0, 0.0],
        ldi: vec![0.0, 0.0],
        lsi: vec![0.0, 0.0],
        erdi: vec![0.0, 0.0],
        classes: vec![LysosomeDependencyClass::Unclassified; 2],
        vulnerabilities: vec![
            vec![
                DrugVulnerabilityTag::AutophagyInhibitionSensitive,
                DrugVulnerabilityTag::LowVulnerability,
            ],
            vec![
                DrugVulnerabilityTag::AutophagyInhibitionSensitive,
                DrugVulnerabilityTag::LysosomalDestabilizationSensitive,
            ],
        ],
        therapy_response: Some(vec![TherapyResponseClass::NoResponse]),
    };

    let out_dir = std::env::temp_dir().join("kira_autolys_vuln_summary");
    write_vulnerabilities(&view, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_vulnerabilities.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!(
        (json
            .get("AutophagyInhibitionSensitive")
            .unwrap()
            .as_f64()
            .unwrap()
            - 1.0)
            .abs()
            < 1e-6
    );
    assert!((json.get("LowVulnerability").unwrap().as_f64().unwrap() - 0.5).abs() < 1e-6);
    let _ = std::fs::remove_dir_all(out_dir);
}
