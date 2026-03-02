use super::*;
use crate::model::ctx::{Ctx, LysosomeClassMetrics, SampleInfo};
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn test_global_fractions() {
    let mut ctx = Ctx::default();
    ctx.obs_ids = vec!["a".to_string(), "b".to_string(), "c".to_string()];
    ctx.lysosome_class = Some(LysosomeClassMetrics {
        class: vec![
            LysosomeDependencyClass::AutophagyAdaptive,
            LysosomeDependencyClass::AutophagyAdaptive,
            LysosomeDependencyClass::NonLysosomal,
        ],
    });

    let out_dir = temp_out_dir("cohort_classes");
    write_cohort_class_summary(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_classes.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    let global = json.get("global").unwrap();
    assert!(
        (global.get("AutophagyAdaptive").unwrap().as_f64().unwrap() - (2.0 / 3.0)).abs() < 1e-6
    );
    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_by_sample_group() {
    let mut ctx = Ctx::default();
    ctx.obs_ids = vec!["a".to_string(), "b".to_string()];
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
            sample_group: "g2".to_string(),
            treatment: None,
            timepoint: 0,
            timepoint_label: None,
        },
    ];
    ctx.lysosome_class = Some(LysosomeClassMetrics {
        class: vec![
            LysosomeDependencyClass::AutophagyAdaptive,
            LysosomeDependencyClass::NonLysosomal,
        ],
    });

    let out_dir = temp_out_dir("cohort_classes_grouped");
    write_cohort_class_summary(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("cohort_classes.json")).unwrap();
    let json: serde_json::Value = serde_json::from_str(&content).unwrap();
    let by_group = json.get("by_sample_group").unwrap();
    assert!(by_group.get("g1").is_some());
    assert!(by_group.get("g2").is_some());
    let _ = std::fs::remove_dir_all(out_dir);
}

fn temp_out_dir(label: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("kira_autolys_{label}_{nanos}"))
}
