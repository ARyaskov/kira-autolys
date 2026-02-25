
use super::*;
use crate::model::ctx::{
    AutophagyMetrics, AutophagySelectivityMetrics, CouplingMetrics, CrossOrganelleMetrics, Ctx,
    DrugVulnerabilityMetrics, DrugVulnerabilityTag, LysosomalDamageMetrics, LysosomeClassMetrics,
    LysosomeDependencyClass, LysosomeMetrics, RegulatoryMetrics, SampleInfo, SurvivalMetrics,
    SurvivalStats,
};

fn base_ctx() -> Ctx {
    let mut ctx = Ctx::default();
    ctx.input_type = Some(crate::model::ctx::InputType::Mtx10x);
    ctx.mode = Some(crate::model::ctx::Mode::Cell);
    ctx.n_obs = 2;
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
            sample_group: "g1".to_string(),
            treatment: None,
            timepoint: 1,
            timepoint_label: None,
        },
    ];
    ctx.autophagy = Some(AutophagyMetrics {
        afp: vec![1.0, 2.0],
        initiation: vec![0.1, 0.2],
        elongation: vec![0.3, 0.4],
        degradation: vec![0.5, 0.6],
        cargo: vec![0.7, 0.8],
        stall: vec![0.0, 0.1],
    });
    ctx.lysosome = Some(LysosomeMetrics {
        lds: vec![2.0, 3.0],
        vatp: vec![0.1, 0.2],
        prot: vec![0.3, 0.4],
        mem: vec![0.5, 0.6],
        global_load: vec![1.0, 1.1],
    });
    ctx.regulatory = Some(RegulatoryMetrics {
        tfeb: vec![0.1, 0.2],
        mtor: vec![0.3, 0.4],
        tfeb_act: vec![0.0, 0.0],
        mtor_supp: vec![0.0, 0.0],
        tfeb_mtor_diff: vec![0.1, -0.2],
        tfeb_mtor_ratio: vec![1.0, 2.0],
        median_tfeb: 0.15,
        median_mtor: 0.35,
    });
    ctx.survival = Some(SurvivalMetrics {
        prolif: vec![0.0, 0.0],
        apop_ready: vec![0.0, 0.0],
        z_afp: vec![0.1, 0.2],
        z_lds: vec![0.3, 0.4],
        z_prolif: vec![0.5, 0.6],
        z_apop: vec![0.7, 0.8],
        ssm: vec![0.9, 1.0],
        flag_ssm_high: vec![false, true],
        flag_afp_high: vec![false, false],
        flag_lds_high: vec![false, true],
        flag_prolif_low: vec![true, false],
        flag_apop_low: vec![false, false],
    });
    ctx.survival_stats = Some(SurvivalStats {
        mean_afp: 1.5,
        std_afp: 0.1,
        mean_lds: 2.5,
        std_lds: 0.2,
        mean_prolif: 0.0,
        std_prolif: 1.0,
        mean_apop: 0.0,
        std_apop: 1.0,
    });
    ctx.lysosomal_damage = Some(LysosomalDamageMetrics {
        lmp: vec![0.0, 0.0],
        stress: vec![0.0, 0.0],
        cathepsin_membrane_imbalance: vec![0.0, 0.0],
        ldi: vec![0.0, 0.0],
    });
    ctx.autophagy_selectivity = Some(AutophagySelectivityMetrics {
        mitophagy: vec![0.0, 0.0],
        aggrephagy: vec![0.0, 0.0],
        erphagy: vec![0.0, 0.0],
        ferritinophagy: vec![0.0, 0.0],
        lipophagy: vec![0.0, 0.0],
        mito_frac: vec![0.2, 0.3],
        aggre_frac: vec![0.2, 0.2],
        er_frac: vec![0.2, 0.2],
        ferr_frac: vec![0.2, 0.2],
        lipo_frac: vec![0.2, 0.1],
        entropy: vec![1.0, 1.1],
    });
    ctx.coupling = Some(CouplingMetrics {
        ampk: vec![0.0, 0.0],
        ragulator: vec![0.0, 0.0],
        ampk_mtor_ratio: vec![1.0, 1.1],
        tfeb_lyso_alignment: vec![0.5, 0.6],
        rag_mtor_mismatch: vec![0.0, 0.0],
        locked_survival_index: vec![0.7, 0.8],
    });
    ctx.cross_organelle = Some(CrossOrganelleMetrics {
        mitophagy_reliance: vec![0.1, 0.2],
        mito_lyso_imbalance: vec![f32::NAN, f32::NAN],
        mito_consumption_risk: vec![f32::NAN, f32::NAN],
        energy_recycling_dependency: vec![0.4, 0.5],
    });
    ctx.lysosome_class = Some(LysosomeClassMetrics {
        class: vec![LysosomeDependencyClass::Unclassified; 2],
    });
    ctx.drug_vulnerability = Some(DrugVulnerabilityMetrics {
        tags: vec![vec![DrugVulnerabilityTag::LowVulnerability]; 2],
    });
    ctx
}

fn temp_out_dir() -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    path.push(format!("kira_autolys_v2_{nanos}"));
    path
}

#[test]
fn test_core_header() {
    let ctx = base_ctx();
    let out_dir = temp_out_dir();
    run_stage14_export_v2(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("core_scores.tsv")).unwrap();
    let header = content.lines().next().unwrap();
    assert_eq!(header, CORE_HEADER);
    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_summary_v2_structure() {
    let ctx = base_ctx();
    let out_dir = temp_out_dir();
    run_stage14_export_v2(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("summary_v2.json")).unwrap();
    let value: serde_json::Value = serde_json::from_str(&content).unwrap();
    let schema = value
        .get("tool")
        .and_then(|v| v.get("schema"))
        .and_then(|v| v.as_str())
        .unwrap();
    assert_eq!(schema, "v2");
    assert!(value.get("observations").is_some());
    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_classification_header() {
    let ctx = base_ctx();
    let out_dir = temp_out_dir();
    run_stage14_export_v2(&ctx, &out_dir).unwrap();
    let content = std::fs::read_to_string(out_dir.join("classification.tsv")).unwrap();
    let header = content.lines().next().unwrap();
    assert_eq!(header, CLASS_HEADER);
    let _ = std::fs::remove_dir_all(out_dir);
}
