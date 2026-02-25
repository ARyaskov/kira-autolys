
use super::*;
use crate::model::cli::RunMode;
use crate::model::ctx::NormalizedExpr;
use crate::model::ctx::{
    AutophagyMetrics, Ctx, LysosomeMetrics, RegulatoryMetrics, SurvivalMetrics, SurvivalStats,
};

struct DummyNorm {
    rows: Vec<Vec<(usize, f32)>>,
}

impl NormalizedExpr for DummyNorm {
    fn n_obs(&self) -> usize {
        self.rows.len()
    }
    fn for_each_in_obs(&self, obs_idx: usize, f: &mut dyn FnMut(usize, f32)) {
        for &(gene, value) in &self.rows[obs_idx] {
            f(gene, value);
        }
    }
}

fn test_ctx() -> Ctx {
    let mut ctx = Ctx::default();
    ctx.input_type = Some(crate::model::ctx::InputType::Mtx10x);
    ctx.mode = Some(crate::model::ctx::Mode::Cell);
    ctx.n_obs = 2;
    ctx.n_vars = 3;
    ctx.obs_ids = vec!["cell1".to_string(), "cell2".to_string()];
    ctx.obs_libsize = vec![100.0, 200.0];
    ctx.obs_nnz = vec![10, 20];
    ctx.obs_expressed_genes = vec![10, 20];
    ctx.gene_symbols = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    ctx.gene_sets = vec![crate::model::ctx::GeneSetResolved {
        name: "PANEL1".to_string(),
        indices: vec![0, 1],
    }];
    ctx.gene_panel_mask = vec![1, 1, 0];
    ctx.missing_genes = vec![crate::model::ctx::MissingGenes {
        panel: "PANEL1".to_string(),
        missing: vec!["Z".to_string()],
    }];
    ctx.normalized = Some(Box::new(DummyNorm {
        rows: vec![vec![(0, 1.0), (1, 2.0)], vec![(0, 0.5), (1, 1.0)]],
    }));
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
        z_afp: vec![0.0, 0.0],
        z_lds: vec![0.0, 0.0],
        z_prolif: vec![0.0, 0.0],
        z_apop: vec![0.0, 0.0],
        ssm: vec![0.0, 1.0],
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
    ctx
}

fn temp_out_dir() -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    path.push(format!("kira_autolys_test_{nanos}"));
    path
}

#[test]
fn test_scores_tsv_header() {
    let ctx = test_ctx();
    let out_dir = temp_out_dir();
    run_stage6(&ctx, &out_dir, true).unwrap();

    let content = std::fs::read_to_string(out_dir.join("scores.tsv")).unwrap();
    let header = content.lines().next().unwrap();
    assert_eq!(header, TSV_HEADER);

    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_summary_json_structure() {
    let ctx = test_ctx();
    let out_dir = temp_out_dir();
    run_stage6(&ctx, &out_dir, true).unwrap();

    let content = std::fs::read_to_string(out_dir.join("summary.json")).unwrap();
    let value: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!(value.get("tool").is_some());
    assert!(value.get("dataset_stats").is_some());
    assert!(value.get("top_findings").is_some());
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out_dir.join("metadata.json")).unwrap())
            .unwrap();
    assert!(meta.get("input").is_some());

    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_deterministic_ordering() {
    let ctx = test_ctx();
    let out_dir = temp_out_dir();
    run_stage6(&ctx, &out_dir, true).unwrap();

    let content = std::fs::read_to_string(out_dir.join("scores.tsv")).unwrap();
    let mut lines = content.lines();
    lines.next();
    let first = lines.next().unwrap();
    assert!(first.starts_with("cell1\t"));

    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_pipeline_autolys_header_order() {
    let mut ctx = test_ctx();
    ctx.run_mode = RunMode::Pipeline;
    let out_dir = temp_out_dir();
    run_stage6(&ctx, &out_dir, false).unwrap();
    let content = std::fs::read_to_string(out_dir.join("autolys.tsv")).unwrap();
    let header = content.lines().next().unwrap();
    assert_eq!(header, PIPELINE_AUTOLYS_HEADER);
    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_pipeline_summary_schema() {
    let mut ctx = test_ctx();
    ctx.run_mode = RunMode::Pipeline;
    let out_dir = temp_out_dir();
    run_stage6(&ctx, &out_dir, false).unwrap();
    let value: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out_dir.join("summary.json")).unwrap())
            .unwrap();
    assert_eq!(value["tool"], "kira-autolys");
    assert!(value.get("input").is_some());
    assert_eq!(value["input"]["mode"], "pipeline");
    assert!(value.get("distributions").is_some());
    assert!(value.get("regimes").is_some());
    assert!(value.get("qc").is_some());
    assert!(value["qc"].get("coverage_median").is_some());
    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_pipeline_step_schema() {
    let mut ctx = test_ctx();
    ctx.run_mode = RunMode::Pipeline;
    let out_dir = temp_out_dir();
    run_stage6(&ctx, &out_dir, false).unwrap();
    let value: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out_dir.join("pipeline_step.json")).unwrap())
            .unwrap();
    assert_eq!(value["tool"], "kira-autolys");
    assert_eq!(value["mode"], "pipeline");
    assert_eq!(value["artifacts"]["summary"], "summary.json");
    assert_eq!(value["artifacts"]["primary_metrics"], "autolys.tsv");
    assert_eq!(value["cell_metrics"]["id_column"], "id");
    assert_eq!(value["cell_metrics"]["regime_column"], "regime");
    assert_eq!(value["cell_metrics"]["confidence_column"], "confidence");
    assert_eq!(value["cell_metrics"]["flag_column"], "flags");
    assert_eq!(value["regimes"]["source"], "autolys.tsv");
    assert_eq!(value["regimes"]["column"], "regime");
    let _ = std::fs::remove_dir_all(out_dir);
}

#[test]
fn test_pipeline_deterministic_outputs() {
    let mut ctx = test_ctx();
    ctx.run_mode = RunMode::Pipeline;
    let out_a = temp_out_dir();
    let out_b = temp_out_dir();
    run_stage6(&ctx, &out_a, false).unwrap();
    run_stage6(&ctx, &out_b, false).unwrap();

    for name in [
        "autolys.tsv",
        "summary.json",
        "panels_report.tsv",
        "pipeline_step.json",
    ] {
        let a = std::fs::read(out_a.join(name)).unwrap();
        let b = std::fs::read(out_b.join(name)).unwrap();
        assert_eq!(a, b, "mismatch in {name}");
    }
    let _ = std::fs::remove_dir_all(out_a);
    let _ = std::fs::remove_dir_all(out_b);
}
