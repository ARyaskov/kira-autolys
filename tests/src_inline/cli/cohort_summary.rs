use super::*;
use std::fs;

fn write_file(path: &Path, content: &str) {
    fs::write(path, content).unwrap();
}

#[test]
fn test_cohort_summary_end_to_end() {
    let tmp = std::env::temp_dir().join("kira_autolys_cohort_cli");
    let _ = fs::remove_dir_all(&tmp);
    fs::create_dir_all(&tmp).unwrap();

    let summary = r#"{
  "tool": {"name":"kira-autolys","version":"0.1.0","schema":"v2"},
  "input": {"type":"mtx","mode":"cell","n_obs":2,
    "samples":[
      {"id":"a","sample_group":"g1","timepoint":0},
      {"id":"b","sample_group":"g1","timepoint":1}
    ]},
  "thresholds": {
    "ssm_high":1.0,"afp_high":1.0,"lds_high":1.0,"prolif_low":-0.5,"apop_low":-0.5,
    "erdi_high":1.0,"lds_low":0.3,"afp_low":0.3,"ldi_high":1.0,"stall_high":0.5,
    "vatp_high":1.0,"cat_imbalance_high":1.2,"mito_reliance_high":1.0,"lsi_high":1.0
  },
  "dataset_stats": {"core": {"mean_afp":0.0,"std_afp":1.0,"mean_lds":0.0,"std_lds":1.0,"mean_prolif":0.0,"std_prolif":1.0,"mean_apop":0.0,"std_apop":1.0},
                    "regulatory": {"median_tfeb":0.0,"median_mtor":0.0}},
  "observations": [{"id":"a"},{"id":"b"}]
}"#;
    write_file(&tmp.join("summary_v2.json"), summary);

    write_file(
        &tmp.join("core_scores.tsv"),
        "id\tAFP\tLDS\tTFEB\tmTOR\tSSM\tz_AFP\tz_LDS\tz_PROLIF\tz_APOP\n\
a\t1.0\t1.0\t0.0\t0.0\t1.2\t0.0\t0.0\t0.0\t0.0\n\
b\t0.5\t0.5\t0.0\t0.0\t0.2\t0.0\t0.0\t0.0\t0.0\n",
    );
    write_file(
        &tmp.join("damage.tsv"),
        "id\tLDI\tLMP\tstress\tcathepsin_membrane_imbalance\n\
a\t1.2\t0.0\t0.0\t0.0\n\
b\t0.2\t0.0\t0.0\t0.0\n",
    );
    write_file(
        &tmp.join("coupling.tsv"),
        "id\tLSI\tampk_mtor_ratio\ttfeb_lyso_alignment\n\
a\t0.5\t0.0\t0.0\n\
b\t0.1\t0.0\t0.0\n",
    );
    write_file(
        &tmp.join("cross_organelle.tsv"),
        "id\tERDI\tmitophagy_reliance\tmito_lyso_imbalance\tmito_consumption_risk\n\
a\t1.2\t0.0\tNaN\tNaN\n\
b\t0.2\t0.0\tNaN\tNaN\n",
    );
    write_file(
        &tmp.join("classification.tsv"),
        "id\tlysosome_dependency_class\n\
a\tAutophagyAdaptive\n\
b\tUnclassified\n",
    );
    write_file(
        &tmp.join("vulnerabilities.tsv"),
        "id\tvulnerability_tags\n\
a\tAutophagyInhibitionSensitive\n\
b\tLowVulnerability\n",
    );

    run_cohort_summary(&tmp, true).unwrap();
    assert!(tmp.join("cohort/cohort_classes.json").is_file());
    assert!(tmp.join("cohort/cohort_metrics.json").is_file());
    assert!(tmp.join("cohort/cohort_signatures.json").is_file());
    assert!(tmp.join("cohort/cohort_vulnerabilities.json").is_file());

    let _ = fs::remove_dir_all(&tmp);
}
