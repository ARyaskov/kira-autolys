use crate::config::thresholds::Thresholds;
use std::path::Path;

use serde::Deserialize;

use crate::model::ctx::{DrugVulnerabilityTag, LysosomeDependencyClass, TherapyResponseClass};
use crate::model::schema_v2::ThresholdsInfo;
use crate::stage_error::StageError;

#[derive(Debug)]
pub struct CohortView {
    pub obs_ids: Vec<String>,
    pub sample_info: Vec<SampleInfo>,
    pub thresholds: Thresholds,
    pub ssm: Vec<f32>,
    pub afp: Vec<f32>,
    pub lds: Vec<f32>,
    pub ldi: Vec<f32>,
    pub lsi: Vec<f32>,
    pub erdi: Vec<f32>,
    pub classes: Vec<LysosomeDependencyClass>,
    pub vulnerabilities: Vec<Vec<DrugVulnerabilityTag>>,
    pub therapy_response: Option<Vec<TherapyResponseClass>>,
}

#[derive(Debug, Clone)]
pub struct SampleInfo {
    pub id: String,
    pub sample_group: String,
    pub timepoint: i32,
}

#[derive(Debug, Deserialize)]
struct SummaryV2 {
    tool: SummaryTool,
    input: SummaryInput,
    thresholds: ThresholdsInfo,
    observations: Vec<ObsStub>,
}

#[derive(Debug, Deserialize)]
struct SummaryTool {
    schema: String,
}

#[derive(Debug, Deserialize)]
struct SummaryInput {
    samples: Vec<SummarySample>,
}

#[derive(Debug, Deserialize)]
struct SummarySample {
    id: String,
    sample_group: String,
    timepoint: i32,
}

#[derive(Debug, Deserialize)]
struct ObsStub {
    id: String,
}

pub fn load_cohort_view(input_dir: &Path) -> Result<CohortView, StageError> {
    let summary_path = input_dir.join("summary_v2.json");
    if !summary_path.is_file() {
        return Err(StageError::Validation(
            "summary_v2.json missing; run schema v2 export first".to_string(),
        ));
    }

    let summary_text = std::fs::read_to_string(&summary_path)?;
    let summary: SummaryV2 = serde_json::from_str(&summary_text)
        .map_err(|err| StageError::Format(format!("invalid summary_v2.json: {err}")))?;
    if summary.tool.schema != "v2" {
        return Err(StageError::Validation(
            "summary_v2.json schema is not v2".to_string(),
        ));
    }

    let obs_ids: Vec<String> = summary.observations.into_iter().map(|o| o.id).collect();
    let sample_info: Vec<SampleInfo> = summary
        .input
        .samples
        .into_iter()
        .map(|s| SampleInfo {
            id: s.id,
            sample_group: s.sample_group,
            timepoint: s.timepoint,
        })
        .collect();

    if !sample_info.is_empty() && sample_info.len() != obs_ids.len() {
        return Err(StageError::Validation(
            "sample metadata length does not match observations".to_string(),
        ));
    }

    let thresholds = thresholds_from_info(&summary.thresholds);

    let ssm = load_metric_col(&input_dir.join("core_scores.tsv"), 5, &obs_ids, "SSM")?;
    let afp = load_metric_col(&input_dir.join("core_scores.tsv"), 1, &obs_ids, "AFP")?;
    let lds = load_metric_col(&input_dir.join("core_scores.tsv"), 2, &obs_ids, "LDS")?;
    let ldi = load_metric_col(&input_dir.join("damage.tsv"), 1, &obs_ids, "LDI")?;
    let lsi = load_metric_col(&input_dir.join("coupling.tsv"), 1, &obs_ids, "LSI")?;
    let erdi = load_metric_col(&input_dir.join("cross_organelle.tsv"), 1, &obs_ids, "ERDI")?;

    let classes = load_classes(&input_dir.join("classification.tsv"), &obs_ids)?;
    let vulnerabilities = load_vulnerabilities(&input_dir.join("vulnerabilities.tsv"), &obs_ids)?;

    let therapy_response = {
        let path = input_dir.join("therapy_delta.tsv");
        if path.is_file() {
            Some(load_therapy_responses(&path)?)
        } else {
            None
        }
    };

    Ok(CohortView {
        obs_ids,
        sample_info,
        thresholds,
        ssm,
        afp,
        lds,
        ldi,
        lsi,
        erdi,
        classes,
        vulnerabilities,
        therapy_response,
    })
}

fn thresholds_from_info(info: &ThresholdsInfo) -> Thresholds {
    Thresholds {
        ssm_high: info.ssm_high,
        afp_high: info.afp_high,
        lds_high: info.lds_high,
        prolif_low: info.prolif_low,
        apop_low: info.apop_low,
        erdi_high: info.erdi_high,
        lds_low: info.lds_low,
        afp_low: info.afp_low,
        ldi_high: info.ldi_high,
        stall_high: info.stall_high,
        vatp_high: info.vatp_high,
        cat_imbalance_high: info.cat_imbalance_high,
        mito_reliance_high: info.mito_reliance_high,
        lsi_high: info.lsi_high,
    }
}

fn load_metric_col(
    path: &Path,
    col_idx: usize,
    obs_ids: &[String],
    label: &str,
) -> Result<Vec<f32>, StageError> {
    let lines = read_tsv(path)?;
    if lines.len().saturating_sub(1) != obs_ids.len() {
        return Err(StageError::Validation(format!(
            "row count mismatch in {}",
            path.display()
        )));
    }
    let mut values = Vec::with_capacity(obs_ids.len());
    for (row_idx, row) in lines.iter().enumerate().skip(1) {
        let id = row.get(0).ok_or_else(|| {
            StageError::Format(format!("missing id column in {}", path.display()))
        })?;
        let expected_id = &obs_ids[row_idx - 1];
        if id != expected_id {
            return Err(StageError::Validation(format!(
                "obs_id mismatch in {}: expected {}, got {}",
                path.display(),
                expected_id,
                id
            )));
        }
        let value = row.get(col_idx).ok_or_else(|| {
            StageError::Format(format!("missing {label} column in {}", path.display()))
        })?;
        values.push(parse_f32(value, path)?);
    }
    Ok(values)
}

fn load_classes(
    path: &Path,
    obs_ids: &[String],
) -> Result<Vec<LysosomeDependencyClass>, StageError> {
    let lines = read_tsv(path)?;
    if lines.len().saturating_sub(1) != obs_ids.len() {
        return Err(StageError::Validation(format!(
            "row count mismatch in {}",
            path.display()
        )));
    }
    let mut classes = Vec::with_capacity(obs_ids.len());
    for (row_idx, row) in lines.iter().enumerate().skip(1) {
        let id = row
            .get(0)
            .ok_or_else(|| StageError::Format("missing id".to_string()))?;
        let expected_id = &obs_ids[row_idx - 1];
        if id != expected_id {
            return Err(StageError::Validation(format!(
                "obs_id mismatch in {}: expected {}, got {}",
                path.display(),
                expected_id,
                id
            )));
        }
        let value = row
            .get(1)
            .ok_or_else(|| StageError::Format("missing class".to_string()))?;
        classes.push(parse_class(value, path)?);
    }
    Ok(classes)
}

fn load_vulnerabilities(
    path: &Path,
    obs_ids: &[String],
) -> Result<Vec<Vec<DrugVulnerabilityTag>>, StageError> {
    let lines = read_tsv(path)?;
    if lines.len().saturating_sub(1) != obs_ids.len() {
        return Err(StageError::Validation(format!(
            "row count mismatch in {}",
            path.display()
        )));
    }
    let mut tags = Vec::with_capacity(obs_ids.len());
    for (row_idx, row) in lines.iter().enumerate().skip(1) {
        let id = row
            .get(0)
            .ok_or_else(|| StageError::Format("missing id".to_string()))?;
        let expected_id = &obs_ids[row_idx - 1];
        if id != expected_id {
            return Err(StageError::Validation(format!(
                "obs_id mismatch in {}: expected {}, got {}",
                path.display(),
                expected_id,
                id
            )));
        }
        let value = row
            .get(1)
            .ok_or_else(|| StageError::Format("missing tags".to_string()))?;
        let mut obs_tags = Vec::new();
        if !value.is_empty() {
            for token in value.split(',') {
                obs_tags.push(parse_vulnerability(token, path)?);
            }
        }
        tags.push(obs_tags);
    }
    Ok(tags)
}

fn load_therapy_responses(path: &Path) -> Result<Vec<TherapyResponseClass>, StageError> {
    let lines = read_tsv(path)?;
    let mut responses = Vec::new();
    for row in lines.iter().skip(1) {
        let value = row
            .get(9)
            .ok_or_else(|| StageError::Format("missing response_class".to_string()))?;
        responses.push(parse_response(value, path)?);
    }
    Ok(responses)
}

fn parse_class(value: &str, path: &Path) -> Result<LysosomeDependencyClass, StageError> {
    match value {
        "LysosomeDependent" => Ok(LysosomeDependencyClass::LysosomeDependent),
        "AutophagyAdaptive" => Ok(LysosomeDependencyClass::AutophagyAdaptive),
        "StalledAutophagy" => Ok(LysosomeDependencyClass::StalledAutophagy),
        "LysosomeOverloaded" => Ok(LysosomeDependencyClass::LysosomeOverloaded),
        "EnergyRecyclingDependent" => Ok(LysosomeDependencyClass::EnergyRecyclingDependent),
        "NonLysosomal" => Ok(LysosomeDependencyClass::NonLysosomal),
        "Unclassified" => Ok(LysosomeDependencyClass::Unclassified),
        _ => Err(StageError::Format(format!(
            "unknown class {} in {}",
            value,
            path.display()
        ))),
    }
}

fn parse_vulnerability(value: &str, path: &Path) -> Result<DrugVulnerabilityTag, StageError> {
    match value {
        "AutophagyInhibitionSensitive" => Ok(DrugVulnerabilityTag::AutophagyInhibitionSensitive),
        "LysosomeAcidificationSensitive" => {
            Ok(DrugVulnerabilityTag::LysosomeAcidificationSensitive)
        }
        "VATPaseInhibitionSensitive" => Ok(DrugVulnerabilityTag::VATPaseInhibitionSensitive),
        "CathepsinInhibitionSensitive" => Ok(DrugVulnerabilityTag::CathepsinInhibitionSensitive),
        "MitophagyDisruptionSensitive" => Ok(DrugVulnerabilityTag::MitophagyDisruptionSensitive),
        "EnergyStressAmplificationSensitive" => {
            Ok(DrugVulnerabilityTag::EnergyStressAmplificationSensitive)
        }
        "LysosomalDestabilizationSensitive" => {
            Ok(DrugVulnerabilityTag::LysosomalDestabilizationSensitive)
        }
        "LowVulnerability" => Ok(DrugVulnerabilityTag::LowVulnerability),
        _ => Err(StageError::Format(format!(
            "unknown vulnerability tag {} in {}",
            value,
            path.display()
        ))),
    }
}

fn parse_response(value: &str, path: &Path) -> Result<TherapyResponseClass, StageError> {
    match value {
        "CYTOTOXIC_RESPONSE" => Ok(TherapyResponseClass::CytotoxicResponse),
        "ADAPTIVE_SURVIVAL" => Ok(TherapyResponseClass::AdaptiveSurvival),
        "LYSO_ESCAPE" => Ok(TherapyResponseClass::LysoEscape),
        "DAMAGE_ACCUMULATION" => Ok(TherapyResponseClass::DamageAccumulation),
        "NO_RESPONSE" => Ok(TherapyResponseClass::NoResponse),
        _ => Err(StageError::Format(format!(
            "unknown therapy response {} in {}",
            value,
            path.display()
        ))),
    }
}

fn parse_f32(value: &str, path: &Path) -> Result<f32, StageError> {
    value
        .parse::<f32>()
        .map_err(|_| StageError::Format(format!("invalid float {} in {}", value, path.display())))
}

fn read_tsv(path: &Path) -> Result<Vec<Vec<String>>, StageError> {
    if !path.is_file() {
        return Err(StageError::Validation(format!(
            "missing required file {}",
            path.display()
        )));
    }
    let content = std::fs::read_to_string(path)?;
    let mut rows = Vec::new();
    for line in content.lines() {
        let cols = line.split('\t').map(|s| s.to_string()).collect::<Vec<_>>();
        rows.push(cols);
    }
    Ok(rows)
}
