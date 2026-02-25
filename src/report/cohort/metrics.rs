use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;

use crate::model::cohort_view::CohortView;
use crate::stage_error::StageError;

pub fn write_metrics(view: &CohortView, out_dir: &Path) -> Result<(), StageError> {
    let mut summary: BTreeMap<&'static str, MetricSummary> = BTreeMap::new();
    summary.insert("SSM", summarize_metric(&view.ssm, view.thresholds.ssm_high));
    summary.insert("AFP", summarize_metric(&view.afp, view.thresholds.afp_high));
    summary.insert("LDS", summarize_metric(&view.lds, view.thresholds.lds_high));
    summary.insert("LDI", summarize_metric(&view.ldi, view.thresholds.ldi_high));
    summary.insert("LSI", summarize_metric(&view.lsi, view.thresholds.lsi_high));
    summary.insert(
        "ERDI",
        summarize_metric(&view.erdi, view.thresholds.erdi_high),
    );

    fs::create_dir_all(out_dir)?;
    let path = out_dir.join("cohort_metrics.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &summary)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;

    Ok(())
}

#[derive(Debug, serde::Serialize)]
struct MetricSummary {
    mean: f32,
    median: f32,
    p90: f32,
    fraction_high: f32,
}

fn summarize_metric(values: &[f32], threshold: f32) -> MetricSummary {
    if values.is_empty() {
        return MetricSummary {
            mean: 0.0,
            median: 0.0,
            p90: 0.0,
            fraction_high: 0.0,
        };
    }
    let mean = values.iter().sum::<f32>() / values.len() as f32;
    let median = percentile(values, 0.5);
    let p90 = percentile(values, 0.9);
    let high = values.iter().filter(|v| **v >= threshold).count() as f32;
    let fraction_high = high / values.len() as f32;
    MetricSummary {
        mean,
        median,
        p90,
        fraction_high,
    }
}

fn percentile(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = ((data.len() - 1) as f32 * p).round() as usize;
    data[idx]
}

#[cfg(test)]
#[path = "../../../tests/src_inline/report/cohort/metrics.rs"]
mod tests;
