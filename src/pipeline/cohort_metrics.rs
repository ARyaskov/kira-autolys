use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;

use crate::config::thresholds::Thresholds;
use crate::model::ctx::Ctx;
use crate::stage_error::StageError;

pub fn write_cohort_metric_summary(ctx: &Ctx, out_dir: &Path) -> Result<(), StageError> {
    let mut summary: BTreeMap<&'static str, MetricSummary> = BTreeMap::new();
    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);

    if let Some(metrics) = ctx.survival.as_ref() {
        summary.insert("SSM", summarize_metric(&metrics.ssm, thresholds.ssm_high));
    }
    if let Some(metrics) = ctx.autophagy.as_ref() {
        summary.insert("AFP", summarize_metric(&metrics.afp, thresholds.afp_high));
    }
    if let Some(metrics) = ctx.lysosome.as_ref() {
        summary.insert("LDS", summarize_metric(&metrics.lds, thresholds.lds_high));
    }
    if let Some(metrics) = ctx.lysosomal_damage.as_ref() {
        summary.insert("LDI", summarize_metric(&metrics.ldi, thresholds.ldi_high));
    }
    if let Some(metrics) = ctx.coupling.as_ref() {
        summary.insert(
            "LSI",
            summarize_metric(&metrics.locked_survival_index, thresholds.lsi_high),
        );
    }
    if let Some(metrics) = ctx.cross_organelle.as_ref() {
        summary.insert(
            "ERDI",
            summarize_metric(&metrics.energy_recycling_dependency, thresholds.erdi_high),
        );
    }

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
#[path = "../../tests/src_inline/pipeline/cohort_metrics.rs"]
mod tests;
