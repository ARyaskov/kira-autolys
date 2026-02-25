use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use serde_json::json;

use crate::config::thresholds::Thresholds;
use crate::genesets::v1;
use crate::model::cli::RunMode;
use crate::model::ctx::{Ctx, SampleInfo};
use crate::stage_error::StageError;

const TSV_HEADER: &str = "id\tAFP\tAFP_initiation\tAFP_elongation\tAFP_degradation\tAFP_cargo\tAFP_stall\tLDS\tLDS_vatp\tLDS_prot\tLDS_mem\tglobal_load\tTFEBp\tmTORp\tTFEB_MTOR\tTFEB_MTOR_diff\tPROLIF\tAPOP_ready\tz_AFP\tz_LDS\tz_PROLIF\tz_APOP\tSSM\tflag_SSM_high\tflag_AFP_high\tflag_LDS_high\tflag_PROLIF_low\tflag_APOP_low";
const PIPELINE_AUTOLYS_HEADER: &str = "id\tbarcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\tlysosomal_load\tautophagy_flux_proxy\tcatabolic_bias\trecycling_efficiency\tstress_autophagy_index\tregime\tflags\tconfidence";
const PIPELINE_PANELS_HEADER: &str = "panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99";
const REGIMES: [&str; 6] = [
    "HomeostaticRecycling",
    "StressInducedAutophagy",
    "CatabolicOverdrive",
    "AutophagyImpaired",
    "DegenerativeLysosomal",
    "Unclassified",
];

pub fn run_stage6(ctx: &Ctx, out_dir: &Path, write_scores: bool) -> Result<(), StageError> {
    if matches!(ctx.run_mode, RunMode::Pipeline) {
        return run_stage6_pipeline(ctx, out_dir);
    }

    let input_type = ctx
        .input_type
        .ok_or_else(|| StageError::Validation("input_type missing".to_string()))?;
    let mode = ctx
        .mode
        .ok_or_else(|| StageError::Validation("mode missing".to_string()))?;

    let n_obs = ctx.n_obs;
    if ctx.obs_ids.len() != n_obs {
        return Err(StageError::Validation(
            "obs_ids length does not match n_obs".to_string(),
        ));
    }

    let autophagy = ctx
        .autophagy
        .as_ref()
        .ok_or_else(|| StageError::Validation("autophagy metrics missing".to_string()))?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    let regulatory = ctx
        .regulatory
        .as_ref()
        .ok_or_else(|| StageError::Validation("regulatory metrics missing".to_string()))?;
    let survival = ctx
        .survival
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival metrics missing".to_string()))?;
    let stats = ctx
        .survival_stats
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival stats missing".to_string()))?;

    fs::create_dir_all(out_dir)?;

    if write_scores {
        write_scores_tsv(out_dir, ctx, autophagy, lysosome, regulatory, survival)?;
    }
    write_panels_report(out_dir, ctx)?;
    write_metadata_json(out_dir, ctx, input_type, mode)?;
    write_summary_json(
        out_dir, ctx, autophagy, lysosome, regulatory, survival, stats,
    )?;

    let summary_text = build_stdout_summary(ctx, autophagy, lysosome, survival);
    print_summary(&summary_text);
    write_report(out_dir, &summary_text)?;
    write_report_summary(out_dir, ctx, autophagy, lysosome, survival)?;

    Ok(())
}

#[derive(Clone)]
struct PipelineRow {
    id: String,
    barcode: String,
    sample: String,
    condition: String,
    species: &'static str,
    libsize: f32,
    nnz: u32,
    expressed_genes: u32,
    lysosomal_load: f32,
    autophagy_flux_proxy: f32,
    catabolic_bias: f32,
    recycling_efficiency: f32,
    stress_autophagy_index: f32,
    regime: &'static str,
    flags: String,
    confidence: f32,
}

fn run_stage6_pipeline(ctx: &Ctx, out_dir: &Path) -> Result<(), StageError> {
    let autophagy = ctx
        .autophagy
        .as_ref()
        .ok_or_else(|| StageError::Validation("autophagy metrics missing".to_string()))?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    let survival = ctx
        .survival
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival metrics missing".to_string()))?;
    let normalized = ctx
        .normalized
        .as_ref()
        .ok_or_else(|| StageError::Validation("normalized expression missing".to_string()))?;

    if ctx.obs_ids.len() != ctx.n_obs
        || ctx.obs_libsize.len() != ctx.n_obs
        || ctx.obs_nnz.len() != ctx.n_obs
        || ctx.obs_expressed_genes.len() != ctx.n_obs
    {
        return Err(StageError::Validation(
            "stage1 per-observation QC vectors are missing or length-mismatched".to_string(),
        ));
    }

    fs::create_dir_all(out_dir)?;

    let sample_map = build_sample_map(ctx);
    let species = infer_species(&ctx.gene_symbols);

    let ldi_values = if let Some(damage) = ctx.lysosomal_damage.as_ref() {
        damage.ldi.clone()
    } else {
        vec![0.0_f32; ctx.n_obs]
    };

    let raw_catabolic: Vec<f32> = (0..ctx.n_obs)
        .map(|i| lysosome.lds[i] + autophagy.afp[i] - survival.prolif[i])
        .collect();
    let raw_recycling: Vec<f32> = (0..ctx.n_obs)
        .map(|i| autophagy.afp[i] - autophagy.stall[i] - ldi_values[i])
        .collect();

    let lyso_norm = quantile_normalizer(&lysosome.lds, 0.10, 0.90);
    let afp_norm = quantile_normalizer(&autophagy.afp, 0.10, 0.90);
    let stress_norm = quantile_normalizer(&survival.ssm, 0.10, 0.90);
    let cat_norm = quantile_normalizer(&raw_catabolic, 0.10, 0.90);
    let rec_norm = quantile_normalizer(&raw_recycling, 0.10, 0.90);

    let mut rows = Vec::with_capacity(ctx.n_obs);
    let mut low_conf_count = 0usize;
    let mut low_lyso_count = 0usize;

    for idx in 0..ctx.n_obs {
        let (sample, condition) = sample_map
            .get(ctx.obs_ids[idx].as_str())
            .cloned()
            .unwrap_or_else(|| (ctx.obs_ids[idx].clone(), "NA".to_string()));

        let lysosomal_load = lyso_norm(lysosome.lds[idx]);
        let autophagy_flux_proxy = afp_norm(autophagy.afp[idx]);
        let catabolic_bias = cat_norm(raw_catabolic[idx]);
        let recycling_efficiency = rec_norm(raw_recycling[idx]);
        let stress_autophagy_index = stress_norm(survival.ssm[idx]);

        let mut flags = Vec::new();
        let mut confidence = 1.0_f32;
        if ctx.obs_libsize[idx] <= 0.0 || ctx.obs_nnz[idx] < 20 {
            flags.push("LOW_CONFIDENCE");
            confidence -= 0.45;
            low_conf_count += 1;
        }
        if lysosomal_load < 0.15 {
            flags.push("LOW_LYSO_SIGNAL");
            confidence -= 0.35;
            low_lyso_count += 1;
        }
        if confidence < 0.0 {
            confidence = 0.0;
        }

        let regime = infer_regime(
            lysosomal_load,
            autophagy_flux_proxy,
            catabolic_bias,
            recycling_efficiency,
            stress_autophagy_index,
            confidence,
        );
        let flags = flags.join(",");

        rows.push(PipelineRow {
            id: ctx.obs_ids[idx].clone(),
            barcode: ctx.obs_ids[idx].clone(),
            sample,
            condition,
            species,
            libsize: ctx.obs_libsize[idx],
            nnz: ctx.obs_nnz[idx],
            expressed_genes: ctx.obs_expressed_genes[idx],
            lysosomal_load,
            autophagy_flux_proxy,
            catabolic_bias,
            recycling_efficiency,
            stress_autophagy_index,
            regime,
            flags,
            confidence,
        });
    }

    rows.sort_by(|a, b| a.id.cmp(&b.id).then_with(|| a.barcode.cmp(&b.barcode)));

    write_autolys_tsv(out_dir, &rows)?;
    let coverage_median = write_panels_report_pipeline(out_dir, ctx, normalized.as_ref())?;
    write_summary_pipeline_json(
        out_dir,
        ctx,
        &rows,
        low_conf_count,
        low_lyso_count,
        coverage_median,
    )?;
    write_pipeline_step_json(out_dir)?;
    Ok(())
}

fn write_autolys_tsv(out_dir: &Path, rows: &[PipelineRow]) -> Result<(), StageError> {
    let mut writer = BufWriter::new(File::create(out_dir.join("autolys.tsv"))?);
    writeln!(writer, "{PIPELINE_AUTOLYS_HEADER}")?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}",
            row.id,
            row.barcode,
            row.sample,
            row.condition,
            row.species,
            row.libsize,
            row.nnz,
            row.expressed_genes,
            row.lysosomal_load,
            row.autophagy_flux_proxy,
            row.catabolic_bias,
            row.recycling_efficiency,
            row.stress_autophagy_index,
            row.regime,
            row.flags,
            row.confidence
        )?;
    }
    Ok(())
}

fn write_panels_report_pipeline(
    out_dir: &Path,
    ctx: &Ctx,
    normalized: &dyn crate::model::ctx::NormalizedExpr,
) -> Result<f32, StageError> {
    let defs = v1::gene_sets();
    let mut def_size_by_name: BTreeMap<&str, usize> = BTreeMap::new();
    for def in defs {
        def_size_by_name.insert(def.name, def.genes.len());
    }
    let n_panels = ctx.gene_sets.len();
    let mut coverage_vals = (0..n_panels)
        .map(|_| Vec::with_capacity(ctx.n_obs))
        .collect::<Vec<_>>();
    let mut sum_vals = (0..n_panels)
        .map(|_| Vec::with_capacity(ctx.n_obs))
        .collect::<Vec<_>>();
    let missing_map: BTreeMap<&str, &Vec<String>> = ctx
        .missing_genes
        .iter()
        .map(|m| (m.panel.as_str(), &m.missing))
        .collect();

    for obs_idx in 0..ctx.n_obs {
        let mut counts = vec![0u32; n_panels];
        let mut sums = vec![0f32; n_panels];
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, norm_value| {
            let mut mask = ctx.gene_panel_mask[gene_idx];
            while mask != 0 {
                let panel_idx = mask.trailing_zeros() as usize;
                if panel_idx < n_panels {
                    counts[panel_idx] += 1;
                    sums[panel_idx] += norm_value;
                }
                mask &= mask - 1;
            }
        });
        for panel_idx in 0..n_panels {
            let panel_name = ctx.gene_sets[panel_idx].name.as_str();
            let missing_count = missing_map.get(panel_name).map(|v| v.len()).unwrap_or(0);
            let defined = def_size_by_name
                .get(panel_name)
                .copied()
                .unwrap_or(ctx.gene_sets[panel_idx].indices.len() + missing_count)
                as f32;
            let coverage = if defined > 0.0 {
                counts[panel_idx] as f32 / defined
            } else {
                0.0
            };
            coverage_vals[panel_idx].push(coverage);
            sum_vals[panel_idx].push(sums[panel_idx]);
        }
    }

    let mut writer = BufWriter::new(File::create(out_dir.join("panels_report.tsv"))?);
    writeln!(writer, "{PIPELINE_PANELS_HEADER}")?;
    let mut panel_coverage_medians = Vec::with_capacity(n_panels);
    for panel_idx in 0..n_panels {
        let panel_name = ctx.gene_sets[panel_idx].name.as_str();
        let panel_group = panel_name.split('_').next().unwrap_or("UNKNOWN");
        let panel_size_defined = def_size_by_name
            .get(panel_name)
            .copied()
            .unwrap_or_else(|| {
                let missing_count = missing_map.get(panel_name).map(|v| v.len()).unwrap_or(0);
                ctx.gene_sets[panel_idx].indices.len() + missing_count
            });
        let panel_size_mappable = ctx.gene_sets[panel_idx].indices.len();
        let missing = missing_map
            .get(panel_name)
            .map(|v| {
                let mut x = (*v).clone();
                x.sort();
                x.join(",")
            })
            .unwrap_or_default();
        let cov = &coverage_vals[panel_idx];
        let sum = &sum_vals[panel_idx];
        let coverage_median = percentile(cov, 0.50);
        panel_coverage_medians.push(coverage_median);
        writeln!(
            writer,
            "P{:03}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            panel_idx + 1,
            panel_name,
            panel_group,
            panel_size_defined,
            panel_size_mappable,
            missing,
            coverage_median,
            percentile(cov, 0.10),
            percentile(sum, 0.50),
            percentile(sum, 0.90),
            percentile(sum, 0.99),
        )?;
    }
    Ok(percentile(&panel_coverage_medians, 0.50))
}

fn write_summary_pipeline_json(
    out_dir: &Path,
    ctx: &Ctx,
    rows: &[PipelineRow],
    low_conf_count: usize,
    low_lyso_count: usize,
    coverage_median: f32,
) -> Result<(), StageError> {
    let mut regime_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for regime in REGIMES {
        regime_counts.insert(regime, 0);
    }
    let mut lyso = Vec::with_capacity(rows.len());
    let mut afp = Vec::with_capacity(rows.len());
    let mut stress = Vec::with_capacity(rows.len());
    for row in rows {
        *regime_counts.get_mut(row.regime).unwrap() += 1;
        lyso.push(row.lysosomal_load);
        afp.push(row.autophagy_flux_proxy);
        stress.push(row.stress_autophagy_index);
    }
    let n = rows.len().max(1) as f32;
    let mut regime_fractions: BTreeMap<&str, f32> = BTreeMap::new();
    for (k, v) in &regime_counts {
        regime_fractions.insert(k, round6(*v as f32 / n));
    }
    let n_samples = rows
        .iter()
        .map(|row| row.sample.as_str())
        .collect::<BTreeSet<_>>()
        .len();

    let summary = json!({
        "tool": "kira-autolys",
        "meta": {
            "version": env!("CARGO_PKG_VERSION"),
            "simd": simd_label(),
        },
        "input": {
            "mode": "pipeline",
            "n_cells": ctx.n_obs,
            "n_samples": n_samples,
            "species": infer_species(&ctx.gene_symbols),
        },
        "distributions": {
            "lysosomal_load": {
                "median": round6(percentile(&lyso, 0.50)),
                "p90": round6(percentile(&lyso, 0.90)),
                "p99": round6(percentile(&lyso, 0.99)),
            },
            "autophagy_flux_proxy": {
                "median": round6(percentile(&afp, 0.50)),
                "p90": round6(percentile(&afp, 0.90)),
                "p99": round6(percentile(&afp, 0.99)),
            },
            "stress_autophagy_index": {
                "median": round6(percentile(&stress, 0.50)),
                "p90": round6(percentile(&stress, 0.90)),
                "p99": round6(percentile(&stress, 0.99)),
            },
        },
        "regimes": regime_fractions,
        "regime_counts": regime_counts,
        "qc": {
            "coverage_median": round6(coverage_median),
            "low_confidence_fraction": round6(low_conf_count as f32 / n),
            "low_lysosomal_signal_fraction": round6(low_lyso_count as f32 / n),
        }
    });
    serde_json::to_writer_pretty(File::create(out_dir.join("summary.json"))?, &summary)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;
    Ok(())
}

fn write_pipeline_step_json(out_dir: &Path) -> Result<(), StageError> {
    let step = json!({
        "tool": "kira-autolys",
        "mode": "pipeline",
        "artifacts": {
            "summary": "summary.json",
            "primary_metrics": "autolys.tsv",
        },
        "cell_metrics": {
            "file": "autolys.tsv",
            "id_column": "id",
            "regime_column": "regime",
            "confidence_column": "confidence",
            "flag_column": "flags",
        },
        "regimes": {
            "source": "autolys.tsv",
            "column": "regime",
        },
    });
    serde_json::to_writer_pretty(File::create(out_dir.join("pipeline_step.json"))?, &step)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;
    Ok(())
}

fn build_sample_map(ctx: &Ctx) -> BTreeMap<&str, (String, String)> {
    let mut map = BTreeMap::new();
    for s in &ctx.samples {
        map.insert(
            s.id.as_str(),
            (
                s.sample_group.clone(),
                s.treatment.clone().unwrap_or_else(|| "NA".to_string()),
            ),
        );
    }
    map
}

fn quantile_normalizer<'a>(values: &'a [f32], p_low: f32, p_high: f32) -> impl Fn(f32) -> f32 + 'a {
    let lo = percentile(values, p_low);
    let hi = percentile(values, p_high);
    move |v: f32| {
        if hi <= lo {
            return 0.5;
        }
        ((v - lo) / (hi - lo)).clamp(0.0, 1.0)
    }
}

fn percentile(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = data.len();
    let pos = ((n - 1) as f32 * p.clamp(0.0, 1.0)).round() as usize;
    data[pos]
}

fn round6(value: f32) -> f32 {
    (value * 1_000_000.0).round() / 1_000_000.0
}

fn infer_regime(
    lysosomal_load: f32,
    autophagy_flux_proxy: f32,
    catabolic_bias: f32,
    recycling_efficiency: f32,
    stress_autophagy_index: f32,
    confidence: f32,
) -> &'static str {
    if confidence < 0.20 {
        return "Unclassified";
    }
    if stress_autophagy_index >= 0.75 && catabolic_bias >= 0.75 {
        return "CatabolicOverdrive";
    }
    if stress_autophagy_index >= 0.60 && recycling_efficiency < 0.30 {
        return "AutophagyImpaired";
    }
    if lysosomal_load < 0.20 && recycling_efficiency < 0.30 {
        return "DegenerativeLysosomal";
    }
    if stress_autophagy_index >= 0.60 || autophagy_flux_proxy >= 0.70 {
        return "StressInducedAutophagy";
    }
    if recycling_efficiency >= 0.40 {
        return "HomeostaticRecycling";
    }
    "Unclassified"
}

fn infer_species(_gene_symbols: &[String]) -> &'static str {
    "unknown"
}

fn simd_label() -> &'static str {
    #[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
    {
        return "avx2";
    }
    #[cfg(all(feature = "simd", target_arch = "aarch64", target_feature = "neon"))]
    {
        return "neon";
    }
    #[cfg(not(any(
        all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"),
        all(feature = "simd", target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        "scalar"
    }
}

fn write_scores_tsv(
    out_dir: &Path,
    ctx: &Ctx,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    regulatory: &crate::model::ctx::RegulatoryMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("scores.tsv");
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "{TSV_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write_f32(&mut writer, autophagy.afp[idx])?;
        write_f32(&mut writer, autophagy.initiation[idx])?;
        write_f32(&mut writer, autophagy.elongation[idx])?;
        write_f32(&mut writer, autophagy.degradation[idx])?;
        write_f32(&mut writer, autophagy.cargo[idx])?;
        write_f32(&mut writer, autophagy.stall[idx])?;
        write_f32(&mut writer, lysosome.lds[idx])?;
        write_f32(&mut writer, lysosome.vatp[idx])?;
        write_f32(&mut writer, lysosome.prot[idx])?;
        write_f32(&mut writer, lysosome.mem[idx])?;
        write_f32(&mut writer, lysosome.global_load[idx])?;
        write_f32(&mut writer, regulatory.tfeb[idx])?;
        write_f32(&mut writer, regulatory.mtor[idx])?;
        write_f32(&mut writer, regulatory.tfeb_mtor_ratio[idx])?;
        write_f32(&mut writer, regulatory.tfeb_mtor_diff[idx])?;
        write_f32(&mut writer, survival.prolif[idx])?;
        write_f32(&mut writer, survival.apop_ready[idx])?;
        write_f32(&mut writer, survival.z_afp[idx])?;
        write_f32(&mut writer, survival.z_lds[idx])?;
        write_f32(&mut writer, survival.z_prolif[idx])?;
        write_f32(&mut writer, survival.z_apop[idx])?;
        write_f32(&mut writer, survival.ssm[idx])?;
        write_bool(&mut writer, survival.flag_ssm_high[idx])?;
        write_bool(&mut writer, survival.flag_afp_high[idx])?;
        write_bool(&mut writer, survival.flag_lds_high[idx])?;
        write_bool(&mut writer, survival.flag_prolif_low[idx])?;
        write_bool(&mut writer, survival.flag_apop_low[idx])?;
        writeln!(writer)?;
    }

    Ok(())
}

fn write_f32<W: Write>(writer: &mut W, value: f32) -> Result<(), StageError> {
    write!(writer, "\t{:.6}", value)?;
    Ok(())
}

fn write_bool<W: Write>(writer: &mut W, value: bool) -> Result<(), StageError> {
    write!(writer, "\t{}", if value { "true" } else { "false" })?;
    Ok(())
}

fn write_panels_report(out_dir: &Path, ctx: &Ctx) -> Result<(), StageError> {
    let path = out_dir.join("panels_report.tsv");
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "panel_name\tpresent_genes\tmissing_genes")?;

    let mut missing_map: BTreeMap<&str, &Vec<String>> = BTreeMap::new();
    for entry in &ctx.missing_genes {
        missing_map.insert(entry.panel.as_str(), &entry.missing);
    }

    for panel in &ctx.gene_sets {
        let mut present: Vec<String> = panel
            .indices
            .iter()
            .filter_map(|&idx| ctx.gene_symbols.get(idx).cloned())
            .collect();
        present.sort();
        let missing = missing_map
            .get(panel.name.as_str())
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        writeln!(
            writer,
            "{}\t{}\t{}",
            panel.name,
            present.join(","),
            missing.join(",")
        )?;
    }

    Ok(())
}

fn write_summary_json(
    out_dir: &Path,
    ctx: &Ctx,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    regulatory: &crate::model::ctx::RegulatoryMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
    stats: &crate::model::ctx::SurvivalStats,
) -> Result<(), StageError> {
    let missing_map = build_missing_map(ctx);
    let thresholds = ctx.thresholds.unwrap_or_default();
    let normalization = ctx
        .normalization
        .clone()
        .unwrap_or_else(|| serde_json::Value::Null);

    let ssm_high_count = survival.flag_ssm_high.iter().filter(|v| **v).count();
    let ssm_high_fraction = if ctx.n_obs == 0 {
        0.0
    } else {
        ssm_high_count as f32 / ctx.n_obs as f32
    };

    let examples = top_examples(ctx, autophagy, lysosome, survival, 3);

    let json = json!({
        "tool": {
            "name": "kira-autolys",
            "version": env!("CARGO_PKG_VERSION"),
            "schema": "v1",
        },
        "normalization": normalization,
        "gene_sets": {
            "version": "v1",
            "missing": missing_map,
        },
        "dataset_stats": {
            "means": {
                "AFP": stats.mean_afp,
                "LDS": stats.mean_lds,
                "PROLIF": stats.mean_prolif,
                "APOP": stats.mean_apop,
            },
            "stds": {
                "AFP": stats.std_afp,
                "LDS": stats.std_lds,
                "PROLIF": stats.std_prolif,
                "APOP": stats.std_apop,
            },
            "medians": {
                "TFEB": regulatory.median_tfeb,
                "MTOR": regulatory.median_mtor,
            },
        },
        "thresholds": {
            "SSM_high": thresholds.ssm_high,
            "AFP_high": thresholds.afp_high,
            "LDS_high": thresholds.lds_high,
            "PROLIF_low": thresholds.prolif_low,
            "APOP_low": thresholds.apop_low,
        },
        "top_findings": {
            "ssm_high_count": ssm_high_count,
            "ssm_high_fraction": ssm_high_fraction,
            "examples": examples,
        }
    });

    let path = out_dir.join("summary.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &json)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;
    Ok(())
}

fn write_metadata_json(
    out_dir: &Path,
    ctx: &Ctx,
    input_type: crate::model::ctx::InputType,
    mode: crate::model::ctx::Mode,
) -> Result<(), StageError> {
    let samples = build_samples(ctx);
    let normalization = ctx
        .normalization
        .clone()
        .unwrap_or_else(|| serde_json::Value::Null);
    let json = json!({
        "tool": {
            "name": "kira-autolys",
            "version": env!("CARGO_PKG_VERSION"),
            "schema": "v1",
        },
        "input": {
            "type": match input_type {
                crate::model::ctx::InputType::Mtx10x => "mtx",
                crate::model::ctx::InputType::DenseTsv => "dense_tsv",
                crate::model::ctx::InputType::H5ad => "h5ad",
            },
            "mode": match mode { crate::model::ctx::Mode::Cell => "cell", crate::model::ctx::Mode::Sample => "sample" },
            "n_obs": ctx.n_obs,
            "n_vars": ctx.n_vars,
            "samples": samples,
        },
        "normalization": normalization,
    });

    let path = out_dir.join("metadata.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &json)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;
    Ok(())
}

fn build_samples(ctx: &Ctx) -> Vec<serde_json::Value> {
    if !ctx.samples.is_empty() {
        return ctx.samples.iter().map(|s| sample_to_json(s)).collect();
    }
    ctx.obs_ids
        .iter()
        .map(|id| {
            sample_to_json(&SampleInfo {
                id: id.clone(),
                sample_group: id.clone(),
                treatment: None,
                timepoint: 0,
                timepoint_label: None,
            })
        })
        .collect()
}

fn sample_to_json(sample: &SampleInfo) -> serde_json::Value {
    let timepoint = sample
        .timepoint_label
        .as_ref()
        .map(|v| serde_json::Value::String(v.clone()))
        .unwrap_or_else(|| serde_json::Value::from(sample.timepoint));
    json!({
        "id": sample.id,
        "treatment": sample.treatment,
        "timepoint": timepoint,
    })
}

fn build_missing_map(ctx: &Ctx) -> serde_json::Value {
    let mut map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for entry in &ctx.missing_genes {
        let mut missing = entry.missing.clone();
        missing.sort();
        map.insert(entry.panel.clone(), missing);
    }
    serde_json::to_value(map).unwrap_or_else(|_| serde_json::Value::Null)
}

fn top_examples(
    ctx: &Ctx,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
    limit: usize,
) -> Vec<serde_json::Value> {
    let mut rows: Vec<(usize, f32, &str)> = ctx
        .obs_ids
        .iter()
        .enumerate()
        .map(|(idx, id)| (idx, survival.ssm[idx], id.as_str()))
        .collect();
    rows.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.2.cmp(b.2))
    });

    rows.into_iter()
        .take(limit)
        .map(|(idx, _, id)| {
            let mut flags = Vec::new();
            if survival.flag_ssm_high[idx] {
                flags.push("SSM_HIGH");
            }
            if survival.flag_afp_high[idx] {
                flags.push("AFP_HIGH");
            }
            if survival.flag_lds_high[idx] {
                flags.push("LDS_HIGH");
            }
            if survival.flag_prolif_low[idx] {
                flags.push("PROLIF_LOW");
            }
            if survival.flag_apop_low[idx] {
                flags.push("APOP_LOW");
            }
            json!({
                "id": id,
                "SSM": survival.ssm[idx],
                "AFP": autophagy.afp[idx],
                "LDS": lysosome.lds[idx],
                "flags": flags,
            })
        })
        .collect()
}

fn build_stdout_summary(
    ctx: &Ctx,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
) -> String {
    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);
    let ssm_high_count = survival.flag_ssm_high.iter().filter(|v| **v).count();
    let ssm_high_fraction = if ctx.n_obs == 0 {
        0.0
    } else {
        ssm_high_count as f32 / ctx.n_obs as f32 * 100.0
    };

    let median_afp = median(&autophagy.afp);
    let median_lds = median(&lysosome.lds);
    let median_ssm = median(&survival.ssm);

    let top_examples = top_examples(ctx, autophagy, lysosome, survival, 3);
    let top_ids: Vec<&str> = top_examples
        .iter()
        .filter_map(|v| v.get("id").and_then(|v| v.as_str()))
        .collect();

    let mode_label = match ctx.mode {
        Some(crate::model::ctx::Mode::Cell) => "cell",
        Some(crate::model::ctx::Mode::Sample) => "sample",
        None => "unknown",
    };
    let mut lines = Vec::new();
    lines.push(format!("kira-autolys v{}", env!("CARGO_PKG_VERSION")));
    lines.push(format!("Input: {} obs, mode={}", ctx.n_obs, mode_label));
    lines.push(format!(
        "SSM_high: {} obs ({:.1}%)",
        ssm_high_count, ssm_high_fraction
    ));
    lines.push(format!(
        "Median AFP={:.2}  LDS={:.2}  SSM={:.2}",
        median_afp, median_lds, median_ssm
    ));
    if !top_ids.is_empty() {
        lines.push(format!("Top SSM obs: {}", top_ids.join(", ")));
    }

    let (stable, overloaded, damaged) = lysosome_state_counts(ctx, lysosome, &thresholds);
    let total = ctx.n_obs.max(1) as f32;
    lines.push(format!(
        "Lysosome state: stable {:.1}%, overloaded {:.1}%, damaged {:.1}%",
        stable as f32 / total * 100.0,
        overloaded as f32 / total * 100.0,
        damaged as f32 / total * 100.0
    ));

    if let Some(classes) = ctx.lysosome_class.as_ref() {
        lines.push("Class distribution:".to_string());
        for (label, count) in class_counts(classes).into_iter() {
            lines.push(format!("  {label}: {count}"));
        }
    } else {
        lines.push("Class distribution: not computed (Stage 12 not run)".to_string());
    }

    lines.join("\n")
}

fn median(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = data.len() / 2;
    if data.len() % 2 == 1 {
        data[mid]
    } else {
        0.5 * (data[mid - 1] + data[mid])
    }
}

fn print_summary(summary: &str) {
    for line in summary.lines() {
        println!("{line}");
    }
}

fn write_report(out_dir: &Path, summary: &str) -> Result<(), StageError> {
    let path = out_dir.join("report.txt");
    let mut file = BufWriter::new(File::create(path)?);
    writeln!(file, "{summary}")?;
    Ok(())
}

fn write_report_summary(
    out_dir: &Path,
    ctx: &Ctx,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
) -> Result<(), StageError> {
    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);
    let (stable, overloaded, damaged) = lysosome_state_counts(ctx, lysosome, &thresholds);
    let total = ctx.n_obs.max(1) as f32;

    let class_dist = ctx
        .lysosome_class
        .as_ref()
        .map(class_counts)
        .unwrap_or_default();

    let json = json!({
        "tool": {
            "name": "kira-autolys",
            "version": env!("CARGO_PKG_VERSION"),
            "schema": "report-v1",
        },
        "input": {
            "n_obs": ctx.n_obs,
            "mode": ctx.mode.map(|m| match m {
                crate::model::ctx::Mode::Cell => "cell",
                crate::model::ctx::Mode::Sample => "sample",
            }).unwrap_or("unknown"),
        },
        "summary": {
            "median_afp": median(&autophagy.afp),
            "median_lds": median(&lysosome.lds),
            "median_ssm": median(&survival.ssm),
            "ssm_high_count": survival.flag_ssm_high.iter().filter(|v| **v).count(),
        },
        "lysosome_state": {
            "stable_count": stable,
            "overloaded_count": overloaded,
            "damaged_count": damaged,
            "stable_fraction": stable as f32 / total,
            "overloaded_fraction": overloaded as f32 / total,
            "damaged_fraction": damaged as f32 / total,
        },
        "class_distribution": class_dist,
    });

    let path = out_dir.join("report_summary.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &json)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;
    Ok(())
}

fn lysosome_state_counts(
    ctx: &Ctx,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    thresholds: &Thresholds,
) -> (usize, usize, usize) {
    let mut stable = 0usize;
    let mut overloaded = 0usize;
    let mut damaged = 0usize;
    if let Some(damage) = ctx.lysosomal_damage.as_ref() {
        for i in 0..ctx.n_obs {
            let lds = lysosome.lds[i];
            let ldi = damage.ldi[i];
            if ldi >= thresholds.ldi_high {
                damaged += 1;
            } else if lds >= thresholds.lds_high {
                overloaded += 1;
            } else {
                stable += 1;
            }
        }
    } else {
        stable = ctx.n_obs;
    }
    (stable, overloaded, damaged)
}

fn class_counts(classes: &crate::model::ctx::LysosomeClassMetrics) -> Vec<(&'static str, usize)> {
    use crate::model::ctx::LysosomeDependencyClass as C;
    let mut counts = vec![
        ("LysosomeDependent", 0usize),
        ("AutophagyAdaptive", 0),
        ("StalledAutophagy", 0),
        ("LysosomeOverloaded", 0),
        ("EnergyRecyclingDependent", 0),
        ("NonLysosomal", 0),
        ("Unclassified", 0),
    ];
    for class in &classes.class {
        let idx = match class {
            C::LysosomeDependent => 0,
            C::AutophagyAdaptive => 1,
            C::StalledAutophagy => 2,
            C::LysosomeOverloaded => 3,
            C::EnergyRecyclingDependent => 4,
            C::NonLysosomal => 5,
            C::Unclassified => 6,
        };
        counts[idx].1 += 1;
    }
    counts
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage6_export.rs"]
mod tests;
