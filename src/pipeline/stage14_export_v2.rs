use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::config::thresholds::Thresholds;
use crate::model::ctx::Ctx;
use crate::model::schema_v2 as schema;
use crate::stage_error::StageError;

const CORE_HEADER: &str = "id\tAFP\tLDS\tTFEB\tmTOR\tSSM\tz_AFP\tz_LDS\tz_PROLIF\tz_APOP";
const DAMAGE_HEADER: &str = "id\tLDI\tLMP\tstress\tcathepsin_membrane_imbalance";
const SELECTIVITY_HEADER: &str =
    "id\tmito_frac\taggre_frac\ter_frac\tferr_frac\tlipo_frac\tentropy";
const COUPLING_HEADER: &str = "id\tLSI\tampk_mtor_ratio\ttfeb_lyso_alignment";
const CROSS_HEADER: &str =
    "id\tERDI\tmitophagy_reliance\tmito_lyso_imbalance\tmito_consumption_risk";
const CLASS_HEADER: &str = "id\tlysosome_dependency_class";
const VULN_HEADER: &str = "id\tvulnerability_tags";
const THERAPY_HEADER: &str = "transition_id\tfrom_id\tto_id\tdelta_SSM\tdelta_AFP\tdelta_LDS\tdelta_LDI\tdelta_ERDI\tASI\tresponse_class";
const IO_BUF_CAPACITY: usize = 1024 * 1024;

pub fn run_stage14_export_v2(ctx: &Ctx, out_dir: &Path) -> Result<(), StageError> {
    let input_type = ctx
        .input_type
        .ok_or_else(|| StageError::Validation("input_type missing".to_string()))?;
    let mode = ctx
        .mode
        .ok_or_else(|| StageError::Validation("mode missing".to_string()))?;

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
    let survival_stats = ctx
        .survival_stats
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival stats missing".to_string()))?;
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;
    let selectivity = ctx
        .autophagy_selectivity
        .as_ref()
        .ok_or_else(|| StageError::Validation("selectivity metrics missing".to_string()))?;
    let coupling = ctx
        .coupling
        .as_ref()
        .ok_or_else(|| StageError::Validation("coupling metrics missing".to_string()))?;
    let cross = ctx
        .cross_organelle
        .as_ref()
        .ok_or_else(|| StageError::Validation("cross-organelle metrics missing".to_string()))?;
    let class_metrics = ctx
        .lysosome_class
        .as_ref()
        .ok_or_else(|| StageError::Validation("classification missing".to_string()))?;
    let vuln_metrics = ctx
        .drug_vulnerability
        .as_ref()
        .ok_or_else(|| StageError::Validation("vulnerability tags missing".to_string()))?;

    fs::create_dir_all(out_dir)?;

    write_core_scores(out_dir, ctx, autophagy, lysosome, regulatory, survival)?;
    write_damage(out_dir, ctx, damage)?;
    write_selectivity(out_dir, ctx, selectivity)?;
    write_coupling(out_dir, ctx, coupling)?;
    write_cross_organelle(out_dir, ctx, cross)?;
    write_classification(out_dir, ctx, class_metrics)?;
    write_vulnerabilities(out_dir, ctx, vuln_metrics)?;

    if let Some(therapy) = ctx.therapy_delta.as_ref() {
        write_therapy_delta(out_dir, ctx, therapy)?;
    }

    write_summary_v2(
        out_dir,
        ctx,
        input_type,
        mode,
        survival_stats,
        regulatory,
        autophagy,
        lysosome,
        survival,
        damage,
        selectivity,
        coupling,
        cross,
        class_metrics,
        vuln_metrics,
    )?;

    Ok(())
}

fn build_samples<'a>(ctx: &'a Ctx) -> Vec<schema::SampleInfo<'a>> {
    if !ctx.samples.is_empty() {
        return ctx
            .samples
            .iter()
            .map(|s| schema::SampleInfo {
                id: s.id.as_str(),
                sample_group: s.sample_group.as_str(),
                timepoint: s.timepoint,
            })
            .collect();
    }

    ctx.obs_ids
        .iter()
        .map(|id| schema::SampleInfo {
            id: id.as_str(),
            sample_group: id.as_str(),
            timepoint: 0,
        })
        .collect()
}

fn write_summary_v2(
    out_dir: &Path,
    ctx: &Ctx,
    input_type: crate::model::ctx::InputType,
    mode: crate::model::ctx::Mode,
    survival_stats: &crate::model::ctx::SurvivalStats,
    regulatory: &crate::model::ctx::RegulatoryMetrics,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
    damage: &crate::model::ctx::LysosomalDamageMetrics,
    selectivity: &crate::model::ctx::AutophagySelectivityMetrics,
    coupling: &crate::model::ctx::CouplingMetrics,
    cross: &crate::model::ctx::CrossOrganelleMetrics,
    class_metrics: &crate::model::ctx::LysosomeClassMetrics,
    vuln_metrics: &crate::model::ctx::DrugVulnerabilityMetrics,
) -> Result<(), StageError> {
    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);
    let tool = schema::ToolInfo {
        name: "kira-autolys",
        version: env!("CARGO_PKG_VERSION"),
        schema: "v2",
    };
    let input = schema::InputInfo {
        r#type: match input_type {
            crate::model::ctx::InputType::Mtx10x => "mtx",
            crate::model::ctx::InputType::DenseTsv => "dense_tsv",
            crate::model::ctx::InputType::H5ad => "h5ad",
        },
        mode: match mode {
            crate::model::ctx::Mode::Cell => "cell",
            crate::model::ctx::Mode::Sample => "sample",
        },
        n_obs: ctx.n_obs,
        samples: build_samples(ctx),
    };

    let dataset_stats = schema::DatasetStats {
        core: schema::CoreStats::from(survival_stats),
        regulatory: schema::RegulatoryStats::from(regulatory),
        _phantom: None,
    };

    let thresholds = schema::ThresholdsInfo::from(&thresholds);

    let mut observations = Vec::with_capacity(ctx.n_obs);
    for idx in 0..ctx.n_obs {
        observations.push(schema::Observation {
            id: ctx.obs_ids[idx].as_str(),
            autophagy: Some(schema::make_autophagy_block(autophagy, idx)),
            lysosome: Some(schema::make_lysosome_block(lysosome, idx)),
            regulatory: Some(schema::make_regulatory_block(regulatory, idx)),
            survival: Some(schema::make_survival_block(survival, idx)),
            damage: Some(schema::make_damage_block(damage, idx)),
            selectivity: Some(schema::make_selectivity_block(selectivity, idx)),
            coupling: Some(schema::make_coupling_block(coupling, idx)),
            cross_organelle: Some(schema::make_cross_organelle_block(cross, idx)),
            classification: Some(schema::make_classification_block(class_metrics, idx)),
            vulnerabilities: Some(schema::make_vulnerability_block(vuln_metrics, idx)),
        });
    }

    let summary = serde_json::json!({
        "tool": tool,
        "input": input,
        "thresholds": thresholds,
        "dataset_stats": dataset_stats,
        "observations": observations,
    });

    let path = out_dir.join("summary_v2.json");
    let writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    // Compact JSON is significantly faster for very large observation arrays.
    serde_json::to_writer(writer, &summary)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;
    Ok(())
}

fn write_core_scores(
    out_dir: &Path,
    ctx: &Ctx,
    autophagy: &crate::model::ctx::AutophagyMetrics,
    lysosome: &crate::model::ctx::LysosomeMetrics,
    regulatory: &crate::model::ctx::RegulatoryMetrics,
    survival: &crate::model::ctx::SurvivalMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("core_scores.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{CORE_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write_f32(&mut writer, autophagy.afp[idx])?;
        write_f32(&mut writer, lysosome.lds[idx])?;
        write_f32(&mut writer, regulatory.tfeb[idx])?;
        write_f32(&mut writer, regulatory.mtor[idx])?;
        write_f32(&mut writer, survival.ssm[idx])?;
        write_f32(&mut writer, survival.z_afp[idx])?;
        write_f32(&mut writer, survival.z_lds[idx])?;
        write_f32(&mut writer, survival.z_prolif[idx])?;
        write_f32(&mut writer, survival.z_apop[idx])?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_damage(
    out_dir: &Path,
    ctx: &Ctx,
    damage: &crate::model::ctx::LysosomalDamageMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("damage.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{DAMAGE_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write_f32(&mut writer, damage.ldi[idx])?;
        write_f32(&mut writer, damage.lmp[idx])?;
        write_f32(&mut writer, damage.stress[idx])?;
        write_f32(&mut writer, damage.cathepsin_membrane_imbalance[idx])?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_selectivity(
    out_dir: &Path,
    ctx: &Ctx,
    selectivity: &crate::model::ctx::AutophagySelectivityMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("selectivity.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{SELECTIVITY_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write_f32(&mut writer, selectivity.mito_frac[idx])?;
        write_f32(&mut writer, selectivity.aggre_frac[idx])?;
        write_f32(&mut writer, selectivity.er_frac[idx])?;
        write_f32(&mut writer, selectivity.ferr_frac[idx])?;
        write_f32(&mut writer, selectivity.lipo_frac[idx])?;
        write_f32(&mut writer, selectivity.entropy[idx])?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_coupling(
    out_dir: &Path,
    ctx: &Ctx,
    coupling: &crate::model::ctx::CouplingMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("coupling.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{COUPLING_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write_f32(&mut writer, coupling.locked_survival_index[idx])?;
        write_f32(&mut writer, coupling.ampk_mtor_ratio[idx])?;
        write_f32(&mut writer, coupling.tfeb_lyso_alignment[idx])?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_cross_organelle(
    out_dir: &Path,
    ctx: &Ctx,
    cross: &crate::model::ctx::CrossOrganelleMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("cross_organelle.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{CROSS_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write_f32(&mut writer, cross.energy_recycling_dependency[idx])?;
        write_f32(&mut writer, cross.mitophagy_reliance[idx])?;
        write_f32(&mut writer, cross.mito_lyso_imbalance[idx])?;
        write_f32(&mut writer, cross.mito_consumption_risk[idx])?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_classification(
    out_dir: &Path,
    ctx: &Ctx,
    class_metrics: &crate::model::ctx::LysosomeClassMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("classification.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{CLASS_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        write!(
            writer,
            "\t{}",
            schema::class_to_str(class_metrics.class[idx])
        )?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_vulnerabilities(
    out_dir: &Path,
    ctx: &Ctx,
    vuln_metrics: &crate::model::ctx::DrugVulnerabilityMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("vulnerabilities.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{VULN_HEADER}")?;

    for idx in 0..ctx.n_obs {
        write!(writer, "{}", ctx.obs_ids[idx])?;
        let tags: Vec<&'static str> = vuln_metrics.tags[idx]
            .iter()
            .copied()
            .map(schema::vulnerability_to_str)
            .collect();
        write!(writer, "\t{}", tags.join(","))?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_therapy_delta(
    out_dir: &Path,
    ctx: &Ctx,
    therapy: &crate::model::ctx::TherapyDeltaMetrics,
) -> Result<(), StageError> {
    let path = out_dir.join("therapy_delta.tsv");
    let mut writer = BufWriter::with_capacity(IO_BUF_CAPACITY, File::create(path)?);
    writeln!(writer, "{THERAPY_HEADER}")?;

    for idx in 0..therapy.delta_ssm.len() {
        let from_id = ctx.obs_ids[therapy.from_obs[idx]].as_str();
        let to_id = ctx.obs_ids[therapy.to_obs[idx]].as_str();
        let transition_id = format!(
            "{}:{}->{}",
            therapy.sample_group[idx], therapy.timepoint0[idx], therapy.timepoint1[idx]
        );
        write!(writer, "{}", transition_id)?;
        write!(writer, "\t{}\t{}", from_id, to_id)?;
        write_f32(&mut writer, therapy.delta_ssm[idx])?;
        write_f32(&mut writer, therapy.delta_afp[idx])?;
        write_f32(&mut writer, therapy.delta_lds[idx])?;
        write_f32(&mut writer, therapy.delta_ldi[idx])?;
        write_f32(&mut writer, therapy.delta_erdi[idx])?;
        write_f32(&mut writer, therapy.asi[idx])?;
        write!(
            writer,
            "\t{}",
            match therapy.response_class[idx] {
                crate::model::ctx::TherapyResponseClass::CytotoxicResponse => "CYTOTOXIC_RESPONSE",
                crate::model::ctx::TherapyResponseClass::AdaptiveSurvival => "ADAPTIVE_SURVIVAL",
                crate::model::ctx::TherapyResponseClass::LysoEscape => "LYSO_ESCAPE",
                crate::model::ctx::TherapyResponseClass::DamageAccumulation =>
                    "DAMAGE_ACCUMULATION",
                crate::model::ctx::TherapyResponseClass::NoResponse => "NO_RESPONSE",
            }
        )?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_f32<W: Write>(writer: &mut W, value: f32) -> Result<(), StageError> {
    write!(writer, "\t{:.6}", value)?;
    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage14_export_v2.rs"]
mod tests;
