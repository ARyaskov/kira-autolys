use std::collections::BTreeMap;
use tracing::warn;

use crate::model::ctx::{Ctx, TherapyDeltaMetrics, TherapyResponseClass};
use crate::stage_error::StageError;

pub fn run_stage11(ctx: &mut Ctx) -> Result<(), StageError> {
    let n_obs = ctx.n_obs;
    if ctx.obs_ids.len() != n_obs {
        return Err(StageError::Validation(
            "obs_ids length does not match n_obs".to_string(),
        ));
    }

    if ctx.samples.len() != n_obs {
        return Err(StageError::Validation(
            "sample metadata missing or incomplete".to_string(),
        ));
    }

    let survival = ctx
        .survival
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival metrics missing".to_string()))?;
    let autophagy = ctx
        .autophagy
        .as_ref()
        .ok_or_else(|| StageError::Validation("autophagy metrics missing".to_string()))?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;
    let cross = ctx
        .cross_organelle
        .as_ref()
        .ok_or_else(|| StageError::Validation("cross-organelle metrics missing".to_string()))?;

    let mut groups: BTreeMap<String, Vec<(i32, usize)>> = BTreeMap::new();
    for (idx, sample) in ctx.samples.iter().enumerate() {
        if sample.sample_group.is_empty() {
            return Err(StageError::Validation(
                "sample_group missing in metadata".to_string(),
            ));
        }
        groups
            .entry(sample.sample_group.clone())
            .or_default()
            .push((sample.timepoint, idx));
    }

    let mut metrics = TherapyDeltaMetrics {
        delta_ssm: Vec::new(),
        delta_afp: Vec::new(),
        delta_lds: Vec::new(),
        delta_ldi: Vec::new(),
        delta_erdi: Vec::new(),
        delta_prolif: Vec::new(),
        asi: Vec::new(),
        response_class: Vec::new(),
        from_obs: Vec::new(),
        to_obs: Vec::new(),
        sample_group: Vec::new(),
        timepoint0: Vec::new(),
        timepoint1: Vec::new(),
    };

    for (group, mut entries) in groups {
        if entries.len() < 2 {
            warn!("sample_group {} has <2 timepoints; skipping", group);
            continue;
        }
        entries.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        for pair in entries.windows(2) {
            let (t0, idx0) = pair[0];
            let (t1, idx1) = pair[1];

            let d_ssm = survival.ssm[idx1] - survival.ssm[idx0];
            let d_afp = autophagy.afp[idx1] - autophagy.afp[idx0];
            let d_lds = lysosome.lds[idx1] - lysosome.lds[idx0];
            let d_ldi = damage.ldi[idx1] - damage.ldi[idx0];
            let d_erdi =
                cross.energy_recycling_dependency[idx1] - cross.energy_recycling_dependency[idx0];
            let d_prolif = survival.prolif[idx1] - survival.prolif[idx0];

            let asi = 0.4 * d_ssm + 0.3 * d_afp + 0.3 * d_lds;

            let response = if d_prolif < 0.0 && d_ssm < 0.0 {
                TherapyResponseClass::CytotoxicResponse
            } else if d_prolif < 0.0 && d_ssm > 0.0 {
                TherapyResponseClass::AdaptiveSurvival
            } else if d_lds > 0.0 && d_afp > 0.0 && d_prolif <= 0.0 {
                TherapyResponseClass::LysoEscape
            } else if d_ldi > 0.0 && d_ssm <= 0.0 {
                TherapyResponseClass::DamageAccumulation
            } else {
                TherapyResponseClass::NoResponse
            };

            metrics.delta_ssm.push(d_ssm);
            metrics.delta_afp.push(d_afp);
            metrics.delta_lds.push(d_lds);
            metrics.delta_ldi.push(d_ldi);
            metrics.delta_erdi.push(d_erdi);
            metrics.delta_prolif.push(d_prolif);
            metrics.asi.push(asi);
            metrics.response_class.push(response);
            metrics.from_obs.push(idx0);
            metrics.to_obs.push(idx1);
            metrics.sample_group.push(group.clone());
            metrics.timepoint0.push(t0);
            metrics.timepoint1.push(t1);
        }
    }

    ctx.therapy_delta = Some(metrics);
    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage11_therapy_delta.rs"]
mod tests;
