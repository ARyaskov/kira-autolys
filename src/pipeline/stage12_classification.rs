use crate::config::thresholds::Thresholds;
use crate::model::ctx::{Ctx, LysosomeClassMetrics, LysosomeDependencyClass};
use crate::stage_error::StageError;

pub fn run_stage12(ctx: &mut Ctx) -> Result<(), StageError> {
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
    let survival = ctx
        .survival
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival metrics missing".to_string()))?;
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;
    let coupling = ctx
        .coupling
        .as_ref()
        .ok_or_else(|| StageError::Validation("coupling metrics missing".to_string()))?;
    let cross = ctx
        .cross_organelle
        .as_ref()
        .ok_or_else(|| StageError::Validation("cross-organelle metrics missing".to_string()))?;

    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);

    let mut classes = Vec::with_capacity(n_obs);

    for idx in 0..n_obs {
        let class = if cross.energy_recycling_dependency[idx] >= thresholds.erdi_high {
            LysosomeDependencyClass::EnergyRecyclingDependent
        } else if lysosome.lds[idx] >= thresholds.lds_high && damage.ldi[idx] >= thresholds.ldi_high
        {
            LysosomeDependencyClass::LysosomeOverloaded
        } else if survival.flag_ssm_high[idx]
            && autophagy.afp[idx] >= thresholds.afp_high
            && damage.ldi[idx] < thresholds.ldi_high
        {
            LysosomeDependencyClass::AutophagyAdaptive
        } else if autophagy.stall[idx] >= thresholds.stall_high {
            LysosomeDependencyClass::StalledAutophagy
        } else if lysosome.lds[idx] >= thresholds.lds_high && !survival.flag_ssm_high[idx] {
            LysosomeDependencyClass::LysosomeDependent
        } else if lysosome.lds[idx] < thresholds.lds_low && autophagy.afp[idx] < thresholds.afp_low
        {
            LysosomeDependencyClass::NonLysosomal
        } else {
            LysosomeDependencyClass::Unclassified
        };

        let _ = coupling; // coupling is required for stage, but not directly used in rules now.
        classes.push(class);
    }

    ctx.lysosome_class = Some(LysosomeClassMetrics { class: classes });
    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage12_classification.rs"]
mod tests;
