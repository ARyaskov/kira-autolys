use std::fs::{self, File};
use std::path::Path;

use crate::config::thresholds::Thresholds;
use crate::model::ctx::{Ctx, LysosomeDependencyClass};
use crate::stage_error::StageError;

pub fn write_cohort_signatures(ctx: &Ctx, out_dir: &Path) -> Result<(), StageError> {
    let n_obs = ctx.obs_ids.len();
    if n_obs == 0 {
        return Err(StageError::Validation("no observations".to_string()));
    }

    let survival = ctx
        .survival
        .as_ref()
        .ok_or_else(|| StageError::Validation("survival metrics missing".to_string()))?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;
    let classes = ctx
        .lysosome_class
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome classes missing".to_string()))?;
    let cross = ctx
        .cross_organelle
        .as_ref()
        .ok_or_else(|| StageError::Validation("cross-organelle metrics missing".to_string()))?;

    if classes.class.len() != n_obs || survival.flag_ssm_high.len() != n_obs {
        return Err(StageError::Validation(
            "metrics length does not match obs_ids".to_string(),
        ));
    }

    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);

    let mut ssm_high = 0usize;
    let mut ldi_lds_high = 0usize;
    let mut adaptive_erdi_high = 0usize;

    for idx in 0..n_obs {
        if survival.flag_ssm_high[idx] {
            ssm_high += 1;
        }
        if damage.ldi[idx] >= thresholds.ldi_high && lysosome.lds[idx] >= thresholds.lds_high {
            ldi_lds_high += 1;
        }
        if classes.class[idx] == LysosomeDependencyClass::AutophagyAdaptive
            && cross.energy_recycling_dependency[idx] >= thresholds.erdi_high
        {
            adaptive_erdi_high += 1;
        }
    }

    let denom = n_obs as f32;
    let json = serde_json::json!({
        "SSM_high": ssm_high as f32 / denom,
        "LDI_LDS_high": ldi_lds_high as f32 / denom,
        "AutophagyAdaptive_ERDI_high": adaptive_erdi_high as f32 / denom,
    });

    fs::create_dir_all(out_dir)?;
    let path = out_dir.join("cohort_signatures.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &json)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;

    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/cohort_signatures.rs"]
mod tests;
