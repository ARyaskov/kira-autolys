use std::fs::{self, File};
use std::path::Path;

use crate::model::cohort_view::CohortView;
use crate::model::ctx::LysosomeDependencyClass;
use crate::stage_error::StageError;

pub fn write_signatures(view: &CohortView, out_dir: &Path) -> Result<(), StageError> {
    let n_obs = view.obs_ids.len();
    if n_obs == 0 {
        return Err(StageError::Validation("no observations".to_string()));
    }

    let mut ssm_high = 0usize;
    let mut ldi_lds_high = 0usize;
    let mut adaptive_erdi_high = 0usize;

    for idx in 0..n_obs {
        if view.ssm[idx] >= view.thresholds.ssm_high {
            ssm_high += 1;
        }
        if view.ldi[idx] >= view.thresholds.ldi_high && view.lds[idx] >= view.thresholds.lds_high {
            ldi_lds_high += 1;
        }
        if view.classes[idx] == LysosomeDependencyClass::AutophagyAdaptive
            && view.erdi[idx] >= view.thresholds.erdi_high
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
#[path = "../../../tests/src_inline/report/cohort/signatures.rs"]
mod tests;
