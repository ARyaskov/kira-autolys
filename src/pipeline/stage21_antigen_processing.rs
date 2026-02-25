use crate::genesets::v1;
use crate::model::ctx::{AntigenProcessingMetrics, Ctx};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct AntigenAcc {
    mhc_sum: f32,
    mhc_cnt: u32,
    protease_sum: f32,
    protease_cnt: u32,
    accessory_sum: f32,
    accessory_cnt: u32,
}

pub fn run_stage21(ctx: &mut Ctx) -> Result<(), StageError> {
    let n_obs = ctx.n_obs;
    if ctx.obs_ids.len() != n_obs {
        return Err(StageError::Validation(
            "obs_ids length does not match n_obs".to_string(),
        ));
    }

    let panel_mask = &ctx.gene_panel_mask;
    if panel_mask.is_empty() {
        return Err(StageError::Validation(
            "gene_panel_mask is empty; run stage0 first".to_string(),
        ));
    }

    let normalized = ctx
        .normalized
        .as_ref()
        .ok_or_else(|| StageError::Validation("normalized expression missing".to_string()))?;

    let bits = panel_bits()?;

    let mut metrics = AntigenProcessingMetrics {
        mhc_expression_load: vec![0.0; n_obs],
        protease_load: vec![0.0; n_obs],
        loading_accessory_load: vec![0.0; n_obs],
        antigen_processing_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = AntigenAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.mhc != 0 {
                acc.mhc_sum += value;
                acc.mhc_cnt += 1;
            }
            if mask & bits.protease != 0 {
                acc.protease_sum += value;
                acc.protease_cnt += 1;
            }
            if mask & bits.accessory != 0 {
                acc.accessory_sum += value;
                acc.accessory_cnt += 1;
            }
        });

        let mhc = mean(acc.mhc_sum, acc.mhc_cnt);
        let protease = mean(acc.protease_sum, acc.protease_cnt);
        let accessory = mean(acc.accessory_sum, acc.accessory_cnt);
        let index = mhc * protease * accessory;

        metrics.mhc_expression_load[obs_idx] = mhc;
        metrics.protease_load[obs_idx] = protease;
        metrics.loading_accessory_load[obs_idx] = accessory;
        metrics.antigen_processing_index[obs_idx] = index;
    }

    ctx.antigen_processing = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    mhc: u64,
    protease: u64,
    accessory: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut mhc = None;
    let mut protease = None;
    let mut accessory = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "MHC_CLASS_II_CORE" => mhc = Some(bit),
            "ANTIGEN_PROCESSING_PROTEASES" => protease = Some(bit),
            "MHC_LOADING_ACCESSORY" => accessory = Some(bit),
            _ => {}
        }
    }

    let mhc =
        mhc.ok_or_else(|| StageError::Validation("missing MHC_CLASS_II_CORE panel".to_string()))?;
    let protease = protease.ok_or_else(|| {
        StageError::Validation("missing ANTIGEN_PROCESSING_PROTEASES panel".to_string())
    })?;
    let accessory = accessory
        .ok_or_else(|| StageError::Validation("missing MHC_LOADING_ACCESSORY panel".to_string()))?;

    Ok(PanelBits {
        mhc,
        protease,
        accessory,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage21_antigen_processing.rs"]
mod tests;
