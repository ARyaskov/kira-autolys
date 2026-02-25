use crate::genesets::v1;
use crate::model::ctx::{Ctx, LipidBufferingMetrics};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct LipidAcc {
    storage_sum: f32,
    storage_cnt: u32,
    util_sum: f32,
    util_cnt: u32,
    context_sum: f32,
    context_cnt: u32,
}

pub fn run_stage22(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = LipidBufferingMetrics {
        storage_buffering_load: vec![0.0; n_obs],
        utilization_load: vec![0.0; n_obs],
        lipophagy_context: vec![0.0; n_obs],
        lipid_buffering_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = LipidAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.storage != 0 {
                acc.storage_sum += value;
                acc.storage_cnt += 1;
            }
            if mask & bits.utilization != 0 {
                acc.util_sum += value;
                acc.util_cnt += 1;
            }
            if mask & bits.context != 0 {
                acc.context_sum += value;
                acc.context_cnt += 1;
            }
        });

        let storage = mean(acc.storage_sum, acc.storage_cnt);
        let utilization = mean(acc.util_sum, acc.util_cnt);
        let context = mean(acc.context_sum, acc.context_cnt);
        let index = (storage * context) / (utilization + EPS);

        metrics.storage_buffering_load[obs_idx] = storage;
        metrics.utilization_load[obs_idx] = utilization;
        metrics.lipophagy_context[obs_idx] = context;
        metrics.lipid_buffering_index[obs_idx] = index;
    }

    ctx.lipid_buffering = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    storage: u64,
    utilization: u64,
    context: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut storage = None;
    let mut utilization = None;
    let mut context = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LIPID_STORAGE_BUFFERING" => storage = Some(bit),
            "LIPID_MOBILIZATION_UTILIZATION" => utilization = Some(bit),
            "LIPID_LIPOPHAGY_CONTEXT" => context = Some(bit),
            _ => {}
        }
    }

    let storage = storage.ok_or_else(|| {
        StageError::Validation("missing LIPID_STORAGE_BUFFERING panel".to_string())
    })?;
    let utilization = utilization.ok_or_else(|| {
        StageError::Validation("missing LIPID_MOBILIZATION_UTILIZATION panel".to_string())
    })?;
    let context = context.ok_or_else(|| {
        StageError::Validation("missing LIPID_LIPOPHAGY_CONTEXT panel".to_string())
    })?;

    Ok(PanelBits {
        storage,
        utilization,
        context,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage22_lipid_buffering.rs"]
mod tests;
