use crate::genesets::v1;
use crate::model::ctx::{Ctx, SecretoryLysosomeMetrics};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct SecretoryAcc {
    trafficking_sum: f32,
    trafficking_cnt: u32,
    fusion_sum: f32,
    fusion_cnt: u32,
    regulation_sum: f32,
    regulation_cnt: u32,
}

pub fn run_stage20(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = SecretoryLysosomeMetrics {
        trafficking_load: vec![0.0; n_obs],
        fusion_load: vec![0.0; n_obs],
        exocytosis_regulation: vec![0.0; n_obs],
        secretory_bias_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = SecretoryAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.trafficking != 0 {
                acc.trafficking_sum += value;
                acc.trafficking_cnt += 1;
            }
            if mask & bits.fusion != 0 {
                acc.fusion_sum += value;
                acc.fusion_cnt += 1;
            }
            if mask & bits.regulation != 0 {
                acc.regulation_sum += value;
                acc.regulation_cnt += 1;
            }
        });

        let trafficking = mean(acc.trafficking_sum, acc.trafficking_cnt);
        let fusion = mean(acc.fusion_sum, acc.fusion_cnt);
        let regulation = mean(acc.regulation_sum, acc.regulation_cnt);
        let index = trafficking * fusion * regulation;

        metrics.trafficking_load[obs_idx] = trafficking;
        metrics.fusion_load[obs_idx] = fusion;
        metrics.exocytosis_regulation[obs_idx] = regulation;
        metrics.secretory_bias_index[obs_idx] = index;
    }

    ctx.secretory_lysosome = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    trafficking: u64,
    fusion: u64,
    regulation: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut trafficking = None;
    let mut fusion = None;
    let mut regulation = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "SECRETORY_LYSOSOME_TRAFFICKING" => trafficking = Some(bit),
            "SECRETORY_FUSION_MACHINERY" => fusion = Some(bit),
            "LYSOSOME_EXOCYTOSIS_REGULATORS" => regulation = Some(bit),
            _ => {}
        }
    }

    let trafficking = trafficking.ok_or_else(|| {
        StageError::Validation("missing SECRETORY_LYSOSOME_TRAFFICKING panel".to_string())
    })?;
    let fusion = fusion.ok_or_else(|| {
        StageError::Validation("missing SECRETORY_FUSION_MACHINERY panel".to_string())
    })?;
    let regulation = regulation.ok_or_else(|| {
        StageError::Validation("missing LYSOSOME_EXOCYTOSIS_REGULATORS panel".to_string())
    })?;

    Ok(PanelBits {
        trafficking,
        fusion,
        regulation,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage20_secretory_lysosome.rs"]
mod tests;
