use crate::genesets::v1;
use crate::model::ctx::{Ctx, LysosomePositioningMetrics};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct PositioningAcc {
    perinuclear_sum: f32,
    perinuclear_cnt: u32,
    peripheral_sum: f32,
    peripheral_cnt: u32,
}

pub fn run_stage16(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = LysosomePositioningMetrics {
        perinuclear_mean: vec![0.0; n_obs],
        peripheral_mean: vec![0.0; n_obs],
        positioning_ratio: vec![0.0; n_obs],
        positioning_bias: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = PositioningAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.perinuclear != 0 {
                acc.perinuclear_sum += value;
                acc.perinuclear_cnt += 1;
            }
            if mask & bits.peripheral != 0 {
                acc.peripheral_sum += value;
                acc.peripheral_cnt += 1;
            }
        });

        let perinuclear_mean = mean(acc.perinuclear_sum, acc.perinuclear_cnt);
        let peripheral_mean = mean(acc.peripheral_sum, acc.peripheral_cnt);
        let ratio = (perinuclear_mean + EPS) / (peripheral_mean + EPS);
        let bias = ratio.log2();

        metrics.perinuclear_mean[obs_idx] = perinuclear_mean;
        metrics.peripheral_mean[obs_idx] = peripheral_mean;
        metrics.positioning_ratio[obs_idx] = ratio;
        metrics.positioning_bias[obs_idx] = bias;
    }

    ctx.lysosome_positioning = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    perinuclear: u64,
    peripheral: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut perinuclear = None;
    let mut peripheral = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "PERINUCLEAR_LYSOSOME" => perinuclear = Some(bit),
            "PERIPHERAL_LYSOSOME" => peripheral = Some(bit),
            _ => {}
        }
    }

    let perinuclear = perinuclear
        .ok_or_else(|| StageError::Validation("missing PERINUCLEAR_LYSOSOME panel".to_string()))?;
    let peripheral = peripheral
        .ok_or_else(|| StageError::Validation("missing PERIPHERAL_LYSOSOME panel".to_string()))?;

    Ok(PanelBits {
        perinuclear,
        peripheral,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage16_positioning.rs"]
mod tests;
