use crate::genesets::v1;
use crate::model::ctx::{Ctx, ERLysosomeContactMetrics};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct ContactAcc {
    tether_sum: f32,
    tether_cnt: u32,
    calcium_sum: f32,
    calcium_cnt: u32,
    adaptor_sum: f32,
    adaptor_cnt: u32,
}

pub fn run_stage19(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = ERLysosomeContactMetrics {
        tethering_load: vec![0.0; n_obs],
        calcium_transfer_load: vec![0.0; n_obs],
        adaptor_stress: vec![0.0; n_obs],
        contact_stress_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = ContactAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.tethering != 0 {
                acc.tether_sum += value;
                acc.tether_cnt += 1;
            }
            if mask & bits.calcium != 0 {
                acc.calcium_sum += value;
                acc.calcium_cnt += 1;
            }
            if mask & bits.adaptor != 0 {
                acc.adaptor_sum += value;
                acc.adaptor_cnt += 1;
            }
        });

        let tether = mean(acc.tether_sum, acc.tether_cnt);
        let calcium = mean(acc.calcium_sum, acc.calcium_cnt);
        let adaptor = mean(acc.adaptor_sum, acc.adaptor_cnt);
        let index = tether * calcium * adaptor;

        metrics.tethering_load[obs_idx] = tether;
        metrics.calcium_transfer_load[obs_idx] = calcium;
        metrics.adaptor_stress[obs_idx] = adaptor;
        metrics.contact_stress_index[obs_idx] = index;
    }

    ctx.er_lysosome_contact = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    tethering: u64,
    calcium: u64,
    adaptor: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut tethering = None;
    let mut calcium = None;
    let mut adaptor = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "ER_LYSO_CONTACT_TETHERING" => tethering = Some(bit),
            "ER_LYSO_CALCIUM_AXIS" => calcium = Some(bit),
            "ER_LYSO_STRESS_ADAPTORS" => adaptor = Some(bit),
            _ => {}
        }
    }

    let tethering = tethering.ok_or_else(|| {
        StageError::Validation("missing ER_LYSO_CONTACT_TETHERING panel".to_string())
    })?;
    let calcium = calcium
        .ok_or_else(|| StageError::Validation("missing ER_LYSO_CALCIUM_AXIS panel".to_string()))?;
    let adaptor = adaptor.ok_or_else(|| {
        StageError::Validation("missing ER_LYSO_STRESS_ADAPTORS panel".to_string())
    })?;

    Ok(PanelBits {
        tethering,
        calcium,
        adaptor,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage19_er_lysosome.rs"]
mod tests;
