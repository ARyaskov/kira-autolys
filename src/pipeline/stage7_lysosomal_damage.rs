use crate::genesets::v1;
use crate::model::ctx::{Ctx, LysosomalDamageMetrics};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct LysosomalDamageAcc {
    lmp_sum: f32,
    lmp_cnt: u32,
    stress_sum: f32,
    stress_cnt: u32,
}

pub fn run_stage7(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;

    let bits = panel_bits()?;

    let mut metrics = LysosomalDamageMetrics {
        lmp: vec![0.0; n_obs],
        stress: vec![0.0; n_obs],
        cathepsin_membrane_imbalance: vec![0.0; n_obs],
        ldi: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = LysosomalDamageAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.lmp != 0 {
                acc.lmp_sum += value;
                acc.lmp_cnt += 1;
            }
            if mask & bits.stress != 0 {
                acc.stress_sum += value;
                acc.stress_cnt += 1;
            }
        });

        let lmp = mean(acc.lmp_sum, acc.lmp_cnt);
        let stress = mean(acc.stress_sum, acc.stress_cnt);
        let prot = lysosome.prot[obs_idx];
        let mem = lysosome.mem[obs_idx];
        let cathepsin_membrane_imbalance = prot / (mem + 1e-6);
        let ldi = 0.4 * lmp + 0.3 * stress + 0.3 * cathepsin_membrane_imbalance;

        metrics.lmp[obs_idx] = lmp;
        metrics.stress[obs_idx] = stress;
        metrics.cathepsin_membrane_imbalance[obs_idx] = cathepsin_membrane_imbalance;
        metrics.ldi[obs_idx] = ldi;
    }

    ctx.lysosomal_damage = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    lmp: u64,
    stress: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut lmp = None;
    let mut stress = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LYSO_LMP" => lmp = Some(bit),
            "LYSO_STRESS" => stress = Some(bit),
            _ => {}
        }
    }

    let lmp = lmp.ok_or_else(|| StageError::Validation("missing LYSO_LMP panel".to_string()))?;
    let stress =
        stress.ok_or_else(|| StageError::Validation("missing LYSO_STRESS panel".to_string()))?;

    Ok(PanelBits { lmp, stress })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage7_lysosomal_damage.rs"]
mod tests;
