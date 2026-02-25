use crate::genesets::v1;
use crate::math::stats::RunningStats;
use crate::model::ctx::{Ctx, MembraneRepairMetrics};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct RepairAcc {
    escrt_sum: f32,
    escrt_cnt: u32,
    ca_sum: f32,
    ca_cnt: u32,
    lpe_sum: f32,
    lpe_cnt: u32,
    msr_sum: f32,
    msr_cnt: u32,
}

pub fn run_stage26(ctx: &mut Ctx) -> Result<(), StageError> {
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

    if ctx.autophagy.is_none() {
        return Err(StageError::Validation(
            "autophagy metrics missing; run stage2 first".to_string(),
        ));
    }
    if ctx.lysosome.is_none() {
        return Err(StageError::Validation(
            "lysosome metrics missing; run stage3 first".to_string(),
        ));
    }
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;

    let bits = panel_bits()?;

    let mut rma = vec![0.0; n_obs];
    let mut lpe = vec![0.0; n_obs];
    let mut msr = vec![0.0; n_obs];
    let mut erc = vec![0.0; n_obs];
    let mut rdr = vec![0.0; n_obs];

    for obs_idx in 0..n_obs {
        let mut acc = RepairAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & bits.escrt != 0 {
                acc.escrt_sum += value;
                acc.escrt_cnt += 1;
            }
            if mask & bits.ca != 0 {
                acc.ca_sum += value;
                acc.ca_cnt += 1;
            }
            if mask & bits.lysophagy != 0 {
                acc.lpe_sum += value;
                acc.lpe_cnt += 1;
            }
            if mask & bits.stabilization != 0 {
                acc.msr_sum += value;
                acc.msr_cnt += 1;
            }
        });

        let rma_sum = acc.escrt_sum + acc.ca_sum;
        let rma_cnt = acc.escrt_cnt + acc.ca_cnt;
        let rma_mean = mean(rma_sum, rma_cnt);
        let lpe_mean = mean(acc.lpe_sum, acc.lpe_cnt);
        let msr_mean = mean(acc.msr_sum, acc.msr_cnt);
        let erc_val = rma_mean + lpe_mean + msr_mean;

        rma[obs_idx] = rma_mean;
        lpe[obs_idx] = lpe_mean;
        msr[obs_idx] = msr_mean;
        erc[obs_idx] = erc_val;
        rdr[obs_idx] = erc_val / (damage.ldi[obs_idx] + EPS);
    }

    let (mean_rdr, std_rdr) = mean_std(&rdr);
    let (mean_erc, std_erc) = mean_std(&erc);
    let (mean_ldi, std_ldi) = mean_std(&damage.ldi);

    let mut lmrci = vec![0.0; n_obs];
    for i in 0..n_obs {
        let z_rdr = zscore(rdr[i], mean_rdr, std_rdr);
        let z_erc = zscore(erc[i], mean_erc, std_erc);
        let z_ldi = zscore(damage.ldi[i], mean_ldi, std_ldi);
        lmrci[i] = z_rdr + z_erc - z_ldi;
    }

    ctx.membrane_repair = Some(MembraneRepairMetrics {
        rma,
        lpe,
        msr,
        erc,
        rdr,
        lmrci,
    });

    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

fn mean_std(values: &[f32]) -> (f32, f32) {
    let mut stats = RunningStats::new();
    for &v in values {
        stats.update(v);
    }
    (stats.mean(), stats.std())
}

fn zscore(value: f32, mean: f32, std: f32) -> f32 {
    (value - mean) / (std + EPS)
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    escrt: u64,
    ca: u64,
    lysophagy: u64,
    stabilization: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut escrt = None;
    let mut ca = None;
    let mut lysophagy = None;
    let mut stabilization = None;
    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LYSO_ESCRT_REPAIR" => escrt = Some(bit),
            "LYSO_CA_RECRUIT" => ca = Some(bit),
            "LYSO_LYSOPHAGY_INIT" => lysophagy = Some(bit),
            "LYSO_MEM_STABILIZATION" => stabilization = Some(bit),
            _ => {}
        }
    }
    let escrt = escrt
        .ok_or_else(|| StageError::Validation("missing LYSO_ESCRT_REPAIR panel".to_string()))?;
    let ca =
        ca.ok_or_else(|| StageError::Validation("missing LYSO_CA_RECRUIT panel".to_string()))?;
    let lysophagy = lysophagy
        .ok_or_else(|| StageError::Validation("missing LYSO_LYSOPHAGY_INIT panel".to_string()))?;
    let stabilization = stabilization.ok_or_else(|| {
        StageError::Validation("missing LYSO_MEM_STABILIZATION panel".to_string())
    })?;
    Ok(PanelBits {
        escrt,
        ca,
        lysophagy,
        stabilization,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage26_membrane_repair.rs"]
mod tests;
