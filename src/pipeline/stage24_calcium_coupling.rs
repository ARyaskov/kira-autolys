use crate::genesets::v1;
use crate::math::stats::RunningStats;
use crate::model::ctx::{CalciumCouplingMetrics, Ctx};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct CalciumAcc {
    lcrc_sum: f32,
    lcrc_cnt: u32,
    mcuc_sum: f32,
    mcuc_cnt: u32,
    csap_sum: f32,
    csap_cnt: u32,
}

pub fn run_stage24(ctx: &mut Ctx) -> Result<(), StageError> {
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
    if ctx.cross_organelle.is_none() {
        return Err(StageError::Validation(
            "cross-organelle metrics missing; run stage10 first".to_string(),
        ));
    }
    let _ = ctx.mito.as_ref(); // optional; present only for interpretation downstream

    let bits = panel_bits()?;

    let mut lcrc = vec![0.0; n_obs];
    let mut mcuc = vec![0.0; n_obs];
    let mut csap = vec![0.0; n_obs];
    let mut ccb = vec![0.0; n_obs];

    for obs_idx in 0..n_obs {
        let mut acc = CalciumAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & bits.lcrc != 0 {
                acc.lcrc_sum += value;
                acc.lcrc_cnt += 1;
            }
            if mask & bits.mcuc != 0 {
                acc.mcuc_sum += value;
                acc.mcuc_cnt += 1;
            }
            if mask & bits.csap != 0 {
                acc.csap_sum += value;
                acc.csap_cnt += 1;
            }
        });

        let lcrc_mean = mean(acc.lcrc_sum, acc.lcrc_cnt);
        let mcuc_mean = mean(acc.mcuc_sum, acc.mcuc_cnt);
        let csap_mean = mean(acc.csap_sum, acc.csap_cnt);

        lcrc[obs_idx] = lcrc_mean;
        mcuc[obs_idx] = mcuc_mean;
        csap[obs_idx] = csap_mean;
        ccb[obs_idx] = lcrc_mean / (mcuc_mean + EPS);
    }

    let (mean_lcrc, std_lcrc) = mean_std(&lcrc);
    let (mean_mcuc, std_mcuc) = mean_std(&mcuc);
    let (mean_csap, std_csap) = mean_std(&csap);
    let (mean_ccb, std_ccb) = mean_std(&ccb);

    let mut lmcci = vec![0.0; n_obs];
    for i in 0..n_obs {
        let z_lcrc = zscore(lcrc[i], mean_lcrc, std_lcrc);
        let z_mcuc = zscore(mcuc[i], mean_mcuc, std_mcuc);
        let z_csap = zscore(csap[i], mean_csap, std_csap);
        let z_ccb = zscore(ccb[i], mean_ccb, std_ccb);
        lmcci[i] = z_lcrc + z_mcuc + z_csap - z_ccb.abs();
    }

    ctx.calcium_coupling = Some(CalciumCouplingMetrics {
        lcrc,
        mcuc,
        csap,
        ccb,
        lmcci,
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
    lcrc: u64,
    mcuc: u64,
    csap: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut lcrc = None;
    let mut mcuc = None;
    let mut csap = None;
    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LYSO_CA_RELEASE" => lcrc = Some(bit),
            "MITO_CA_UPTAKE" => mcuc = Some(bit),
            "CA_SIGNALING_ADAPTORS" => csap = Some(bit),
            "CA_METABOLIC_ACTIVATION" => {
                if let Some(existing) = csap {
                    csap = Some(existing | bit);
                } else {
                    csap = Some(bit);
                }
            }
            _ => {}
        }
    }
    let lcrc =
        lcrc.ok_or_else(|| StageError::Validation("missing LYSO_CA_RELEASE panel".to_string()))?;
    let mcuc =
        mcuc.ok_or_else(|| StageError::Validation("missing MITO_CA_UPTAKE panel".to_string()))?;
    let csap = csap
        .ok_or_else(|| StageError::Validation("missing CA_SIGNALING_ADAPTORS panel".to_string()))?;
    Ok(PanelBits { lcrc, mcuc, csap })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage24_calcium_coupling.rs"]
mod tests;
