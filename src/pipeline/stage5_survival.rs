use crate::config::thresholds::Thresholds;
use crate::genesets::v1;
use crate::math::simd::{SurvivalFlags, compute_ssm_and_flags, compute_z_scores};
use crate::math::stats::RunningStats;
use crate::model::ctx::{Ctx, SurvivalMetrics, SurvivalStats};
use crate::stage_error::StageError;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(feature = "parallel")]
const BLOCK_SIZE: usize = 8192;

#[cfg(feature = "parallel")]
#[derive(Clone, Copy)]
struct SendPtr<T>(*mut T);

#[cfg(feature = "parallel")]
unsafe impl<T> Send for SendPtr<T> {}

#[cfg(feature = "parallel")]
unsafe impl<T> Sync for SendPtr<T> {}

#[cfg(feature = "parallel")]
impl<T> SendPtr<T> {
    #[inline]
    fn ptr(self) -> *mut T {
        self.0
    }
}
#[derive(Debug, Clone, Copy, Default)]
struct CellFateAcc {
    prolif_sum: f32,
    prolif_cnt: u32,
    apop_pro_sum: f32,
    apop_pro_cnt: u32,
    apop_anti_sum: f32,
    apop_anti_cnt: u32,
}

pub fn run_stage5(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let autophagy = ctx
        .autophagy
        .as_ref()
        .ok_or_else(|| StageError::Validation("autophagy metrics missing".to_string()))?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    if ctx.regulatory.is_none() {
        return Err(StageError::Validation(
            "regulatory metrics missing".to_string(),
        ));
    }

    let normalized = ctx
        .normalized
        .as_ref()
        .ok_or_else(|| StageError::Validation("normalized expression missing".to_string()))?;

    let bits = panel_bits()?;

    let mut prolif = Vec::with_capacity(n_obs);
    let mut apop_ready = Vec::with_capacity(n_obs);

    for obs_idx in 0..n_obs {
        let mut acc = CellFateAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.prolif != 0 {
                acc.prolif_sum += value;
                acc.prolif_cnt += 1;
            }
            if mask & bits.apop_pro != 0 {
                acc.apop_pro_sum += value;
                acc.apop_pro_cnt += 1;
            }
            if mask & bits.apop_anti != 0 {
                acc.apop_anti_sum += value;
                acc.apop_anti_cnt += 1;
            }
        });
        let prolif_mean = mean(acc.prolif_sum, acc.prolif_cnt);
        let apop_pro_mean = mean(acc.apop_pro_sum, acc.apop_pro_cnt);
        let apop_anti_mean = mean(acc.apop_anti_sum, acc.apop_anti_cnt);
        prolif.push(prolif_mean);
        apop_ready.push(apop_pro_mean - apop_anti_mean);
    }

    let mut stats_afp = RunningStats::new();
    let mut stats_lds = RunningStats::new();
    let mut stats_prolif = RunningStats::new();
    let mut stats_apop = RunningStats::new();

    for idx in 0..n_obs {
        stats_afp.update(autophagy.afp[idx]);
        stats_lds.update(lysosome.lds[idx]);
        stats_prolif.update(prolif[idx]);
        stats_apop.update(apop_ready[idx]);
    }

    let stats = SurvivalStats {
        mean_afp: stats_afp.mean(),
        std_afp: stats_afp.std(),
        mean_lds: stats_lds.mean(),
        std_lds: stats_lds.std(),
        mean_prolif: stats_prolif.mean(),
        std_prolif: stats_prolif.std(),
        mean_apop: stats_apop.mean(),
        std_apop: stats_apop.std(),
    };

    let mut metrics = SurvivalMetrics {
        prolif,
        apop_ready,
        z_afp: Vec::with_capacity(n_obs),
        z_lds: Vec::with_capacity(n_obs),
        z_prolif: Vec::with_capacity(n_obs),
        z_apop: Vec::with_capacity(n_obs),
        ssm: Vec::with_capacity(n_obs),
        flag_ssm_high: Vec::with_capacity(n_obs),
        flag_afp_high: Vec::with_capacity(n_obs),
        flag_lds_high: Vec::with_capacity(n_obs),
        flag_prolif_low: Vec::with_capacity(n_obs),
        flag_apop_low: Vec::with_capacity(n_obs),
    };

    metrics.z_afp.resize(n_obs, 0.0);
    metrics.z_lds.resize(n_obs, 0.0);
    metrics.z_prolif.resize(n_obs, 0.0);
    metrics.z_apop.resize(n_obs, 0.0);
    metrics.ssm.resize(n_obs, 0.0);
    metrics.flag_ssm_high.resize(n_obs, false);
    metrics.flag_afp_high.resize(n_obs, false);
    metrics.flag_lds_high.resize(n_obs, false);
    metrics.flag_prolif_low.resize(n_obs, false);
    metrics.flag_apop_low.resize(n_obs, false);

    let thresholds = ctx.thresholds.unwrap_or_else(Thresholds::default);
    let mut flags = SurvivalFlags {
        ssm_high: &mut metrics.flag_ssm_high,
        afp_high: &mut metrics.flag_afp_high,
        lds_high: &mut metrics.flag_lds_high,
        prolif_low: &mut metrics.flag_prolif_low,
        apop_low: &mut metrics.flag_apop_low,
    };

    if should_parallel(n_obs) {
        compute_z_scores_parallel(
            &autophagy.afp,
            stats.mean_afp,
            stats.std_afp,
            &mut metrics.z_afp,
        );
        compute_z_scores_parallel(
            &lysosome.lds,
            stats.mean_lds,
            stats.std_lds,
            &mut metrics.z_lds,
        );
        compute_z_scores_parallel(
            &metrics.prolif,
            stats.mean_prolif,
            stats.std_prolif,
            &mut metrics.z_prolif,
        );
        compute_z_scores_parallel(
            &metrics.apop_ready,
            stats.mean_apop,
            stats.std_apop,
            &mut metrics.z_apop,
        );
        compute_ssm_and_flags_parallel(
            &metrics.z_afp,
            &metrics.z_lds,
            &metrics.z_prolif,
            &metrics.z_apop,
            &thresholds,
            &mut metrics.ssm,
            &mut flags,
        );
    } else {
        let range = 0..n_obs;
        compute_z_scores_block(
            &autophagy.afp,
            stats.mean_afp,
            stats.std_afp,
            &mut metrics.z_afp,
            range.clone(),
        );
        compute_z_scores_block(
            &lysosome.lds,
            stats.mean_lds,
            stats.std_lds,
            &mut metrics.z_lds,
            range.clone(),
        );
        compute_z_scores_block(
            &metrics.prolif,
            stats.mean_prolif,
            stats.std_prolif,
            &mut metrics.z_prolif,
            range.clone(),
        );
        compute_z_scores_block(
            &metrics.apop_ready,
            stats.mean_apop,
            stats.std_apop,
            &mut metrics.z_apop,
            range.clone(),
        );
        compute_ssm_and_flags_block(
            &metrics.z_afp,
            &metrics.z_lds,
            &metrics.z_prolif,
            &metrics.z_apop,
            &thresholds,
            &mut metrics.ssm,
            &mut flags,
            range,
        );
    }

    ctx.survival = Some(metrics);
    ctx.survival_stats = Some(stats);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

fn should_parallel(n_obs: usize) -> bool {
    #[cfg(feature = "parallel")]
    {
        n_obs >= BLOCK_SIZE
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = n_obs;
        false
    }
}

fn compute_z_scores_block(
    values: &[f32],
    mean: f32,
    std: f32,
    out: &mut [f32],
    range: std::ops::Range<usize>,
) {
    compute_z_scores(&values[range.clone()], mean, std, &mut out[range]);
}

fn compute_ssm_and_flags_block(
    z_afp: &[f32],
    z_lds: &[f32],
    z_prolif: &[f32],
    z_apop: &[f32],
    thresholds: &Thresholds,
    out_ssm: &mut [f32],
    flags: &mut SurvivalFlags<'_>,
    range: std::ops::Range<usize>,
) {
    let mut sub_flags = SurvivalFlags {
        ssm_high: &mut flags.ssm_high[range.clone()],
        afp_high: &mut flags.afp_high[range.clone()],
        lds_high: &mut flags.lds_high[range.clone()],
        prolif_low: &mut flags.prolif_low[range.clone()],
        apop_low: &mut flags.apop_low[range.clone()],
    };

    compute_ssm_and_flags(
        &z_afp[range.clone()],
        &z_lds[range.clone()],
        &z_prolif[range.clone()],
        &z_apop[range.clone()],
        thresholds,
        &mut out_ssm[range],
        &mut sub_flags,
    );
}

fn compute_z_scores_parallel(values: &[f32], mean: f32, std: f32, out: &mut [f32]) {
    #[cfg(feature = "parallel")]
    {
        let n = values.len();
        let out_ptr = SendPtr(out.as_mut_ptr());
        (0..n)
            .step_by(BLOCK_SIZE)
            .collect::<Vec<_>>()
            .into_par_iter()
            .for_each({
                let out_ptr = out_ptr;
                move |start| {
                    let end = (start + BLOCK_SIZE).min(n);
                    // Safety: ranges are disjoint per block.
                    let out_slice = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr.ptr().add(start), end - start)
                    };
                    compute_z_scores(&values[start..end], mean, std, out_slice);
                }
            });
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = (values, mean, std, out);
    }
}

fn compute_ssm_and_flags_parallel(
    z_afp: &[f32],
    z_lds: &[f32],
    z_prolif: &[f32],
    z_apop: &[f32],
    thresholds: &Thresholds,
    out_ssm: &mut [f32],
    flags: &mut SurvivalFlags<'_>,
) {
    #[cfg(feature = "parallel")]
    {
        let n = out_ssm.len();
        let ssm_ptr = SendPtr(out_ssm.as_mut_ptr());
        let flag_ssm_ptr = SendPtr(flags.ssm_high.as_mut_ptr());
        let flag_afp_ptr = SendPtr(flags.afp_high.as_mut_ptr());
        let flag_lds_ptr = SendPtr(flags.lds_high.as_mut_ptr());
        let flag_pro_ptr = SendPtr(flags.prolif_low.as_mut_ptr());
        let flag_apop_ptr = SendPtr(flags.apop_low.as_mut_ptr());

        (0..n)
            .step_by(BLOCK_SIZE)
            .collect::<Vec<_>>()
            .into_par_iter()
            .for_each({
                let ssm_ptr = ssm_ptr;
                let flag_ssm_ptr = flag_ssm_ptr;
                let flag_afp_ptr = flag_afp_ptr;
                let flag_lds_ptr = flag_lds_ptr;
                let flag_pro_ptr = flag_pro_ptr;
                let flag_apop_ptr = flag_apop_ptr;
                move |start| {
                    let end = (start + BLOCK_SIZE).min(n);
                    let len = end - start;
                    // Safety: ranges are disjoint per block.
                    let out_slice =
                        unsafe { std::slice::from_raw_parts_mut(ssm_ptr.ptr().add(start), len) };
                    let mut block_flags = SurvivalFlags {
                        ssm_high: unsafe {
                            std::slice::from_raw_parts_mut(flag_ssm_ptr.ptr().add(start), len)
                        },
                        afp_high: unsafe {
                            std::slice::from_raw_parts_mut(flag_afp_ptr.ptr().add(start), len)
                        },
                        lds_high: unsafe {
                            std::slice::from_raw_parts_mut(flag_lds_ptr.ptr().add(start), len)
                        },
                        prolif_low: unsafe {
                            std::slice::from_raw_parts_mut(flag_pro_ptr.ptr().add(start), len)
                        },
                        apop_low: unsafe {
                            std::slice::from_raw_parts_mut(flag_apop_ptr.ptr().add(start), len)
                        },
                    };

                    compute_ssm_and_flags(
                        &z_afp[start..end],
                        &z_lds[start..end],
                        &z_prolif[start..end],
                        &z_apop[start..end],
                        thresholds,
                        out_slice,
                        &mut block_flags,
                    );
                }
            });
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = (z_afp, z_lds, z_prolif, z_apop, thresholds, out_ssm, flags);
    }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    prolif: u64,
    apop_pro: u64,
    apop_anti: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut prolif = None;
    let mut apop_pro = None;
    let mut apop_anti = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "PROLIF" => prolif = Some(bit),
            "APOP_PRO" => apop_pro = Some(bit),
            "APOP_ANTI" => apop_anti = Some(bit),
            _ => {}
        }
    }

    let prolif =
        prolif.ok_or_else(|| StageError::Validation("missing PROLIF panel".to_string()))?;
    let apop_pro =
        apop_pro.ok_or_else(|| StageError::Validation("missing APOP_PRO panel".to_string()))?;
    let apop_anti =
        apop_anti.ok_or_else(|| StageError::Validation("missing APOP_ANTI panel".to_string()))?;

    Ok(PanelBits {
        prolif,
        apop_pro,
        apop_anti,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage5_survival.rs"]
mod tests;
