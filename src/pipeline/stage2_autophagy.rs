use crate::genesets::v1;
use crate::model::ctx::{AutophagyMetrics, Ctx, NormalizedExpr};
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
struct AutophagyAcc {
    init_sum: f32,
    init_cnt: u32,
    elong_sum: f32,
    elong_cnt: u32,
    late_sum: f32,
    late_cnt: u32,
    cargo_sum: f32,
    cargo_cnt: u32,
}

pub fn run_stage2(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = AutophagyMetrics {
        afp: vec![0.0; n_obs],
        initiation: vec![0.0; n_obs],
        elongation: vec![0.0; n_obs],
        degradation: vec![0.0; n_obs],
        cargo: vec![0.0; n_obs],
        stall: vec![0.0; n_obs],
    };

    if should_parallel(n_obs) {
        stage2_compute_parallel(
            normalized.as_ref(),
            panel_mask,
            bits.init,
            bits.elong,
            bits.late,
            bits.cargo,
            &mut metrics,
        )?;
    } else {
        stage2_compute_block(
            normalized.as_ref(),
            panel_mask,
            bits.init,
            bits.elong,
            bits.late,
            bits.cargo,
            &mut metrics,
            0..n_obs,
        )?;
    }

    ctx.autophagy = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    init: u64,
    elong: u64,
    late: u64,
    cargo: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut init = None;
    let mut elong = None;
    let mut late = None;
    let mut cargo = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "AUTOPHAGY_INIT" => init = Some(bit),
            "AUTOPHAGY_ELONG" => elong = Some(bit),
            "AUTOPHAGY_LATE" => late = Some(bit),
            "AUTOPHAGY_CARGO" => cargo = Some(bit),
            _ => {}
        }
    }

    let init =
        init.ok_or_else(|| StageError::Validation("missing AUTOPHAGY_INIT panel".to_string()))?;
    let elong =
        elong.ok_or_else(|| StageError::Validation("missing AUTOPHAGY_ELONG panel".to_string()))?;
    let late =
        late.ok_or_else(|| StageError::Validation("missing AUTOPHAGY_LATE panel".to_string()))?;
    let cargo =
        cargo.ok_or_else(|| StageError::Validation("missing AUTOPHAGY_CARGO panel".to_string()))?;

    Ok(PanelBits {
        init,
        elong,
        late,
        cargo,
    })
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

fn stage2_compute_block(
    normalized: &dyn NormalizedExpr,
    gene_panel_mask: &[u64],
    init_bit: u64,
    elong_bit: u64,
    late_bit: u64,
    cargo_bit: u64,
    out: &mut AutophagyMetrics,
    range: std::ops::Range<usize>,
) -> Result<(), StageError> {
    for obs_idx in range {
        let mut acc = AutophagyAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = gene_panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & init_bit != 0 {
                acc.init_sum += value;
                acc.init_cnt += 1;
            }
            if mask & elong_bit != 0 {
                acc.elong_sum += value;
                acc.elong_cnt += 1;
            }
            if mask & late_bit != 0 {
                acc.late_sum += value;
                acc.late_cnt += 1;
            }
            if mask & cargo_bit != 0 {
                acc.cargo_sum += value;
                acc.cargo_cnt += 1;
            }
        });

        let init = mean(acc.init_sum, acc.init_cnt);
        let elong = mean(acc.elong_sum, acc.elong_cnt);
        let late = mean(acc.late_sum, acc.late_cnt);
        let cargo = mean(acc.cargo_sum, acc.cargo_cnt);

        let axis_early = 0.5 * (init + elong);
        let axis_late = late;
        let stall = (axis_early - axis_late).max(0.0);
        let afp_raw = 0.5 * axis_early + 0.3 * axis_late + 0.2 * cargo;
        let afp = afp_raw - 0.25 * stall;

        out.afp[obs_idx] = afp;
        out.initiation[obs_idx] = init;
        out.elongation[obs_idx] = elong;
        out.degradation[obs_idx] = late;
        out.cargo[obs_idx] = cargo;
        out.stall[obs_idx] = stall;
    }
    Ok(())
}

fn stage2_compute_parallel(
    normalized: &dyn NormalizedExpr,
    gene_panel_mask: &[u64],
    init_bit: u64,
    elong_bit: u64,
    late_bit: u64,
    cargo_bit: u64,
    out: &mut AutophagyMetrics,
) -> Result<(), StageError> {
    #[cfg(feature = "parallel")]
    {
        let n_obs = out.afp.len();
        let out_ptr_afp = SendPtr(out.afp.as_mut_ptr());
        let out_ptr_init = SendPtr(out.initiation.as_mut_ptr());
        let out_ptr_elong = SendPtr(out.elongation.as_mut_ptr());
        let out_ptr_late = SendPtr(out.degradation.as_mut_ptr());
        let out_ptr_cargo = SendPtr(out.cargo.as_mut_ptr());
        let out_ptr_stall = SendPtr(out.stall.as_mut_ptr());

        (0..n_obs)
            .step_by(BLOCK_SIZE)
            .collect::<Vec<_>>()
            .into_par_iter()
            .for_each({
                let out_ptr_afp = out_ptr_afp;
                let out_ptr_init = out_ptr_init;
                let out_ptr_elong = out_ptr_elong;
                let out_ptr_late = out_ptr_late;
                let out_ptr_cargo = out_ptr_cargo;
                let out_ptr_stall = out_ptr_stall;
                move |start| {
                    let end = (start + BLOCK_SIZE).min(n_obs);
                    // Safety: ranges are disjoint per block.
                    let afp = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr_afp.ptr().add(start), end - start)
                    };
                    let init = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr_init.ptr().add(start), end - start)
                    };
                    let elong = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr_elong.ptr().add(start), end - start)
                    };
                    let late = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr_late.ptr().add(start), end - start)
                    };
                    let cargo = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr_cargo.ptr().add(start), end - start)
                    };
                    let stall = unsafe {
                        std::slice::from_raw_parts_mut(out_ptr_stall.ptr().add(start), end - start)
                    };

                    for (offset, obs_idx) in (start..end).enumerate() {
                        let mut acc = AutophagyAcc::default();
                        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
                            let mask = gene_panel_mask[gene_idx];
                            if mask == 0 {
                                return;
                            }
                            if mask & init_bit != 0 {
                                acc.init_sum += value;
                                acc.init_cnt += 1;
                            }
                            if mask & elong_bit != 0 {
                                acc.elong_sum += value;
                                acc.elong_cnt += 1;
                            }
                            if mask & late_bit != 0 {
                                acc.late_sum += value;
                                acc.late_cnt += 1;
                            }
                            if mask & cargo_bit != 0 {
                                acc.cargo_sum += value;
                                acc.cargo_cnt += 1;
                            }
                        });

                        let init_v = mean(acc.init_sum, acc.init_cnt);
                        let elong_v = mean(acc.elong_sum, acc.elong_cnt);
                        let late_v = mean(acc.late_sum, acc.late_cnt);
                        let cargo_v = mean(acc.cargo_sum, acc.cargo_cnt);

                        let axis_early = 0.5 * (init_v + elong_v);
                        let axis_late = late_v;
                        let stall_v = (axis_early - axis_late).max(0.0);
                        let afp_raw = 0.5 * axis_early + 0.3 * axis_late + 0.2 * cargo_v;
                        let afp_v = afp_raw - 0.25 * stall_v;

                        afp[offset] = afp_v;
                        init[offset] = init_v;
                        elong[offset] = elong_v;
                        late[offset] = late_v;
                        cargo[offset] = cargo_v;
                        stall[offset] = stall_v;
                    }
                }
            });
        Ok(())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = (
            normalized,
            gene_panel_mask,
            init_bit,
            elong_bit,
            late_bit,
            cargo_bit,
            out,
        );
        Ok(())
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage2_autophagy.rs"]
mod tests;
