use crate::genesets::v1;
use crate::model::ctx::{Ctx, LysosomeMetrics, NormalizedExpr};
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
struct LysosomeAcc {
    vatp_sum: f32,
    vatp_cnt: u32,
    prot_sum: f32,
    prot_cnt: u32,
    mem_sum: f32,
    mem_cnt: u32,
    global_sum: f32,
    global_cnt: u32,
}

pub fn run_stage3(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = LysosomeMetrics {
        lds: vec![0.0; n_obs],
        vatp: vec![0.0; n_obs],
        prot: vec![0.0; n_obs],
        mem: vec![0.0; n_obs],
        global_load: vec![0.0; n_obs],
    };

    if should_parallel(n_obs) {
        stage3_compute_parallel(
            normalized.as_ref(),
            panel_mask,
            bits.vatp,
            bits.prot,
            bits.mem,
            &mut metrics,
        )?;
    } else {
        stage3_compute_block(
            normalized.as_ref(),
            panel_mask,
            bits.vatp,
            bits.prot,
            bits.mem,
            &mut metrics,
            0..n_obs,
        )?;
    }

    ctx.lysosome = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    vatp: u64,
    prot: u64,
    mem: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut vatp = None;
    let mut prot = None;
    let mut mem = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LYSO_VATP" => vatp = Some(bit),
            "LYSO_PROT" => prot = Some(bit),
            "LYSO_MEM" => mem = Some(bit),
            _ => {}
        }
    }

    let vatp = vatp.ok_or_else(|| StageError::Validation("missing LYSO_VATP panel".to_string()))?;
    let prot = prot.ok_or_else(|| StageError::Validation("missing LYSO_PROT panel".to_string()))?;
    let mem = mem.ok_or_else(|| StageError::Validation("missing LYSO_MEM panel".to_string()))?;

    Ok(PanelBits { vatp, prot, mem })
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

fn stage3_compute_block(
    normalized: &dyn NormalizedExpr,
    gene_panel_mask: &[u64],
    vatp_bit: u64,
    prot_bit: u64,
    mem_bit: u64,
    out: &mut LysosomeMetrics,
    range: std::ops::Range<usize>,
) -> Result<(), StageError> {
    for obs_idx in range {
        let mut acc = LysosomeAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            acc.global_sum += value;
            acc.global_cnt += 1;

            let mask = gene_panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & vatp_bit != 0 {
                acc.vatp_sum += value;
                acc.vatp_cnt += 1;
            }
            if mask & prot_bit != 0 {
                acc.prot_sum += value;
                acc.prot_cnt += 1;
            }
            if mask & mem_bit != 0 {
                acc.mem_sum += value;
                acc.mem_cnt += 1;
            }
        });

        let vatp = mean(acc.vatp_sum, acc.vatp_cnt);
        let prot = mean(acc.prot_sum, acc.prot_cnt);
        let mem = mean(acc.mem_sum, acc.mem_cnt);
        let global_load = mean(acc.global_sum, acc.global_cnt);

        let lds_raw = 0.5 * vatp + 0.3 * prot + 0.2 * mem;
        let lds = lds_raw / (global_load + 1e-6);

        out.lds[obs_idx] = lds;
        out.vatp[obs_idx] = vatp;
        out.prot[obs_idx] = prot;
        out.mem[obs_idx] = mem;
        out.global_load[obs_idx] = global_load;
    }
    Ok(())
}

fn stage3_compute_parallel(
    normalized: &dyn NormalizedExpr,
    gene_panel_mask: &[u64],
    vatp_bit: u64,
    prot_bit: u64,
    mem_bit: u64,
    out: &mut LysosomeMetrics,
) -> Result<(), StageError> {
    #[cfg(feature = "parallel")]
    {
        let n_obs = out.lds.len();
        let lds_ptr = SendPtr(out.lds.as_mut_ptr());
        let vatp_ptr = SendPtr(out.vatp.as_mut_ptr());
        let prot_ptr = SendPtr(out.prot.as_mut_ptr());
        let mem_ptr = SendPtr(out.mem.as_mut_ptr());
        let load_ptr = SendPtr(out.global_load.as_mut_ptr());

        (0..n_obs)
            .step_by(BLOCK_SIZE)
            .collect::<Vec<_>>()
            .into_par_iter()
            .for_each({
                let lds_ptr = lds_ptr;
                let vatp_ptr = vatp_ptr;
                let prot_ptr = prot_ptr;
                let mem_ptr = mem_ptr;
                let load_ptr = load_ptr;
                move |start| {
                    let end = (start + BLOCK_SIZE).min(n_obs);
                    let len = end - start;
                    let lds =
                        unsafe { std::slice::from_raw_parts_mut(lds_ptr.ptr().add(start), len) };
                    let vatp =
                        unsafe { std::slice::from_raw_parts_mut(vatp_ptr.ptr().add(start), len) };
                    let prot =
                        unsafe { std::slice::from_raw_parts_mut(prot_ptr.ptr().add(start), len) };
                    let mem =
                        unsafe { std::slice::from_raw_parts_mut(mem_ptr.ptr().add(start), len) };
                    let load =
                        unsafe { std::slice::from_raw_parts_mut(load_ptr.ptr().add(start), len) };

                    for (offset, obs_idx) in (start..end).enumerate() {
                        let mut acc = LysosomeAcc::default();
                        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
                            acc.global_sum += value;
                            acc.global_cnt += 1;

                            let mask = gene_panel_mask[gene_idx];
                            if mask == 0 {
                                return;
                            }
                            if mask & vatp_bit != 0 {
                                acc.vatp_sum += value;
                                acc.vatp_cnt += 1;
                            }
                            if mask & prot_bit != 0 {
                                acc.prot_sum += value;
                                acc.prot_cnt += 1;
                            }
                            if mask & mem_bit != 0 {
                                acc.mem_sum += value;
                                acc.mem_cnt += 1;
                            }
                        });

                        let vatp_v = mean(acc.vatp_sum, acc.vatp_cnt);
                        let prot_v = mean(acc.prot_sum, acc.prot_cnt);
                        let mem_v = mean(acc.mem_sum, acc.mem_cnt);
                        let load_v = mean(acc.global_sum, acc.global_cnt);

                        let lds_raw = 0.5 * vatp_v + 0.3 * prot_v + 0.2 * mem_v;
                        let lds_v = lds_raw / (load_v + 1e-6);

                        lds[offset] = lds_v;
                        vatp[offset] = vatp_v;
                        prot[offset] = prot_v;
                        mem[offset] = mem_v;
                        load[offset] = load_v;
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
            vatp_bit,
            prot_bit,
            mem_bit,
            out,
        );
        Ok(())
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage3_lysosome.rs"]
mod tests;
