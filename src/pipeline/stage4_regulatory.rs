use crate::genesets::v1;
use crate::model::ctx::{Ctx, RegulatoryMetrics};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct RegulatoryAcc {
    tfeb_sum: f32,
    tfeb_cnt: u32,
    mtor_sum: f32,
    mtor_cnt: u32,
}

pub fn run_stage4(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut tfeb = Vec::with_capacity(n_obs);
    let mut mtor = Vec::with_capacity(n_obs);

    for obs_idx in 0..n_obs {
        let mut acc = RegulatoryAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & bits.tfeb != 0 {
                acc.tfeb_sum += value;
                acc.tfeb_cnt += 1;
            }
            if mask & bits.mtor != 0 {
                acc.mtor_sum += value;
                acc.mtor_cnt += 1;
            }
        });
        tfeb.push(mean(acc.tfeb_sum, acc.tfeb_cnt));
        mtor.push(mean(acc.mtor_sum, acc.mtor_cnt));
    }

    let median_tfeb = median(&tfeb);
    let median_mtor = median(&mtor);

    let mut metrics = RegulatoryMetrics {
        tfeb: tfeb.clone(),
        mtor: mtor.clone(),
        tfeb_act: Vec::with_capacity(n_obs),
        mtor_supp: Vec::with_capacity(n_obs),
        tfeb_mtor_diff: Vec::with_capacity(n_obs),
        tfeb_mtor_ratio: Vec::with_capacity(n_obs),
        median_tfeb,
        median_mtor,
    };

    for idx in 0..n_obs {
        let tfeb_p = tfeb[idx];
        let mtor_p = mtor[idx];
        let tfeb_act = tfeb_p - median_tfeb;
        let mtor_supp = (median_mtor - mtor_p).max(0.0);
        let diff = tfeb_p - mtor_p;
        let ratio = (tfeb_act + 1e-6) / ((mtor_p - median_mtor).abs() + 1e-6);

        metrics.tfeb_act.push(tfeb_act);
        metrics.mtor_supp.push(mtor_supp);
        metrics.tfeb_mtor_diff.push(diff);
        metrics.tfeb_mtor_ratio.push(ratio);
    }

    ctx.regulatory = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

fn median(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = data.len() / 2;
    if data.len() % 2 == 1 {
        data[mid]
    } else {
        0.5 * (data[mid - 1] + data[mid])
    }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    tfeb: u64,
    mtor: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut tfeb = None;
    let mut mtor = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "REG_TFEB" => tfeb = Some(bit),
            "REG_MTOR" => mtor = Some(bit),
            _ => {}
        }
    }

    let tfeb = tfeb.ok_or_else(|| StageError::Validation("missing REG_TFEB panel".to_string()))?;
    let mtor = mtor.ok_or_else(|| StageError::Validation("missing REG_MTOR panel".to_string()))?;

    Ok(PanelBits { tfeb, mtor })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage4_regulatory.rs"]
mod tests;
