use crate::genesets::v1;
use crate::model::ctx::{CouplingMetrics, Ctx};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct CouplingAcc {
    ampk_sum: f32,
    ampk_cnt: u32,
    rag_sum: f32,
    rag_cnt: u32,
}

pub fn run_stage9(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let regulatory = ctx
        .regulatory
        .as_ref()
        .ok_or_else(|| StageError::Validation("regulatory metrics missing".to_string()))?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;

    let bits = panel_bits()?;

    let mut metrics = CouplingMetrics {
        ampk: vec![0.0; n_obs],
        ragulator: vec![0.0; n_obs],
        ampk_mtor_ratio: vec![0.0; n_obs],
        tfeb_lyso_alignment: vec![0.0; n_obs],
        rag_mtor_mismatch: vec![0.0; n_obs],
        locked_survival_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = CouplingAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.ampk != 0 {
                acc.ampk_sum += value;
                acc.ampk_cnt += 1;
            }
            if mask & bits.ragulator != 0 {
                acc.rag_sum += value;
                acc.rag_cnt += 1;
            }
        });

        let ampk = mean(acc.ampk_sum, acc.ampk_cnt);
        let rag = mean(acc.rag_sum, acc.rag_cnt);
        let mtor = regulatory.mtor[obs_idx];
        let tfeb = regulatory.tfeb[obs_idx];
        let lds = lysosome.lds[obs_idx];

        let ampk_mtor_ratio = ampk / (mtor + 1e-6);
        let tfeb_lyso_alignment = tfeb * lds;
        let rag_mtor_mismatch = rag - mtor;
        let locked_survival_index =
            0.35 * ampk_mtor_ratio + 0.35 * tfeb_lyso_alignment + 0.30 * rag_mtor_mismatch.max(0.0);

        metrics.ampk[obs_idx] = ampk;
        metrics.ragulator[obs_idx] = rag;
        metrics.ampk_mtor_ratio[obs_idx] = ampk_mtor_ratio;
        metrics.tfeb_lyso_alignment[obs_idx] = tfeb_lyso_alignment;
        metrics.rag_mtor_mismatch[obs_idx] = rag_mtor_mismatch;
        metrics.locked_survival_index[obs_idx] = locked_survival_index;
    }

    ctx.coupling = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    ampk: u64,
    ragulator: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut ampk = None;
    let mut ragulator = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "REG_AMPK" => ampk = Some(bit),
            "REG_RAGULATOR" => ragulator = Some(bit),
            _ => {}
        }
    }

    let ampk = ampk.ok_or_else(|| StageError::Validation("missing REG_AMPK panel".to_string()))?;
    let ragulator = ragulator
        .ok_or_else(|| StageError::Validation("missing REG_RAGULATOR panel".to_string()))?;

    Ok(PanelBits { ampk, ragulator })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage9_coupling.rs"]
mod tests;
