use crate::genesets::v1;
use crate::model::ctx::{AutophagySelectivityMetrics, Ctx};
use crate::stage_error::StageError;

#[derive(Debug, Clone, Copy, Default)]
struct SelectiveAutophagyAcc {
    mito_sum: f32,
    mito_cnt: u32,
    aggre_sum: f32,
    aggre_cnt: u32,
    er_sum: f32,
    er_cnt: u32,
    ferr_sum: f32,
    ferr_cnt: u32,
    lipo_sum: f32,
    lipo_cnt: u32,
}

pub fn run_stage8(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = AutophagySelectivityMetrics {
        mitophagy: vec![0.0; n_obs],
        aggrephagy: vec![0.0; n_obs],
        erphagy: vec![0.0; n_obs],
        ferritinophagy: vec![0.0; n_obs],
        lipophagy: vec![0.0; n_obs],
        mito_frac: vec![0.0; n_obs],
        aggre_frac: vec![0.0; n_obs],
        er_frac: vec![0.0; n_obs],
        ferr_frac: vec![0.0; n_obs],
        lipo_frac: vec![0.0; n_obs],
        entropy: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = SelectiveAutophagyAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.mito != 0 {
                acc.mito_sum += value;
                acc.mito_cnt += 1;
            }
            if mask & bits.aggre != 0 {
                acc.aggre_sum += value;
                acc.aggre_cnt += 1;
            }
            if mask & bits.er != 0 {
                acc.er_sum += value;
                acc.er_cnt += 1;
            }
            if mask & bits.ferr != 0 {
                acc.ferr_sum += value;
                acc.ferr_cnt += 1;
            }
            if mask & bits.lipo != 0 {
                acc.lipo_sum += value;
                acc.lipo_cnt += 1;
            }
        });

        let mito = mean(acc.mito_sum, acc.mito_cnt);
        let aggre = mean(acc.aggre_sum, acc.aggre_cnt);
        let er = mean(acc.er_sum, acc.er_cnt);
        let ferr = mean(acc.ferr_sum, acc.ferr_cnt);
        let lipo = mean(acc.lipo_sum, acc.lipo_cnt);

        let total = mito + aggre + er + ferr + lipo + 1e-6;

        let mito_frac = mito / total;
        let aggre_frac = aggre / total;
        let er_frac = er / total;
        let ferr_frac = ferr / total;
        let lipo_frac = lipo / total;

        let entropy = -(mito_frac * (mito_frac + 1e-6).ln()
            + aggre_frac * (aggre_frac + 1e-6).ln()
            + er_frac * (er_frac + 1e-6).ln()
            + ferr_frac * (ferr_frac + 1e-6).ln()
            + lipo_frac * (lipo_frac + 1e-6).ln());

        metrics.mitophagy[obs_idx] = mito;
        metrics.aggrephagy[obs_idx] = aggre;
        metrics.erphagy[obs_idx] = er;
        metrics.ferritinophagy[obs_idx] = ferr;
        metrics.lipophagy[obs_idx] = lipo;
        metrics.mito_frac[obs_idx] = mito_frac;
        metrics.aggre_frac[obs_idx] = aggre_frac;
        metrics.er_frac[obs_idx] = er_frac;
        metrics.ferr_frac[obs_idx] = ferr_frac;
        metrics.lipo_frac[obs_idx] = lipo_frac;
        metrics.entropy[obs_idx] = entropy;
    }

    ctx.autophagy_selectivity = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    mito: u64,
    aggre: u64,
    er: u64,
    ferr: u64,
    lipo: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut mito = None;
    let mut aggre = None;
    let mut er = None;
    let mut ferr = None;
    let mut lipo = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "AUTO_MITOPHAGY" => mito = Some(bit),
            "AUTO_AGGREPHAGY" => aggre = Some(bit),
            "AUTO_ERPHAGY" => er = Some(bit),
            "AUTO_FERRITINOPHAGY" => ferr = Some(bit),
            "AUTO_LIPOPHAGY" => lipo = Some(bit),
            _ => {}
        }
    }

    let mito =
        mito.ok_or_else(|| StageError::Validation("missing AUTO_MITOPHAGY panel".to_string()))?;
    let aggre =
        aggre.ok_or_else(|| StageError::Validation("missing AUTO_AGGREPHAGY panel".to_string()))?;
    let er = er.ok_or_else(|| StageError::Validation("missing AUTO_ERPHAGY panel".to_string()))?;
    let ferr = ferr
        .ok_or_else(|| StageError::Validation("missing AUTO_FERRITINOPHAGY panel".to_string()))?;
    let lipo =
        lipo.ok_or_else(|| StageError::Validation("missing AUTO_LIPOPHAGY panel".to_string()))?;

    Ok(PanelBits {
        mito,
        aggre,
        er,
        ferr,
        lipo,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage8_autophagy_selectivity.rs"]
mod tests;
