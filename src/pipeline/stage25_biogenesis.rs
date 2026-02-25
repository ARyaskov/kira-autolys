use crate::genesets::v1;
use crate::math::stats::RunningStats;
use crate::model::ctx::{BiogenesisMetrics, Ctx};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct BiogenesisAcc {
    bd_sum: f32,
    bd_cnt: u32,
    my_sum: f32,
    my_cnt: u32,
}

pub fn run_stage25(ctx: &mut Ctx) -> Result<(), StageError> {
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

    if ctx.regulatory.is_none() {
        return Err(StageError::Validation(
            "regulatory metrics missing; run stage4 first".to_string(),
        ));
    }
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;

    let bits = panel_bits()?;

    let mut biogenesis_drive = vec![0.0; n_obs];
    let mut maturation_yield = vec![0.0; n_obs];

    for obs_idx in 0..n_obs {
        let mut acc = BiogenesisAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & bits.biogenesis != 0 {
                acc.bd_sum += value;
                acc.bd_cnt += 1;
            }
            if mask & bits.maturation != 0 {
                acc.my_sum += value;
                acc.my_cnt += 1;
            }
        });
        biogenesis_drive[obs_idx] = mean(acc.bd_sum, acc.bd_cnt);
        maturation_yield[obs_idx] = mean(acc.my_sum, acc.my_cnt);
    }

    let (mean_ldi, std_ldi) = mean_std(&damage.ldi);
    let (mean_lasi, std_lasi, has_lasi) = match ctx.acidification.as_ref() {
        Some(acid) => {
            let (m, s) = mean_std(&acid.lasi);
            (m, s, true)
        }
        None => (0.0, 1.0, false),
    };

    let mut functional_capacity = vec![0.0; n_obs];
    let mut bpr = vec![0.0; n_obs];

    for i in 0..n_obs {
        let z_ldi = zscore(damage.ldi[i], mean_ldi, std_ldi);
        let z_lasi = if has_lasi {
            let lasi = ctx.acidification.as_ref().map(|a| a.lasi[i]).unwrap_or(0.0);
            zscore(lasi, mean_lasi, std_lasi)
        } else {
            0.0
        };
        let flc = lysosome.lds[i] * (1.0 - z_ldi) * (1.0 + z_lasi);
        functional_capacity[i] = flc;
        bpr[i] = biogenesis_drive[i] / (flc + EPS);
    }

    let (mean_bd, std_bd) = mean_std(&biogenesis_drive);
    let (mean_my, std_my) = mean_std(&maturation_yield);
    let (mean_bpr, std_bpr) = mean_std(&bpr);

    let mut lbpi = vec![0.0; n_obs];
    for i in 0..n_obs {
        let z_bd = zscore(biogenesis_drive[i], mean_bd, std_bd);
        let z_my = zscore(maturation_yield[i], mean_my, std_my);
        let z_bpr = zscore(bpr[i], mean_bpr, std_bpr);
        lbpi[i] = z_bpr + z_bd - z_my;
    }

    ctx.biogenesis = Some(BiogenesisMetrics {
        biogenesis_drive,
        maturation_yield,
        functional_capacity,
        bpr,
        lbpi,
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
    biogenesis: u64,
    maturation: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut biogenesis = None;
    let mut maturation = None;
    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LYSO_BIOGENESIS_PROGRAM" => biogenesis = Some(bit),
            "LYSO_ASSEMBLY_MATURATION" => maturation = Some(bit),
            _ => {}
        }
    }
    let biogenesis = biogenesis.ok_or_else(|| {
        StageError::Validation("missing LYSO_BIOGENESIS_PROGRAM panel".to_string())
    })?;
    let maturation = maturation.ok_or_else(|| {
        StageError::Validation("missing LYSO_ASSEMBLY_MATURATION panel".to_string())
    })?;
    Ok(PanelBits {
        biogenesis,
        maturation,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage25_biogenesis.rs"]
mod tests;
