use crate::genesets::v1;
use crate::model::ctx::{Ctx, FerroptosisMetrics};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct FerroptosisAcc {
    ferritin_sum: f32,
    ferritin_cnt: u32,
    defense_sum: f32,
    defense_cnt: u32,
    lipid_sum: f32,
    lipid_cnt: u32,
    tfrc_sum: f32,
    tfrc_cnt: u32,
    slc40a1_sum: f32,
    slc40a1_cnt: u32,
}

pub fn run_stage17(ctx: &mut Ctx) -> Result<(), StageError> {
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

    if ctx.gene_symbols.is_empty() {
        return Err(StageError::Validation(
            "gene_symbols missing; run stage0 first".to_string(),
        ));
    }

    let bits = panel_bits()?;

    let mut metrics = FerroptosisMetrics {
        ferritinophagy_load: vec![0.0; n_obs],
        iron_import_bias: vec![0.0; n_obs],
        ferroptosis_defense: vec![0.0; n_obs],
        lipid_peroxidation_context: vec![0.0; n_obs],
        ferroptotic_pressure_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = FerroptosisAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.ferritin != 0 {
                acc.ferritin_sum += value;
                acc.ferritin_cnt += 1;
            }
            if mask & bits.defense != 0 {
                acc.defense_sum += value;
                acc.defense_cnt += 1;
            }
            if mask & bits.lipid != 0 {
                acc.lipid_sum += value;
                acc.lipid_cnt += 1;
            }
            if mask & bits.iron != 0 {
                if let Some(symbol) = ctx.gene_symbols.get(gene_idx) {
                    if symbol == "TFRC" {
                        acc.tfrc_sum += value;
                        acc.tfrc_cnt += 1;
                    } else if symbol == "SLC40A1" {
                        acc.slc40a1_sum += value;
                        acc.slc40a1_cnt += 1;
                    }
                }
            }
        });

        let ferritin = mean(acc.ferritin_sum, acc.ferritin_cnt);
        let defense = mean(acc.defense_sum, acc.defense_cnt);
        let lipid = mean(acc.lipid_sum, acc.lipid_cnt);
        let tfrc = mean(acc.tfrc_sum, acc.tfrc_cnt);
        let slc40a1 = mean(acc.slc40a1_sum, acc.slc40a1_cnt);

        let iron_import_bias = tfrc / (slc40a1 + EPS);
        let pressure = (ferritin * iron_import_bias * lipid) / (defense + EPS);

        metrics.ferritinophagy_load[obs_idx] = ferritin;
        metrics.iron_import_bias[obs_idx] = iron_import_bias;
        metrics.ferroptosis_defense[obs_idx] = defense;
        metrics.lipid_peroxidation_context[obs_idx] = lipid;
        metrics.ferroptotic_pressure_index[obs_idx] = pressure;
    }

    ctx.ferroptosis = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    ferritin: u64,
    defense: u64,
    lipid: u64,
    iron: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut ferritin = None;
    let mut defense = None;
    let mut lipid = None;
    let mut iron = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "FERRITINOPHAGY" => ferritin = Some(bit),
            "FERROPTOSIS_DEFENSE" => defense = Some(bit),
            "LIPID_PEROXIDATION_CONTEXT" => lipid = Some(bit),
            "IRON_IMPORT_EXPORT" => iron = Some(bit),
            _ => {}
        }
    }

    let ferritin = ferritin
        .ok_or_else(|| StageError::Validation("missing FERRITINOPHAGY panel".to_string()))?;
    let defense = defense
        .ok_or_else(|| StageError::Validation("missing FERROPTOSIS_DEFENSE panel".to_string()))?;
    let lipid = lipid.ok_or_else(|| {
        StageError::Validation("missing LIPID_PEROXIDATION_CONTEXT panel".to_string())
    })?;
    let iron =
        iron.ok_or_else(|| StageError::Validation("missing IRON_IMPORT_EXPORT panel".to_string()))?;

    Ok(PanelBits {
        ferritin,
        defense,
        lipid,
        iron,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage17_ferroptosis.rs"]
mod tests;
