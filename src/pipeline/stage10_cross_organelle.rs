use crate::model::ctx::{CrossOrganelleMetrics, Ctx};
use crate::stage_error::StageError;

pub fn run_stage10(ctx: &mut Ctx) -> Result<(), StageError> {
    let n_obs = ctx.n_obs;
    if ctx.obs_ids.len() != n_obs {
        return Err(StageError::Validation(
            "obs_ids length does not match n_obs".to_string(),
        ));
    }

    let selectivity = ctx.autophagy_selectivity.as_ref().ok_or_else(|| {
        StageError::Validation("autophagy selectivity metrics missing".to_string())
    })?;
    let lysosome = ctx
        .lysosome
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome metrics missing".to_string()))?;
    let damage = ctx
        .lysosomal_damage
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosomal damage metrics missing".to_string()))?;
    let coupling = ctx
        .coupling
        .as_ref()
        .ok_or_else(|| StageError::Validation("coupling metrics missing".to_string()))?;

    let mito = ctx.mito.as_ref();

    let mut metrics = CrossOrganelleMetrics {
        mitophagy_reliance: vec![0.0; n_obs],
        mito_lyso_imbalance: vec![f32::NAN; n_obs],
        mito_consumption_risk: vec![f32::NAN; n_obs],
        energy_recycling_dependency: vec![0.0; n_obs],
    };

    for idx in 0..n_obs {
        let mitophagy = selectivity.mitophagy[idx];
        let lds = lysosome.lds[idx];
        let mitophagy_reliance = mitophagy * lds;

        let mut mito_lyso_imbalance = f32::NAN;
        let mut mito_consumption_risk = f32::NAN;
        if let Some(mito) = mito {
            let mito_stress = mito.mito_stress[idx];
            let mito_mass = mito.mito_mass_proxy[idx];
            mito_lyso_imbalance = mito_stress / (lds + 1e-6);
            mito_consumption_risk = mitophagy_reliance / (mito_mass + 1e-6);
        }

        let locked = coupling.locked_survival_index[idx].max(0.0);
        let ldi = damage.ldi[idx];
        let erdi = 0.4 * mitophagy_reliance + 0.3 * locked + 0.3 * ldi;

        metrics.mitophagy_reliance[idx] = mitophagy_reliance;
        metrics.mito_lyso_imbalance[idx] = mito_lyso_imbalance;
        metrics.mito_consumption_risk[idx] = mito_consumption_risk;
        metrics.energy_recycling_dependency[idx] = erdi;
    }

    ctx.cross_organelle = Some(metrics);
    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage10_cross_organelle.rs"]
mod tests;
