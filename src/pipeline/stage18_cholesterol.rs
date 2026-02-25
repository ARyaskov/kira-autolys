use crate::genesets::v1;
use crate::model::ctx::{CholesterolTraffickingMetrics, Ctx};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct CholesterolAcc {
    export_sum: f32,
    export_cnt: u32,
    import_sum: f32,
    import_cnt: u32,
    efflux_sum: f32,
    efflux_cnt: u32,
}

pub fn run_stage18(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut metrics = CholesterolTraffickingMetrics {
        export_capacity: vec![0.0; n_obs],
        import_pressure: vec![0.0; n_obs],
        efflux_capacity: vec![0.0; n_obs],
        cholesterol_trap_index: vec![0.0; n_obs],
    };

    for obs_idx in 0..n_obs {
        let mut acc = CholesterolAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask & bits.export != 0 {
                acc.export_sum += value;
                acc.export_cnt += 1;
            }
            if mask & bits.import != 0 {
                acc.import_sum += value;
                acc.import_cnt += 1;
            }
            if mask & bits.efflux != 0 {
                acc.efflux_sum += value;
                acc.efflux_cnt += 1;
            }
        });

        let export = mean(acc.export_sum, acc.export_cnt);
        let import = mean(acc.import_sum, acc.import_cnt);
        let efflux = mean(acc.efflux_sum, acc.efflux_cnt);
        let trap = import / (export + efflux + EPS);

        metrics.export_capacity[obs_idx] = export;
        metrics.import_pressure[obs_idx] = import;
        metrics.efflux_capacity[obs_idx] = efflux;
        metrics.cholesterol_trap_index[obs_idx] = trap;
    }

    ctx.cholesterol = Some(metrics);
    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    export: u64,
    import: u64,
    efflux: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut export = None;
    let mut import = None;
    let mut efflux = None;

    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "CHOLESTEROL_EXPORT" => export = Some(bit),
            "CHOLESTEROL_IMPORT" => import = Some(bit),
            "CHOLESTEROL_EFFLUX" => efflux = Some(bit),
            _ => {}
        }
    }

    let export = export
        .ok_or_else(|| StageError::Validation("missing CHOLESTEROL_EXPORT panel".to_string()))?;
    let import = import
        .ok_or_else(|| StageError::Validation("missing CHOLESTEROL_IMPORT panel".to_string()))?;
    let efflux = efflux
        .ok_or_else(|| StageError::Validation("missing CHOLESTEROL_EFFLUX panel".to_string()))?;

    Ok(PanelBits {
        export,
        import,
        efflux,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage18_cholesterol.rs"]
mod tests;
