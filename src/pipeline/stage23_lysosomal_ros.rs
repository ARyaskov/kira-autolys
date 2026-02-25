use crate::genesets::v1;
use crate::model::ctx::{Ctx, LysosomalROSMetrics};
use crate::stage_error::StageError;

const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, Default)]
struct RosAcc {
    ros_sum: f32,
    ros_cnt: u32,
    iron_sum: f32,
    iron_cnt: u32,
    anti_sum: f32,
    anti_cnt: u32,
    damage_sum: f32,
    damage_cnt: u32,
}

pub fn run_stage23(ctx: &mut Ctx) -> Result<(), StageError> {
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

    let mut ros_generation_load = vec![0.0; n_obs];
    let mut iron_redox_context = vec![0.0; n_obs];
    let mut antioxidant_capacity = vec![0.0; n_obs];
    let mut damage_response_load = vec![0.0; n_obs];
    let mut lysosomal_ros_stress_index = vec![0.0; n_obs];

    for obs_idx in 0..n_obs {
        let mut acc = RosAcc::default();
        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            let mask = panel_mask[gene_idx];
            if mask == 0 {
                return;
            }
            if mask & bits.ros != 0 {
                acc.ros_sum += value;
                acc.ros_cnt += 1;
            }
            if mask & bits.iron != 0 {
                acc.iron_sum += value;
                acc.iron_cnt += 1;
            }
            if mask & bits.anti != 0 {
                acc.anti_sum += value;
                acc.anti_cnt += 1;
            }
            if mask & bits.damage != 0 {
                acc.damage_sum += value;
                acc.damage_cnt += 1;
            }
        });

        let ros_mean = mean(acc.ros_sum, acc.ros_cnt);
        let iron_mean = mean(acc.iron_sum, acc.iron_cnt);
        let anti_mean = mean(acc.anti_sum, acc.anti_cnt);
        let damage_mean = mean(acc.damage_sum, acc.damage_cnt);
        let stress_index = (ros_mean * iron_mean) / (anti_mean + EPS);

        ros_generation_load[obs_idx] = ros_mean;
        iron_redox_context[obs_idx] = iron_mean;
        antioxidant_capacity[obs_idx] = anti_mean;
        damage_response_load[obs_idx] = damage_mean;
        lysosomal_ros_stress_index[obs_idx] = stress_index;
    }

    ctx.lysosomal_ros = Some(LysosomalROSMetrics {
        ros_generation_load,
        iron_redox_context,
        antioxidant_capacity,
        damage_response_load,
        lysosomal_ros_stress_index,
    });

    Ok(())
}

fn mean(sum: f32, cnt: u32) -> f32 {
    if cnt == 0 { 0.0 } else { sum / cnt as f32 }
}

#[derive(Debug, Clone, Copy)]
struct PanelBits {
    ros: u64,
    iron: u64,
    anti: u64,
    damage: u64,
}

fn panel_bits() -> Result<PanelBits, StageError> {
    let mut ros = None;
    let mut iron = None;
    let mut anti = None;
    let mut damage = None;
    for (idx, def) in v1::gene_sets().iter().enumerate() {
        let bit = 1u64
            .checked_shl(idx as u32)
            .ok_or_else(|| StageError::Validation("panel bit overflow".to_string()))?;
        match def.name {
            "LYSO_ROS_SOURCES" => ros = Some(bit),
            "IRON_REDOX_CONTEXT" => iron = Some(bit),
            "ANTIOXIDANT_BUFFERING" => anti = Some(bit),
            "ROS_DAMAGE_RESPONDERS" => damage = Some(bit),
            _ => {}
        }
    }
    let ros =
        ros.ok_or_else(|| StageError::Validation("missing LYSO_ROS_SOURCES panel".to_string()))?;
    let iron =
        iron.ok_or_else(|| StageError::Validation("missing IRON_REDOX_CONTEXT panel".to_string()))?;
    let anti = anti
        .ok_or_else(|| StageError::Validation("missing ANTIOXIDANT_BUFFERING panel".to_string()))?;
    let damage = damage
        .ok_or_else(|| StageError::Validation("missing ROS_DAMAGE_RESPONDERS panel".to_string()))?;
    Ok(PanelBits {
        ros,
        iron,
        anti,
        damage,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage23_lysosomal_ros.rs"]
mod tests;
