use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;

use crate::model::ctx::{Ctx, TherapyResponseClass};
use crate::stage_error::StageError;

pub fn write_cohort_therapy_summary(ctx: &Ctx, out_dir: &Path) -> Result<(), StageError> {
    let therapy = match ctx.therapy_delta.as_ref() {
        Some(value) => value,
        None => return Ok(()),
    };

    let total = therapy.response_class.len();
    if total == 0 {
        return Ok(());
    }

    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for key in all_response_keys() {
        counts.insert(key, 0);
    }
    for class in &therapy.response_class {
        let key = response_to_str(*class);
        if let Some(entry) = counts.get_mut(key) {
            *entry += 1;
        }
    }

    let total_f = total as f32;
    let mut fractions: BTreeMap<&'static str, f32> = BTreeMap::new();
    for (key, count) in counts {
        fractions.insert(key, count as f32 / total_f);
    }

    fs::create_dir_all(out_dir)?;
    let path = out_dir.join("cohort_therapy_response.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &fractions)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;

    Ok(())
}

fn response_to_str(class: TherapyResponseClass) -> &'static str {
    match class {
        TherapyResponseClass::CytotoxicResponse => "CYTOTOXIC_RESPONSE",
        TherapyResponseClass::AdaptiveSurvival => "ADAPTIVE_SURVIVAL",
        TherapyResponseClass::LysoEscape => "LYSO_ESCAPE",
        TherapyResponseClass::DamageAccumulation => "DAMAGE_ACCUMULATION",
        TherapyResponseClass::NoResponse => "NO_RESPONSE",
    }
}

fn all_response_keys() -> [&'static str; 5] {
    [
        "ADAPTIVE_SURVIVAL",
        "CYTOTOXIC_RESPONSE",
        "DAMAGE_ACCUMULATION",
        "LYSO_ESCAPE",
        "NO_RESPONSE",
    ]
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/cohort_therapy.rs"]
mod tests;
