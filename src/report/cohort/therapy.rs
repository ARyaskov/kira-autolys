use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;
use tracing::warn;

use crate::model::cohort_view::CohortView;
use crate::model::ctx::TherapyResponseClass;
use crate::stage_error::StageError;

pub fn write_therapy(view: &CohortView, out_dir: &Path) -> Result<(), StageError> {
    let responses = match view.therapy_response.as_ref() {
        Some(value) => value,
        None => {
            warn!("therapy_delta.tsv missing; skipping cohort therapy summary");
            return Ok(());
        }
    };
    if responses.is_empty() {
        return Ok(());
    }

    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for key in all_response_keys() {
        counts.insert(key, 0);
    }
    for class in responses {
        let key = response_to_str(*class);
        if let Some(entry) = counts.get_mut(key) {
            *entry += 1;
        }
    }

    let total = responses.len() as f32;
    let mut fractions: BTreeMap<&'static str, f32> = BTreeMap::new();
    for (key, count) in counts {
        fractions.insert(key, count as f32 / total);
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
#[path = "../../../tests/src_inline/report/cohort/therapy.rs"]
mod tests;
