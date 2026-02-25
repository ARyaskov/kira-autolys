use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;

use crate::model::cohort_view::CohortView;
use crate::model::ctx::DrugVulnerabilityTag;
use crate::stage_error::StageError;

pub fn write_vulnerabilities(view: &CohortView, out_dir: &Path) -> Result<(), StageError> {
    let n_obs = view.vulnerabilities.len();
    if n_obs == 0 {
        return Err(StageError::Validation("no vulnerability tags".to_string()));
    }

    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for key in all_tag_keys() {
        counts.insert(key, 0);
    }

    for tags in &view.vulnerabilities {
        for tag in tags {
            let key = tag_to_str(*tag);
            if let Some(entry) = counts.get_mut(key) {
                *entry += 1;
            }
        }
    }

    let denom = n_obs as f32;
    let mut fractions: BTreeMap<&'static str, f32> = BTreeMap::new();
    for (key, count) in counts {
        fractions.insert(key, count as f32 / denom);
    }

    fs::create_dir_all(out_dir)?;
    let path = out_dir.join("cohort_vulnerabilities.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &fractions)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;

    Ok(())
}

fn tag_to_str(tag: DrugVulnerabilityTag) -> &'static str {
    use crate::model::ctx::DrugVulnerabilityTag::*;
    match tag {
        AutophagyInhibitionSensitive => "AutophagyInhibitionSensitive",
        LysosomeAcidificationSensitive => "LysosomeAcidificationSensitive",
        VATPaseInhibitionSensitive => "VATPaseInhibitionSensitive",
        CathepsinInhibitionSensitive => "CathepsinInhibitionSensitive",
        MitophagyDisruptionSensitive => "MitophagyDisruptionSensitive",
        EnergyStressAmplificationSensitive => "EnergyStressAmplificationSensitive",
        LysosomalDestabilizationSensitive => "LysosomalDestabilizationSensitive",
        LowVulnerability => "LowVulnerability",
    }
}

fn all_tag_keys() -> [&'static str; 8] {
    [
        "AutophagyInhibitionSensitive",
        "CathepsinInhibitionSensitive",
        "EnergyStressAmplificationSensitive",
        "LowVulnerability",
        "LysosomalDestabilizationSensitive",
        "LysosomeAcidificationSensitive",
        "MitophagyDisruptionSensitive",
        "VATPaseInhibitionSensitive",
    ]
}

#[cfg(test)]
#[path = "../../../tests/src_inline/report/cohort/vulnerabilities.rs"]
mod tests;
