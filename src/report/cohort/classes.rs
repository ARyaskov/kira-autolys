use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;

use crate::model::cohort_view::CohortView;
use crate::model::ctx::LysosomeDependencyClass;
use crate::stage_error::StageError;

pub fn write_classes(view: &CohortView, out_dir: &Path, by_group: bool) -> Result<(), StageError> {
    let global_counts = count_classes(&view.classes);
    let global = fractions_from_counts(&global_counts, view.classes.len());

    let mut by_sample_group = BTreeMap::new();
    if by_group && !view.sample_info.is_empty() {
        let mut group_map: BTreeMap<String, Vec<LysosomeDependencyClass>> = BTreeMap::new();
        for (idx, sample) in view.sample_info.iter().enumerate() {
            if sample.sample_group.is_empty() {
                continue;
            }
            group_map
                .entry(sample.sample_group.clone())
                .or_default()
                .push(view.classes[idx]);
        }
        for (group, values) in group_map {
            let counts = count_classes(&values);
            let fractions = fractions_from_counts(&counts, values.len());
            by_sample_group.insert(group, fractions);
        }
    }

    let json = serde_json::json!({
        "global": global,
        "by_sample_group": if by_sample_group.is_empty() { serde_json::Value::Null } else { serde_json::to_value(by_sample_group).unwrap_or(serde_json::Value::Null) }
    });

    fs::create_dir_all(out_dir)?;
    let path = out_dir.join("cohort_classes.json");
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &json)
        .map_err(|err| StageError::Format(format!("json serialization failed: {err}")))?;

    Ok(())
}

fn count_classes(classes: &[LysosomeDependencyClass]) -> BTreeMap<&'static str, usize> {
    let mut map = BTreeMap::new();
    for key in all_class_keys() {
        map.insert(key, 0usize);
    }
    for class in classes {
        let key = class_to_str(*class);
        if let Some(entry) = map.get_mut(key) {
            *entry += 1;
        }
    }
    map
}

fn fractions_from_counts(
    counts: &BTreeMap<&'static str, usize>,
    total: usize,
) -> BTreeMap<&'static str, f32> {
    let mut map = BTreeMap::new();
    let total_f = if total == 0 { 1.0 } else { total as f32 };
    for (key, count) in counts {
        map.insert(*key, *count as f32 / total_f);
    }
    map
}

fn class_to_str(class: LysosomeDependencyClass) -> &'static str {
    match class {
        LysosomeDependencyClass::LysosomeDependent => "LysosomeDependent",
        LysosomeDependencyClass::AutophagyAdaptive => "AutophagyAdaptive",
        LysosomeDependencyClass::StalledAutophagy => "StalledAutophagy",
        LysosomeDependencyClass::LysosomeOverloaded => "LysosomeOverloaded",
        LysosomeDependencyClass::EnergyRecyclingDependent => "EnergyRecyclingDependent",
        LysosomeDependencyClass::NonLysosomal => "NonLysosomal",
        LysosomeDependencyClass::Unclassified => "Unclassified",
    }
}

fn all_class_keys() -> [&'static str; 7] {
    [
        "AutophagyAdaptive",
        "EnergyRecyclingDependent",
        "LysosomeDependent",
        "LysosomeOverloaded",
        "NonLysosomal",
        "StalledAutophagy",
        "Unclassified",
    ]
}

#[cfg(test)]
#[path = "../../../tests/src_inline/report/cohort/classes.rs"]
mod tests;
