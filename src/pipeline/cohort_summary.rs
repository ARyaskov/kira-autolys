use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;

use serde_json::Value;

use crate::model::ctx::{Ctx, LysosomeDependencyClass};
use crate::stage_error::StageError;

pub fn write_cohort_class_summary(ctx: &Ctx, out_dir: &Path) -> Result<(), StageError> {
    let classes = ctx
        .lysosome_class
        .as_ref()
        .ok_or_else(|| StageError::Validation("lysosome classes missing".to_string()))?;

    if ctx.obs_ids.len() != classes.class.len() {
        return Err(StageError::Validation(
            "obs_ids length does not match class length".to_string(),
        ));
    }

    let global_counts = count_classes(&classes.class);
    let global = fractions_from_counts(&global_counts, classes.class.len());

    let mut by_sample_group = BTreeMap::new();
    if ctx.samples.len() == classes.class.len() {
        let mut group_map: BTreeMap<String, Vec<LysosomeDependencyClass>> = BTreeMap::new();
        for (idx, sample) in ctx.samples.iter().enumerate() {
            if sample.sample_group.is_empty() {
                continue;
            }
            group_map
                .entry(sample.sample_group.clone())
                .or_default()
                .push(classes.class[idx]);
        }

        for (group, values) in group_map {
            let counts = count_classes(&values);
            let fractions = fractions_from_counts(&counts, values.len());
            by_sample_group.insert(group, fractions);
        }
    }

    let json = serde_json::json!({
        "global": global,
        "by_sample_group": if by_sample_group.is_empty() { Value::Null } else { serde_json::to_value(by_sample_group).unwrap_or(Value::Null) }
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
    for class in all_class_keys() {
        map.insert(class, 0usize);
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
#[path = "../../tests/src_inline/pipeline/cohort_summary.rs"]
mod tests;
