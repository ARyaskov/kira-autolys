use std::collections::BTreeMap;

use serde_json::{Map, Value, json};

use crate::metrics::autolys_extension::panels::AUTOLYS_EXTENSION_PANEL_V1;
use crate::metrics::autolys_extension::scores::median_from_slice;
use crate::model::ctx::{AutolysExtensionMetrics, Ctx};

const INITIATION_HIGH: f32 = 2.0;
const FORMATION_HIGH: f32 = 2.0;
const LYSOSOME_HIGH: f32 = 2.0;
const BOTTLENECK_RISK: f32 = 1.5;
const CLEARANCE_DEFICIT: f32 = 1.5;
const AUTOLYS_STRESS_MODE: f32 = 2.0;
const MITOPHAGY_HIGH: f32 = 2.0;

pub fn build_summary_block(ctx: &Ctx, metrics: &AutolysExtensionMetrics) -> Value {
    json!({
        "panel_version": AUTOLYS_EXTENSION_PANEL_V1,
        "thresholds": {
            "initiation_high": INITIATION_HIGH,
            "formation_high": FORMATION_HIGH,
            "lysosome_high": LYSOSOME_HIGH,
            "bottleneck_risk": BOTTLENECK_RISK,
            "clearance_deficit": CLEARANCE_DEFICIT,
            "autolys_stress_mode": AUTOLYS_STRESS_MODE,
            "mitophagy_high": MITOPHAGY_HIGH,
        },
        "global_stats": global_stats(metrics),
        "cluster_stats": cluster_stats(ctx, metrics),
        "missingness": missingness(metrics),
    })
}

fn global_stats(metrics: &AutolysExtensionMetrics) -> Value {
    let mut map = Map::new();
    insert_metric_stats(&mut map, "AIS", &metrics.ais);
    insert_metric_stats(&mut map, "AFS", &metrics.afs);
    insert_metric_stats(&mut map, "LCI", &metrics.lci);
    insert_metric_stats(&mut map, "AFP", &metrics.afp);
    insert_metric_stats(&mut map, "CDS", &metrics.cds);
    insert_metric_stats(&mut map, "ASM", &metrics.asm);
    insert_metric_stats(&mut map, "MitophagyScore", &metrics.mitophagy_score);
    Value::Object(map)
}

fn insert_metric_stats(map: &mut Map<String, Value>, key: &str, values: &[f32]) {
    let finite: Vec<f32> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    let median = median_from_slice(&finite);
    let mad = mad(&finite, median);
    map.insert(
        key.to_string(),
        json!({
            "median": round6(median),
            "mad": round6(mad),
        }),
    );
}

fn cluster_stats(ctx: &Ctx, metrics: &AutolysExtensionMetrics) -> Value {
    let groups = build_groups(ctx);
    let mut cluster_map = Map::new();
    let mut cds_rank = Vec::new();
    let mut asm_rank = Vec::new();

    for (group, indices) in groups {
        let ais = values_for_indices(&metrics.ais, &indices);
        let afs = values_for_indices(&metrics.afs, &indices);
        let lci = values_for_indices(&metrics.lci, &indices);
        let afp = values_for_indices(&metrics.afp, &indices);
        let cds = values_for_indices(&metrics.cds, &indices);
        let asm = values_for_indices(&metrics.asm, &indices);
        let mito = values_for_indices(&metrics.mitophagy_score, &indices);

        let cluster_entry = json!({
            "n_obs": indices.len(),
            "metrics": {
                "AIS": metric_distribution(&ais),
                "AFS": metric_distribution(&afs),
                "LCI": metric_distribution(&lci),
                "AFP": metric_distribution(&afp),
                "CDS": metric_distribution(&cds),
                "ASM": metric_distribution(&asm),
                "MitophagyScore": metric_distribution(&mito),
            },
            "flag_fractions": {
                "initiation_high": flag_fraction(&metrics.initiation_high, &indices),
                "formation_high": flag_fraction(&metrics.formation_high, &indices),
                "lysosome_high": flag_fraction(&metrics.lysosome_high, &indices),
                "bottleneck_risk": flag_fraction(&metrics.bottleneck_risk, &indices),
                "clearance_deficit": flag_fraction(&metrics.clearance_deficit, &indices),
                "autolys_stress_mode": flag_fraction(&metrics.autolys_stress_mode, &indices),
                "mitophagy_high": flag_fraction(&metrics.mitophagy_high, &indices),
            }
        });

        cds_rank.push((group.clone(), percentile(&cds, 0.5)));
        asm_rank.push((group.clone(), percentile(&asm, 0.5)));
        cluster_map.insert(group, cluster_entry);
    }

    cds_rank.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    asm_rank.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    json!({
        "by_cluster": cluster_map,
        "top_clusters": {
            "clearance_deficit": cds_rank.iter().take(3).map(|(name, value)| {
                json!({"cluster": name, "median_cds": round6(*value)})
            }).collect::<Vec<_>>(),
            "autolys_stress_mode": asm_rank.iter().take(3).map(|(name, value)| {
                json!({"cluster": name, "median_asm": round6(*value)})
            }).collect::<Vec<_>>(),
        }
    })
}

fn metric_distribution(values: &[f32]) -> Value {
    json!({
        "median": round6(percentile(values, 0.5)),
        "p10": round6(percentile(values, 0.1)),
        "p90": round6(percentile(values, 0.9)),
    })
}

fn missingness(metrics: &AutolysExtensionMetrics) -> Value {
    json!({
        "init_core": missing_fraction(&metrics.init_core),
        "form_core": missing_fraction(&metrics.form_core),
        "lyso_core": missing_fraction(&metrics.lyso_core),
        "acid_core": missing_fraction(&metrics.acid_core),
        "cargo_core": missing_fraction(&metrics.cargo_core),
        "AIS": missing_fraction(&metrics.ais),
        "AFS": missing_fraction(&metrics.afs),
        "LCI": missing_fraction(&metrics.lci),
        "AFP": missing_fraction(&metrics.afp),
        "CDS": missing_fraction(&metrics.cds),
        "ASM": missing_fraction(&metrics.asm),
        "mito_core": missing_fraction(&metrics.mito_core),
        "MitophagyScore": missing_fraction(&metrics.mitophagy_score),
    })
}

fn missing_fraction(values: &[f32]) -> Value {
    let missing = values.iter().filter(|value| !value.is_finite()).count();
    let n = values.len().max(1) as f32;
    json!({
        "missing_count": missing,
        "missing_fraction": round6(missing as f32 / n),
    })
}

fn build_groups(ctx: &Ctx) -> BTreeMap<String, Vec<usize>> {
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    if ctx.samples.is_empty() {
        groups.insert("all".to_string(), (0..ctx.n_obs).collect());
        return groups;
    }
    let mut sample_group_by_id: BTreeMap<&str, &str> = BTreeMap::new();
    for sample in &ctx.samples {
        sample_group_by_id.insert(sample.id.as_str(), sample.sample_group.as_str());
    }
    for idx in 0..ctx.n_obs {
        let key = sample_group_by_id
            .get(ctx.obs_ids[idx].as_str())
            .copied()
            .unwrap_or("unassigned");
        groups.entry(key.to_string()).or_default().push(idx);
    }
    groups
}

fn values_for_indices(values: &[f32], indices: &[usize]) -> Vec<f32> {
    let mut out = Vec::with_capacity(indices.len());
    for idx in indices {
        let value = values.get(*idx).copied().unwrap_or(f32::NAN);
        if value.is_finite() {
            out.push(value);
        }
    }
    out
}

fn flag_fraction(flags: &[bool], indices: &[usize]) -> f32 {
    if indices.is_empty() {
        return 0.0;
    }
    let count = indices
        .iter()
        .copied()
        .filter(|idx| flags.get(*idx).copied().unwrap_or(false))
        .count();
    round6(count as f32 / indices.len() as f32)
}

fn mad(values: &[f32], median: f32) -> f32 {
    if values.is_empty() || !median.is_finite() {
        return f32::NAN;
    }
    let mut deviations = Vec::with_capacity(values.len());
    for value in values {
        deviations.push((value - median).abs());
    }
    median_from_slice(&deviations)
}

fn percentile(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let pos = ((data.len() - 1) as f32 * p.clamp(0.0, 1.0)).round() as usize;
    data[pos]
}

fn round6(value: f32) -> f32 {
    if value.is_finite() {
        (value * 1_000_000.0).round() / 1_000_000.0
    } else {
        value
    }
}
