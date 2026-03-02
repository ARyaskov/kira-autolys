use std::collections::BTreeMap;

use crate::metrics::autolys_extension::panels::{
    ACIDIFICATION_PANEL, CARGO_ADAPTOR_PANEL, FORMATION_PANEL, INITIATION_PANEL, LYSOSOME_PANEL,
    MIN_GENES, MITOPHAGY_PANEL,
};
use crate::model::ctx::{AutolysExtensionMetrics, Ctx, NormalizedExpr};

const PANEL_COUNT: usize = 6;
const EPS: f32 = 1e-6;

#[derive(Clone)]
struct ResolvedPanels {
    panel_indices: [Vec<usize>; PANEL_COUNT],
    gene_to_panel_mask: Vec<u8>,
    mito_enabled: bool,
}

pub fn compute(ctx: &Ctx) -> AutolysExtensionMetrics {
    let Some(normalized) = ctx.normalized.as_ref() else {
        return AutolysExtensionMetrics::default_for_n_obs(ctx.n_obs);
    };
    let resolved = resolve_panels(ctx);
    let mut metrics = AutolysExtensionMetrics::default_for_n_obs(ctx.n_obs);

    compute_panel_cores(
        normalized.as_ref(),
        &resolved,
        &mut metrics.init_core,
        &mut metrics.form_core,
        &mut metrics.lyso_core,
        &mut metrics.acid_core,
        &mut metrics.cargo_core,
        &mut metrics.mito_core,
    );

    let z_init = robust_z(&metrics.init_core);
    let z_form = robust_z(&metrics.form_core);
    let z_lyso = robust_z(&metrics.lyso_core);
    let z_acid = robust_z(&metrics.acid_core);
    let z_cargo = robust_z(&metrics.cargo_core);
    let z_mito = robust_z(&metrics.mito_core);

    for idx in 0..ctx.n_obs {
        metrics.ais[idx] = z_init[idx];
        metrics.afs[idx] = z_form[idx];
        metrics.lci[idx] = weighted2(z_lyso[idx], 0.6, z_acid[idx], 0.4);
        let drive = weighted2(metrics.ais[idx], 0.5, metrics.afs[idx], 0.5);
        metrics.afp[idx] = nan_sub(drive, metrics.lci[idx]);
        metrics.cds[idx] = nonnegative_or_nan(nan_sub(z_cargo[idx], metrics.lci[idx]));
        let asm = weighted3(
            nonnegative_or_nan(metrics.afp[idx]),
            0.4,
            metrics.cds[idx],
            0.3,
            nonnegative_or_nan(metrics.ais[idx]),
            0.3,
        );
        metrics.asm[idx] = nonnegative_or_nan(asm);
        if resolved.mito_enabled {
            metrics.mitophagy_score[idx] = z_mito[idx];
        }

        metrics.initiation_high[idx] = metrics.ais[idx].is_finite() && metrics.ais[idx] >= 2.0;
        metrics.formation_high[idx] = metrics.afs[idx].is_finite() && metrics.afs[idx] >= 2.0;
        metrics.lysosome_high[idx] = metrics.lci[idx].is_finite() && metrics.lci[idx] >= 2.0;
        metrics.bottleneck_risk[idx] = metrics.afp[idx].is_finite() && metrics.afp[idx] >= 1.5;
        metrics.clearance_deficit[idx] = metrics.cds[idx].is_finite() && metrics.cds[idx] >= 1.5;
        metrics.autolys_stress_mode[idx] = metrics.asm[idx].is_finite() && metrics.asm[idx] >= 2.0;
        metrics.mitophagy_high[idx] =
            metrics.mitophagy_score[idx].is_finite() && metrics.mitophagy_score[idx] >= 2.0;
    }

    metrics
}

fn resolve_panels(ctx: &Ctx) -> ResolvedPanels {
    let mut symbol_to_idx = BTreeMap::new();
    for (symbol, idx) in &ctx.gene_index {
        symbol_to_idx.insert(symbol.to_ascii_uppercase(), *idx);
    }

    let panel_defs = [
        INITIATION_PANEL,
        FORMATION_PANEL,
        LYSOSOME_PANEL,
        ACIDIFICATION_PANEL,
        CARGO_ADAPTOR_PANEL,
        MITOPHAGY_PANEL,
    ];
    let mut panel_indices: [Vec<usize>; PANEL_COUNT] = std::array::from_fn(|_| Vec::new());
    for panel_idx in 0..PANEL_COUNT {
        for aliases in panel_defs[panel_idx].genes {
            if let Some(found_idx) = aliases
                .iter()
                .find_map(|symbol| symbol_to_idx.get(&symbol.to_ascii_uppercase()).copied())
            {
                panel_indices[panel_idx].push(found_idx);
            }
        }
        panel_indices[panel_idx].sort_unstable();
        panel_indices[panel_idx].dedup();
    }

    let mut gene_to_panel_mask = vec![0u8; ctx.n_vars];
    for (panel_idx, indices) in panel_indices.iter().enumerate() {
        let bit = 1u8.checked_shl(panel_idx as u32).unwrap_or(0);
        for &gene_idx in indices {
            if gene_idx < gene_to_panel_mask.len() {
                gene_to_panel_mask[gene_idx] |= bit;
            }
        }
    }

    let mito_enabled = !panel_indices[5].is_empty();
    ResolvedPanels {
        panel_indices,
        gene_to_panel_mask,
        mito_enabled,
    }
}

fn compute_panel_cores(
    normalized: &dyn NormalizedExpr,
    resolved: &ResolvedPanels,
    init_core: &mut [f32],
    form_core: &mut [f32],
    lyso_core: &mut [f32],
    acid_core: &mut [f32],
    cargo_core: &mut [f32],
    mito_core: &mut [f32],
) {
    let mut buffers: [Vec<f32>; PANEL_COUNT] =
        std::array::from_fn(|idx| Vec::with_capacity(resolved.panel_indices[idx].len().max(2)));
    for obs_idx in 0..normalized.n_obs() {
        for panel in &mut buffers {
            panel.clear();
        }

        normalized.for_each_in_obs(obs_idx, &mut |gene_idx, value| {
            if gene_idx >= resolved.gene_to_panel_mask.len() {
                return;
            }
            let mut mask = resolved.gene_to_panel_mask[gene_idx];
            while mask != 0 {
                let panel_idx = mask.trailing_zeros() as usize;
                buffers[panel_idx].push(value);
                mask &= mask - 1;
            }
        });

        init_core[obs_idx] = trimmed_mean(&mut buffers[0], MIN_GENES);
        form_core[obs_idx] = trimmed_mean(&mut buffers[1], MIN_GENES);
        lyso_core[obs_idx] = trimmed_mean(&mut buffers[2], MIN_GENES);
        acid_core[obs_idx] = trimmed_mean(&mut buffers[3], MIN_GENES);
        cargo_core[obs_idx] = trimmed_mean(&mut buffers[4], MIN_GENES);
        mito_core[obs_idx] = trimmed_mean(&mut buffers[5], MIN_GENES);
    }
}

fn trimmed_mean(values: &mut [f32], min_genes: usize) -> f32 {
    if values.len() < min_genes {
        return f32::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let trim = ((values.len() as f32) * 0.10).floor() as usize;
    let start = trim.min(values.len());
    let end = values.len().saturating_sub(trim);
    if start >= end {
        return f32::NAN;
    }
    let mut sum = 0.0_f32;
    for value in values.iter().take(end).skip(start) {
        sum += *value;
    }
    sum / (end - start) as f32
}

pub fn robust_z(values: &[f32]) -> Vec<f32> {
    let finite: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return vec![f32::NAN; values.len()];
    }
    let median = median_from_slice(&finite);
    let mut abs_dev = Vec::with_capacity(finite.len());
    for value in &finite {
        abs_dev.push((value - median).abs());
    }
    let mad = median_from_slice(&abs_dev);
    if mad <= EPS {
        return values
            .iter()
            .map(|value| if value.is_finite() { 0.0 } else { f32::NAN })
            .collect();
    }
    let denom = 1.4826_f32 * mad + EPS;
    values
        .iter()
        .map(|value| {
            if value.is_finite() {
                (value - median) / denom
            } else {
                f32::NAN
            }
        })
        .collect()
}

pub fn median_from_slice(values: &[f32]) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) * 0.5
    } else {
        sorted[mid]
    }
}

fn nan_sub(lhs: f32, rhs: f32) -> f32 {
    if lhs.is_finite() && rhs.is_finite() {
        lhs - rhs
    } else {
        f32::NAN
    }
}

fn weighted2(a: f32, wa: f32, b: f32, wb: f32) -> f32 {
    if a.is_finite() && b.is_finite() {
        wa * a + wb * b
    } else {
        f32::NAN
    }
}

fn weighted3(a: f32, wa: f32, b: f32, wb: f32, c: f32, wc: f32) -> f32 {
    if a.is_finite() && b.is_finite() && c.is_finite() {
        wa * a + wb * b + wc * c
    } else {
        f32::NAN
    }
}

fn nonnegative_or_nan(value: f32) -> f32 {
    if value.is_finite() {
        value.max(0.0)
    } else {
        f32::NAN
    }
}
