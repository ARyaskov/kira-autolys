# kira-autolys Metrics Specification

This document defines canonical metrics, formulas, and threshold constants used by `kira-autolys`.

Scope:
- core run output (`scores.tsv`, `summary.json`, `autolys.tsv`, `pipeline_step.json`)
- full/v2 output (`summary_v2.json` + `core_scores.tsv`, `damage.tsv`, `selectivity.tsv`, `coupling.tsv`, `cross_organelle.tsv`, `classification.tsv`, `vulnerabilities.tsv`, optional `therapy_delta.tsv`)
- cohort aggregates (`cohort/*.json`)

## Canonical Conventions

1. Sample axis
- All metrics are computed per observation (`obs`): a cell/barcode in `--mode cell`, or an aggregated sample in `--mode sample`.

2. Normalized expression
- For nonzero library size:
`E[g, obs] = ln(1 + (count[g, obs] / libsize[obs]) * scale)`
- `scale = 1e4`.
- Zero-libsize observations emit no per-gene values for stage computations.

3. Panel means
- For any panel `P`: `mean_P = arithmetic mean of E[g] over genes in P`.
- If no genes from panel are present: `mean_P = 0`.

4. Reused helpers/constants
- `EPS = 1e-6` where used in divisions/log safety.
- z-score:
`z(x) = (x - mean(x_dataset)) / (std(x_dataset) + EPS)`
- `std` is population standard deviation (`sqrt(m2 / n)`), not sample (`n-1`).
- Percentile helper uses sorted values with rounded index:
`idx = round((n-1)*p)`.

## Core Metrics (Stages 2-10)

### Autophagy (Stage 2)
- `init = mean(AUTOPHAGY_INIT)`
- `elong = mean(AUTOPHAGY_ELONG)`
- `late = mean(AUTOPHAGY_LATE)`
- `cargo = mean(AUTOPHAGY_CARGO)`
- `axis_early = 0.5*(init + elong)`
- `axis_late = late`
- `stall = max(axis_early - axis_late, 0)`
- `afp_raw = 0.5*axis_early + 0.3*axis_late + 0.2*cargo`
- `AFP = afp_raw - 0.25*stall`

### Lysosome Dependency (Stage 3)
- `vatp = mean(LYSO_VATP)`
- `prot = mean(LYSO_PROT)`
- `mem = mean(LYSO_MEM)`
- `global_load = mean(all expressed genes for obs)`
- `lds_raw = 0.5*vatp + 0.3*prot + 0.2*mem`
- `LDS = lds_raw / (global_load + EPS)`

### Regulatory Axis (Stage 4)
- `TFEB = mean(REG_TFEB)`
- `mTOR = mean(REG_MTOR)`
- Dataset medians: `median_tfeb`, `median_mtor`
- `tfeb_act = TFEB - median_tfeb`
- `mtor_supp = max(median_mtor - mTOR, 0)`
- `tfeb_mtor_diff = TFEB - mTOR`
- `tfeb_mtor_ratio = (tfeb_act + EPS) / (abs(mTOR - median_mtor) + EPS)`

### Survival Stress Mode (Stage 5)
- `prolif = mean(PROLIF)`
- `apop_ready = mean(APOP_PRO) - mean(APOP_ANTI)`
- `z_AFP = z(AFP)`, `z_LDS = z(LDS)`, `z_PROLIF = z(prolif)`, `z_APOP = z(apop_ready)`
- `SSM = 0.35*z_AFP + 0.35*z_LDS - 0.15*z_PROLIF - 0.15*z_APOP`

Flags:
- `flag_SSM_high = SSM >= ssm_high`
- `flag_AFP_high = z_AFP >= afp_high`
- `flag_LDS_high = z_LDS >= lds_high`
- `flag_PROLIF_low = z_PROLIF <= prolif_low`
- `flag_APOP_low = z_APOP <= apop_low`

### Lysosomal Damage (Stage 7)
- `LMP = mean(LYSO_LMP)`
- `stress = mean(LYSO_STRESS)`
- `cathepsin_membrane_imbalance = prot / (mem + EPS)`
- `LDI = 0.4*LMP + 0.3*stress + 0.3*cathepsin_membrane_imbalance`

### Selective Autophagy (Stage 8)
- Raw means:
`mitophagy`, `aggrephagy`, `erphagy`, `ferritinophagy`, `lipophagy`
- `total = mito + aggre + er + ferr + lipo + EPS`
- Fractions: each component `/ total`
- `entropy = -sum_i(frac_i * ln(frac_i + EPS))`

### Coupling (Stage 9)
- `ampk = mean(REG_AMPK)`
- `ragulator = mean(REG_RAGULATOR)`
- `ampk_mtor_ratio = ampk / (mTOR + EPS)`
- `tfeb_lyso_alignment = TFEB * LDS`
- `rag_mtor_mismatch = ragulator - mTOR`
- `LSI = 0.35*ampk_mtor_ratio + 0.35*tfeb_lyso_alignment + 0.30*max(rag_mtor_mismatch, 0)`

### Cross-organelle (Stage 10)
- `mitophagy_reliance = mitophagy * LDS`
- If mito block exists:
`mito_lyso_imbalance = mito_stress / (LDS + EPS)`, `mito_consumption_risk = mitophagy_reliance / (mito_mass_proxy + EPS)`;
otherwise both are `NaN`.
- `ERDI = 0.4*mitophagy_reliance + 0.3*max(LSI,0) + 0.3*LDI`

## Classification and Vulnerabilities

### Lysosome Dependency Class (Stage 12, ordered rules)
1. `ERDI >= erdi_high` -> `EnergyRecyclingDependent`
2. `LDS >= lds_high && LDI >= ldi_high` -> `LysosomeOverloaded`
3. `flag_SSM_high && AFP >= afp_high && LDI < ldi_high` -> `AutophagyAdaptive`
4. `stall >= stall_high` -> `StalledAutophagy`
5. `LDS >= lds_high && !flag_SSM_high` -> `LysosomeDependent`
6. `LDS < lds_low && AFP < afp_low` -> `NonLysosomal`
7. else `Unclassified`

### Drug Vulnerability Tags (Stage 13)
- `AutophagyInhibitionSensitive`: class `AutophagyAdaptive` OR (`flag_SSM_high` and `AFP >= afp_high`)
- `LysosomeAcidificationSensitive` + `VATPaseInhibitionSensitive`: `vatp >= vatp_high && LDS >= lds_high`
- `CathepsinInhibitionSensitive`: `cathepsin_membrane_imbalance >= cat_imbalance_high`
- `MitophagyDisruptionSensitive`: `mitophagy_reliance >= mito_reliance_high`
- `EnergyStressAmplificationSensitive`: `LSI >= lsi_high`
- `LysosomalDestabilizationSensitive`: `LDI >= ldi_high && LDS >= lds_high`
- If no tags matched: `LowVulnerability`

## Timecourse Metric (Stage 11, optional)

For consecutive timepoints inside each `sample_group`:
- `delta_* = value_t1 - value_t0` for `SSM`, `AFP`, `LDS`, `LDI`, `ERDI`, `prolif`
- `ASI = 0.4*delta_SSM + 0.3*delta_AFP + 0.3*delta_LDS`

Response classes:
- `CytotoxicResponse`: `delta_prolif < 0 && delta_SSM < 0`
- `AdaptiveSurvival`: `delta_prolif < 0 && delta_SSM > 0`
- `LysoEscape`: `delta_LDS > 0 && delta_AFP > 0 && delta_prolif <= 0`
- `DamageAccumulation`: `delta_LDI > 0 && delta_SSM <= 0`
- else `NoResponse`

## Extended Metrics (Stages 16-26, non-fast mode)

- Stage 16 (positioning):
`positioning_ratio = (perinuclear_mean+EPS)/(peripheral_mean+EPS)`, `positioning_bias = log2(positioning_ratio)`.
- Stage 17 (ferroptosis):
`iron_import_bias = TFRC/(SLC40A1+EPS)`,
`ferroptotic_pressure_index = (ferritinophagy_load * iron_import_bias * lipid_peroxidation_context)/(ferroptosis_defense+EPS)`.
- Stage 18 (cholesterol):
`cholesterol_trap_index = import_pressure/(export_capacity + efflux_capacity + EPS)`.
- Stage 19 (ER-lysosome):
`contact_stress_index = tethering_load * calcium_transfer_load * adaptor_stress`.
- Stage 20 (secretory lysosome):
`secretory_bias_index = trafficking_load * fusion_load * exocytosis_regulation`.
- Stage 21 (antigen processing):
`antigen_processing_index = mhc_expression_load * protease_load * loading_accessory_load`.
- Stage 22 (lipid buffering):
`lipid_buffering_index = (storage_buffering_load * lipophagy_context)/(utilization_load + EPS)`.
- Stage 23 (lysosomal ROS):
`lysosomal_ros_stress_index = (ros_generation_load * iron_redox_context)/(antioxidant_capacity + EPS)`.
- Stage 24 (calcium coupling):
`CCB = LCRC/(MCUC+EPS)`,
`LMCCI = z(LCRC) + z(MCUC) + z(CSAP) - abs(z(CCB))`.
- Stage 25 (biogenesis):
`functional_capacity = LDS * (1 - z(LDI)) * (1 + z(LASI))` (if LASI unavailable, `z(LASI)=0`),
`BPR = biogenesis_drive/(functional_capacity+EPS)`,
`LBPI = z(BPR) + z(biogenesis_drive) - z(maturation_yield)`.
- Stage 26 (membrane repair):
`RMA = mean(ESCRT + Ca-recruit panels)`,
`ERC = RMA + LPE + MSR`,
`RDR = ERC/(LDI+EPS)`,
`LMRCI = z(RDR) + z(ERC) - z(LDI)`.

## Autolys Extension Metrics (Stage 2 additive extension)

The autolys extension derives single-sample compatible, deterministic proxy metrics from one normalized matrix (`E[g,obs]`), using concise mechanism-grounded panels:
- initiation (ULK/AMPK gate)
- autophagosome formation (ATG conjugation)
- lysosomal identity and acidification
- cargo adaptor load
- optional mitophagy context

Panel summarization:
- 10% trimmed mean per observation/panel.
- If expressed genes in panel `< MIN_GENES` (2), panel core is `NaN`.

Robust normalization:
- `Z = (x - median) / (1.4826*MAD + EPS)`.
- If `MAD == 0`, finite values map to `0`.

Per-observation extension metrics:
- `AIS = Z(init_core)` (initiation engagement)
- `AFS = Z(form_core)` (autophagosome machinery engagement)
- `LCI = 0.6*Z(lyso_core) + 0.4*Z(acid_core)` (lysosomal capacity)
- `Drive = 0.5*AIS + 0.5*AFS`
- `Throughput = LCI`
- `AFP_ext = Drive - Throughput` (autophagy flux proxy imbalance)
- `CDS = max(0, Z(cargo_core) - LCI)` (clearance deficit proxy)
- `ASM = max(0, 0.4*max(0,AFP_ext) + 0.3*CDS + 0.3*max(0,AIS))` (autolys stress mode)

Optional mitophagy axis:
- `MitophagyScore = Z(mito_core)` when panel genes are mappable; otherwise `NaN`.

Flags:
- `initiation_high`: `AIS >= 2.0`
- `formation_high`: `AFS >= 2.0`
- `lysosome_high`: `LCI >= 2.0`
- `bottleneck_risk`: `AFP_ext >= 1.5`
- `clearance_deficit`: `CDS >= 1.5`
- `autolys_stress_mode`: `ASM >= 2.0`
- `mitophagy_high`: `MitophagyScore >= 2.0`

NaN policy:
- NaN scores do not trigger flags.
- Missingness is reported in `summary.json` / `summary_v2.json` under `autolys_extension`.

Caveats:
- `AFP_ext` is a transcriptomic proxy and does **not** measure biochemical flux directly.
- Interpret jointly with proteostasis/clearance outputs (for example `kira-proteoqc`) and mitochondrial stress/damage context (for example `kira-mitoqc`).

## Pipeline Summary Metrics (`autolys.tsv`, Stage 6)

Raw composites:
- `raw_catabolic = LDS + AFP - prolif`
- `raw_recycling = AFP - stall - LDI` (if damage missing in fast path, `LDI=0`)

Quantile-normalized (`p10..p90`, clamped to `[0,1]`, fallback `0.5` when `hi<=lo`):
- `lysosomal_load = Qnorm(LDS)`
- `autophagy_flux_proxy = Qnorm(AFP)`
- `catabolic_bias = Qnorm(raw_catabolic)`
- `recycling_efficiency = Qnorm(raw_recycling)`
- `stress_autophagy_index = Qnorm(SSM)`

Confidence:
- start `1.0`
- `-0.45` if `libsize<=0` or `nnz<20` (`LOW_CONFIDENCE`)
- `-0.35` if `lysosomal_load<0.15` (`LOW_LYSO_SIGNAL`)
- lower-bounded by `0`

Regime rules:
1. `confidence < 0.20` -> `Unclassified`
2. `stress_autophagy_index >= 0.75 && catabolic_bias >= 0.75` -> `CatabolicOverdrive`
3. `stress_autophagy_index >= 0.60 && recycling_efficiency < 0.30` -> `AutophagyImpaired`
4. `lysosomal_load < 0.20 && recycling_efficiency < 0.30` -> `DegenerativeLysosomal`
5. `stress_autophagy_index >= 0.60 || autophagy_flux_proxy >= 0.70` -> `StressInducedAutophagy`
6. `recycling_efficiency >= 0.40` -> `HomeostaticRecycling`
7. else `Unclassified`

## Cohort Aggregates

- `cohort_metrics.json` (`SSM`, `AFP`, `LDS`, `LDI`, `LSI`, `ERDI`):
`mean`, `median` (`p50` via rounded-index percentile), `p90`, `fraction_high = count(value>=threshold)/n`.
- `cohort_signatures.json`:
`SSM_high`, `LDI_LDS_high`, `AutophagyAdaptive_ERDI_high` as fractions of observations.
- `cohort_classes.json`:
class fractions globally and optionally by `sample_group`.
- `cohort_vulnerabilities.json`:
per-tag frequency normalized by number of observations.
- `cohort_therapy_response.json`:
therapy response-class fractions.

## Threshold Constants (Default)

From `src/config/thresholds.rs`:
- `ssm_high = 1.0`
- `afp_high = 0.8`
- `lds_high = 0.8`
- `prolif_low = -0.5`
- `apop_low = -0.5`
- `erdi_high = 1.0`
- `lds_low = 0.3`
- `afp_low = 0.3`
- `ldi_high = 1.0`
- `stall_high = 0.5`
- `vatp_high = 1.0`
- `cat_imbalance_high = 1.2`
- `mito_reliance_high = 1.0`
- `lsi_high = 1.0`
