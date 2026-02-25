use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::model::cli::RunMode;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputType {
    Mtx10x,
    DenseTsv,
    H5ad,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Cell,
    Sample,
}

#[derive(Debug, Clone)]
pub struct GeneSetResolved {
    pub name: String,
    pub indices: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct MissingGenes {
    pub panel: String,
    pub missing: Vec<String>,
}

pub use crate::expr::normalized::NormalizedExpr;

#[derive(Debug, Clone)]
pub struct AutophagyMetrics {
    pub afp: Vec<f32>,
    pub initiation: Vec<f32>,
    pub elongation: Vec<f32>,
    pub degradation: Vec<f32>,
    pub cargo: Vec<f32>,
    pub stall: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct LysosomeMetrics {
    pub lds: Vec<f32>,
    pub vatp: Vec<f32>,
    pub prot: Vec<f32>,
    pub mem: Vec<f32>,
    pub global_load: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct LysosomalDamageMetrics {
    pub lmp: Vec<f32>,
    pub stress: Vec<f32>,
    pub cathepsin_membrane_imbalance: Vec<f32>,
    pub ldi: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct AutophagySelectivityMetrics {
    pub mitophagy: Vec<f32>,
    pub aggrephagy: Vec<f32>,
    pub erphagy: Vec<f32>,
    pub ferritinophagy: Vec<f32>,
    pub lipophagy: Vec<f32>,
    pub mito_frac: Vec<f32>,
    pub aggre_frac: Vec<f32>,
    pub er_frac: Vec<f32>,
    pub ferr_frac: Vec<f32>,
    pub lipo_frac: Vec<f32>,
    pub entropy: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct CouplingMetrics {
    pub ampk: Vec<f32>,
    pub ragulator: Vec<f32>,
    pub ampk_mtor_ratio: Vec<f32>,
    pub tfeb_lyso_alignment: Vec<f32>,
    pub rag_mtor_mismatch: Vec<f32>,
    pub locked_survival_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct MitoMetrics {
    pub mito_stress: Vec<f32>,
    pub mito_decay: Vec<f32>,
    pub mito_mass_proxy: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct CrossOrganelleMetrics {
    pub mitophagy_reliance: Vec<f32>,
    pub mito_lyso_imbalance: Vec<f32>,
    pub mito_consumption_risk: Vec<f32>,
    pub energy_recycling_dependency: Vec<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TherapyResponseClass {
    CytotoxicResponse,
    AdaptiveSurvival,
    LysoEscape,
    DamageAccumulation,
    NoResponse,
}

#[derive(Debug, Clone)]
pub struct TherapyDeltaMetrics {
    pub delta_ssm: Vec<f32>,
    pub delta_afp: Vec<f32>,
    pub delta_lds: Vec<f32>,
    pub delta_ldi: Vec<f32>,
    pub delta_erdi: Vec<f32>,
    pub delta_prolif: Vec<f32>,
    pub asi: Vec<f32>,
    pub response_class: Vec<TherapyResponseClass>,
    pub from_obs: Vec<usize>,
    pub to_obs: Vec<usize>,
    pub sample_group: Vec<String>,
    pub timepoint0: Vec<i32>,
    pub timepoint1: Vec<i32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LysosomeDependencyClass {
    LysosomeDependent,
    AutophagyAdaptive,
    StalledAutophagy,
    LysosomeOverloaded,
    EnergyRecyclingDependent,
    NonLysosomal,
    Unclassified,
}

#[derive(Debug, Clone)]
pub struct LysosomeClassMetrics {
    pub class: Vec<LysosomeDependencyClass>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DrugVulnerabilityTag {
    AutophagyInhibitionSensitive,
    LysosomeAcidificationSensitive,
    VATPaseInhibitionSensitive,
    CathepsinInhibitionSensitive,
    MitophagyDisruptionSensitive,
    EnergyStressAmplificationSensitive,
    LysosomalDestabilizationSensitive,
    LowVulnerability,
}

#[derive(Debug, Clone)]
pub struct DrugVulnerabilityMetrics {
    pub tags: Vec<Vec<DrugVulnerabilityTag>>,
}

#[derive(Debug, Clone)]
pub struct LysosomePositioningMetrics {
    pub perinuclear_mean: Vec<f32>,
    pub peripheral_mean: Vec<f32>,
    pub positioning_ratio: Vec<f32>,
    pub positioning_bias: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct FerroptosisMetrics {
    pub ferritinophagy_load: Vec<f32>,
    pub iron_import_bias: Vec<f32>,
    pub ferroptosis_defense: Vec<f32>,
    pub lipid_peroxidation_context: Vec<f32>,
    pub ferroptotic_pressure_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct CholesterolTraffickingMetrics {
    pub export_capacity: Vec<f32>,
    pub import_pressure: Vec<f32>,
    pub efflux_capacity: Vec<f32>,
    pub cholesterol_trap_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct ERLysosomeContactMetrics {
    pub tethering_load: Vec<f32>,
    pub calcium_transfer_load: Vec<f32>,
    pub adaptor_stress: Vec<f32>,
    pub contact_stress_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct SecretoryLysosomeMetrics {
    pub trafficking_load: Vec<f32>,
    pub fusion_load: Vec<f32>,
    pub exocytosis_regulation: Vec<f32>,
    pub secretory_bias_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct AntigenProcessingMetrics {
    pub mhc_expression_load: Vec<f32>,
    pub protease_load: Vec<f32>,
    pub loading_accessory_load: Vec<f32>,
    pub antigen_processing_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct LipidBufferingMetrics {
    pub storage_buffering_load: Vec<f32>,
    pub utilization_load: Vec<f32>,
    pub lipophagy_context: Vec<f32>,
    pub lipid_buffering_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct AcidificationMetrics {
    pub acap: Vec<f32>,
    pub pap: Vec<f32>,
    pub acl: Vec<f32>,
    pub air: Vec<f32>,
    pub lasi: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct BiogenesisMetrics {
    pub biogenesis_drive: Vec<f32>,
    pub maturation_yield: Vec<f32>,
    pub functional_capacity: Vec<f32>,
    pub bpr: Vec<f32>,
    pub lbpi: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct MembraneRepairMetrics {
    pub rma: Vec<f32>,
    pub lpe: Vec<f32>,
    pub msr: Vec<f32>,
    pub erc: Vec<f32>,
    pub rdr: Vec<f32>,
    pub lmrci: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct CalciumCouplingMetrics {
    pub lcrc: Vec<f32>,
    pub mcuc: Vec<f32>,
    pub csap: Vec<f32>,
    pub ccb: Vec<f32>,
    pub lmcci: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct LysosomalROSMetrics {
    pub ros_generation_load: Vec<f32>,
    pub iron_redox_context: Vec<f32>,
    pub antioxidant_capacity: Vec<f32>,
    pub damage_response_load: Vec<f32>,
    pub lysosomal_ros_stress_index: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct RegulatoryMetrics {
    pub tfeb: Vec<f32>,
    pub mtor: Vec<f32>,
    pub tfeb_act: Vec<f32>,
    pub mtor_supp: Vec<f32>,
    pub tfeb_mtor_diff: Vec<f32>,
    pub tfeb_mtor_ratio: Vec<f32>,
    pub median_tfeb: f32,
    pub median_mtor: f32,
}

#[derive(Debug, Clone)]
pub struct SurvivalMetrics {
    pub prolif: Vec<f32>,
    pub apop_ready: Vec<f32>,
    pub z_afp: Vec<f32>,
    pub z_lds: Vec<f32>,
    pub z_prolif: Vec<f32>,
    pub z_apop: Vec<f32>,
    pub ssm: Vec<f32>,
    pub flag_ssm_high: Vec<bool>,
    pub flag_afp_high: Vec<bool>,
    pub flag_lds_high: Vec<bool>,
    pub flag_prolif_low: Vec<bool>,
    pub flag_apop_low: Vec<bool>,
}

#[derive(Debug, Clone, Copy)]
pub struct SurvivalStats {
    pub mean_afp: f32,
    pub std_afp: f32,
    pub mean_lds: f32,
    pub std_lds: f32,
    pub mean_prolif: f32,
    pub std_prolif: f32,
    pub mean_apop: f32,
    pub std_apop: f32,
}

#[derive(Debug, Clone)]
pub struct SampleInfo {
    pub id: String,
    pub sample_group: String,
    pub treatment: Option<String>,
    pub timepoint: i32,
    pub timepoint_label: Option<String>,
}

#[derive(Default)]
pub struct Ctx {
    pub run_mode: RunMode,
    pub input_type: Option<InputType>,
    pub mode: Option<Mode>,
    pub input_path: Option<PathBuf>,
    pub cache_path_override: Option<PathBuf>,
    pub obs_ids: Vec<String>,
    pub gene_symbols: Vec<String>,
    pub gene_index: BTreeMap<String, usize>,
    pub n_obs: usize,
    pub n_vars: usize,
    pub obs_libsize: Vec<f32>,
    pub obs_nnz: Vec<u32>,
    pub obs_expressed_genes: Vec<u32>,
    pub gene_sets: Vec<GeneSetResolved>,
    pub gene_panel_mask: Vec<u64>,
    pub missing_genes: Vec<MissingGenes>,
    pub normalized: Option<Box<dyn NormalizedExpr>>,
    pub autophagy: Option<AutophagyMetrics>,
    pub lysosome: Option<LysosomeMetrics>,
    pub lysosomal_damage: Option<LysosomalDamageMetrics>,
    pub autophagy_selectivity: Option<AutophagySelectivityMetrics>,
    pub coupling: Option<CouplingMetrics>,
    pub mito: Option<MitoMetrics>,
    pub cross_organelle: Option<CrossOrganelleMetrics>,
    pub therapy_delta: Option<TherapyDeltaMetrics>,
    pub lysosome_class: Option<LysosomeClassMetrics>,
    pub drug_vulnerability: Option<DrugVulnerabilityMetrics>,
    pub lysosome_positioning: Option<LysosomePositioningMetrics>,
    pub ferroptosis: Option<FerroptosisMetrics>,
    pub cholesterol: Option<CholesterolTraffickingMetrics>,
    pub er_lysosome_contact: Option<ERLysosomeContactMetrics>,
    pub secretory_lysosome: Option<SecretoryLysosomeMetrics>,
    pub antigen_processing: Option<AntigenProcessingMetrics>,
    pub lipid_buffering: Option<LipidBufferingMetrics>,
    pub acidification: Option<AcidificationMetrics>,
    pub biogenesis: Option<BiogenesisMetrics>,
    pub membrane_repair: Option<MembraneRepairMetrics>,
    pub calcium_coupling: Option<CalciumCouplingMetrics>,
    pub lysosomal_ros: Option<LysosomalROSMetrics>,
    pub regulatory: Option<RegulatoryMetrics>,
    pub survival: Option<SurvivalMetrics>,
    pub survival_stats: Option<SurvivalStats>,
    pub samples: Vec<SampleInfo>,
    pub thresholds: Option<crate::config::thresholds::Thresholds>,
    pub normalization: Option<serde_json::Value>,
}
