use crate::config::thresholds::Thresholds;
use serde::{Deserialize, Serialize};

use crate::model::ctx::{
    AutophagyMetrics, AutophagySelectivityMetrics, CrossOrganelleMetrics, DrugVulnerabilityMetrics,
    LysosomalDamageMetrics, LysosomeClassMetrics, LysosomeDependencyClass, LysosomeMetrics,
    RegulatoryMetrics, SurvivalMetrics, SurvivalStats,
};

#[derive(Debug, Serialize)]
pub struct ToolInfo<'a> {
    pub name: &'a str,
    pub version: &'a str,
    pub schema: &'a str,
}

#[derive(Debug, Serialize)]
pub struct InputInfo<'a> {
    pub r#type: &'a str,
    pub mode: &'a str,
    pub n_obs: usize,
    pub samples: Vec<SampleInfo<'a>>,
}

#[derive(Debug, Serialize)]
pub struct SampleInfo<'a> {
    pub id: &'a str,
    pub sample_group: &'a str,
    pub timepoint: i32,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ThresholdsInfo {
    pub ssm_high: f32,
    pub afp_high: f32,
    pub lds_high: f32,
    pub prolif_low: f32,
    pub apop_low: f32,
    pub erdi_high: f32,
    pub lds_low: f32,
    pub afp_low: f32,
    pub ldi_high: f32,
    pub stall_high: f32,
    pub vatp_high: f32,
    pub cat_imbalance_high: f32,
    pub mito_reliance_high: f32,
    pub lsi_high: f32,
}

impl From<&Thresholds> for ThresholdsInfo {
    fn from(value: &Thresholds) -> Self {
        Self {
            ssm_high: value.ssm_high,
            afp_high: value.afp_high,
            lds_high: value.lds_high,
            prolif_low: value.prolif_low,
            apop_low: value.apop_low,
            erdi_high: value.erdi_high,
            lds_low: value.lds_low,
            afp_low: value.afp_low,
            ldi_high: value.ldi_high,
            stall_high: value.stall_high,
            vatp_high: value.vatp_high,
            cat_imbalance_high: value.cat_imbalance_high,
            mito_reliance_high: value.mito_reliance_high,
            lsi_high: value.lsi_high,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct DatasetStats<'a> {
    pub core: CoreStats,
    pub regulatory: RegulatoryStats,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub _phantom: Option<&'a str>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CoreStats {
    pub mean_afp: f32,
    pub std_afp: f32,
    pub mean_lds: f32,
    pub std_lds: f32,
    pub mean_prolif: f32,
    pub std_prolif: f32,
    pub mean_apop: f32,
    pub std_apop: f32,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct RegulatoryStats {
    pub median_tfeb: f32,
    pub median_mtor: f32,
}

impl From<&SurvivalStats> for CoreStats {
    fn from(value: &SurvivalStats) -> Self {
        Self {
            mean_afp: value.mean_afp,
            std_afp: value.std_afp,
            mean_lds: value.mean_lds,
            std_lds: value.std_lds,
            mean_prolif: value.mean_prolif,
            std_prolif: value.std_prolif,
            mean_apop: value.mean_apop,
            std_apop: value.std_apop,
        }
    }
}

impl From<&RegulatoryMetrics> for RegulatoryStats {
    fn from(value: &RegulatoryMetrics) -> Self {
        Self {
            median_tfeb: value.median_tfeb,
            median_mtor: value.median_mtor,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct Observation<'a> {
    pub id: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub autophagy: Option<AutophagyBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lysosome: Option<LysosomeBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub regulatory: Option<RegulatoryBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub survival: Option<SurvivalBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub damage: Option<DamageBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub selectivity: Option<SelectivityBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coupling: Option<CouplingBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cross_organelle: Option<CrossOrganelleBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub classification: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vulnerabilities: Option<Vec<&'a str>>,
}

#[derive(Debug, Serialize)]
pub struct AutophagyBlock {
    pub afp: f32,
    pub components: AutophagyComponents,
}

#[derive(Debug, Serialize)]
pub struct AutophagyComponents {
    pub init: f32,
    pub elong: f32,
    pub late: f32,
    pub cargo: f32,
    pub stall: f32,
}

#[derive(Debug, Serialize)]
pub struct LysosomeBlock {
    pub lds: f32,
    pub vatp: f32,
    pub protease: f32,
    pub membrane: f32,
}

#[derive(Debug, Serialize)]
pub struct RegulatoryBlock {
    pub tfeb: f32,
    pub mtor: f32,
    pub imbalance: f32,
}

#[derive(Debug, Serialize)]
pub struct SurvivalBlock {
    pub ssm: f32,
    pub z_scores: ZScores,
}

#[derive(Debug, Serialize)]
pub struct ZScores {
    pub afp: f32,
    pub lds: f32,
    pub prolif: f32,
    pub apop: f32,
}

#[derive(Debug, Serialize)]
pub struct DamageBlock {
    pub ldi: f32,
    pub lmp: f32,
    pub stress: f32,
    pub cathepsin_membrane_imbalance: f32,
}

#[derive(Debug, Serialize)]
pub struct SelectivityBlock {
    pub fractions: SelectivityFractions,
    pub entropy: f32,
}

#[derive(Debug, Serialize)]
pub struct SelectivityFractions {
    pub mito: f32,
    pub aggre: f32,
    pub er: f32,
    pub ferr: f32,
    pub lipo: f32,
}

#[derive(Debug, Serialize)]
pub struct CouplingBlock {
    pub lsi: f32,
    pub ampk_mtor_ratio: f32,
}

#[derive(Debug, Serialize)]
pub struct CrossOrganelleBlock {
    pub erdi: f32,
    pub mitophagy_reliance: f32,
}

pub fn class_to_str(class: LysosomeDependencyClass) -> &'static str {
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

pub fn vulnerability_to_str(tag: crate::model::ctx::DrugVulnerabilityTag) -> &'static str {
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

pub fn make_autophagy_block(metrics: &AutophagyMetrics, idx: usize) -> AutophagyBlock {
    AutophagyBlock {
        afp: metrics.afp[idx],
        components: AutophagyComponents {
            init: metrics.initiation[idx],
            elong: metrics.elongation[idx],
            late: metrics.degradation[idx],
            cargo: metrics.cargo[idx],
            stall: metrics.stall[idx],
        },
    }
}

pub fn make_lysosome_block(metrics: &LysosomeMetrics, idx: usize) -> LysosomeBlock {
    LysosomeBlock {
        lds: metrics.lds[idx],
        vatp: metrics.vatp[idx],
        protease: metrics.prot[idx],
        membrane: metrics.mem[idx],
    }
}

pub fn make_regulatory_block(metrics: &RegulatoryMetrics, idx: usize) -> RegulatoryBlock {
    RegulatoryBlock {
        tfeb: metrics.tfeb[idx],
        mtor: metrics.mtor[idx],
        imbalance: metrics.tfeb_mtor_diff[idx],
    }
}

pub fn make_survival_block(metrics: &SurvivalMetrics, idx: usize) -> SurvivalBlock {
    SurvivalBlock {
        ssm: metrics.ssm[idx],
        z_scores: ZScores {
            afp: metrics.z_afp[idx],
            lds: metrics.z_lds[idx],
            prolif: metrics.z_prolif[idx],
            apop: metrics.z_apop[idx],
        },
    }
}

pub fn make_damage_block(metrics: &LysosomalDamageMetrics, idx: usize) -> DamageBlock {
    DamageBlock {
        ldi: metrics.ldi[idx],
        lmp: metrics.lmp[idx],
        stress: metrics.stress[idx],
        cathepsin_membrane_imbalance: metrics.cathepsin_membrane_imbalance[idx],
    }
}

pub fn make_selectivity_block(
    metrics: &AutophagySelectivityMetrics,
    idx: usize,
) -> SelectivityBlock {
    SelectivityBlock {
        fractions: SelectivityFractions {
            mito: metrics.mito_frac[idx],
            aggre: metrics.aggre_frac[idx],
            er: metrics.er_frac[idx],
            ferr: metrics.ferr_frac[idx],
            lipo: metrics.lipo_frac[idx],
        },
        entropy: metrics.entropy[idx],
    }
}

pub fn make_coupling_block(
    metrics: &crate::model::ctx::CouplingMetrics,
    idx: usize,
) -> CouplingBlock {
    CouplingBlock {
        lsi: metrics.locked_survival_index[idx],
        ampk_mtor_ratio: metrics.ampk_mtor_ratio[idx],
    }
}

pub fn make_cross_organelle_block(
    metrics: &CrossOrganelleMetrics,
    idx: usize,
) -> CrossOrganelleBlock {
    CrossOrganelleBlock {
        erdi: metrics.energy_recycling_dependency[idx],
        mitophagy_reliance: metrics.mitophagy_reliance[idx],
    }
}

pub fn make_classification_block(metrics: &LysosomeClassMetrics, idx: usize) -> &'static str {
    class_to_str(metrics.class[idx])
}

pub fn make_vulnerability_block(
    metrics: &DrugVulnerabilityMetrics,
    idx: usize,
) -> Vec<&'static str> {
    metrics.tags[idx]
        .iter()
        .copied()
        .map(vulnerability_to_str)
        .collect()
}
