/// Biological thresholds used for classification, flags, and cohort summaries.
///
/// These values are stable contracts and must not be changed without a
/// documented version bump.
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    /// SSM z-score threshold marking a high stress-survival mode.
    pub ssm_high: f32,
    /// AFP z-score threshold marking high autophagy flux proxy.
    pub afp_high: f32,
    /// LDS z-score threshold marking high lysosome dependency.
    pub lds_high: f32,
    /// Proliferation z-score threshold marking low proliferation.
    pub prolif_low: f32,
    /// Apoptosis readiness z-score threshold marking low apoptotic readiness.
    pub apop_low: f32,
    /// ERDI threshold marking high energy recycling dependency.
    pub erdi_high: f32,
    /// LDS absolute threshold marking low lysosome dependency.
    pub lds_low: f32,
    /// AFP absolute threshold marking low autophagy activity.
    pub afp_low: f32,
    /// LDI threshold marking high lysosomal damage.
    pub ldi_high: f32,
    /// AFP stall threshold marking stalled autophagy.
    pub stall_high: f32,
    /// VATP threshold marking high acidification machinery signal.
    pub vatp_high: f32,
    /// Cathepsin/membrane imbalance threshold for protease-dominant stress.
    pub cat_imbalance_high: f32,
    /// Mitophagy reliance threshold for cross-organelle dependency.
    pub mito_reliance_high: f32,
    /// Locked survival index threshold for regulatory lock-in.
    pub lsi_high: f32,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            ssm_high: 1.0,
            afp_high: 0.8,
            lds_high: 0.8,
            prolif_low: -0.5,
            apop_low: -0.5,
            erdi_high: 1.0,
            lds_low: 0.3,
            afp_low: 0.3,
            ldi_high: 1.0,
            stall_high: 0.5,
            vatp_high: 1.0,
            cat_imbalance_high: 1.2,
            mito_reliance_high: 1.0,
            lsi_high: 1.0,
        }
    }
}
