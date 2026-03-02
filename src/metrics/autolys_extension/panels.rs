pub const AUTOLYS_EXTENSION_PANEL_V1: &str = "AUTOLYS_EXTENSION_PANEL_V1";
pub const MIN_GENES: usize = 2;

#[derive(Debug, Clone, Copy)]
pub struct PanelDef {
    pub name: &'static str,
    pub genes: &'static [&'static [&'static str]],
}

pub const INITIATION_PANEL: PanelDef = PanelDef {
    name: "AUTOLYS_INITIATION",
    genes: &[
        &["ULK1"],
        &["ULK2"],
        &["ATG13"],
        &["RB1CC1", "FIP200"],
        &["PRKAA1"],
        &["PRKAA2"],
    ],
};

pub const FORMATION_PANEL: PanelDef = PanelDef {
    name: "AUTOLYS_FORMATION",
    genes: &[
        &["BECN1"],
        &["ATG5"],
        &["ATG7"],
        &["ATG10"],
        &["ATG12"],
        &["MAP1LC3B", "LC3B"],
        &["WIPI1"],
        &["WIPI2"],
    ],
};

pub const LYSOSOME_PANEL: PanelDef = PanelDef {
    name: "AUTOLYS_LYSOSOME",
    genes: &[&["TFEB"], &["LAMP1"], &["LAMP2"], &["CTSD"], &["CTSB"]],
};

pub const ACIDIFICATION_PANEL: PanelDef = PanelDef {
    name: "AUTOLYS_ACIDIFICATION",
    genes: &[
        &["ATP6V1A"],
        &["ATP6V1B2"],
        &["ATP6V1E1"],
        &["ATP6V0A1"],
        &["ATP6V0D1"],
    ],
};

pub const CARGO_ADAPTOR_PANEL: PanelDef = PanelDef {
    name: "AUTOLYS_CARGO_ADAPTOR",
    genes: &[&["SQSTM1", "P62"], &["OPTN"], &["NBR1"]],
};

pub const MITOPHAGY_PANEL: PanelDef = PanelDef {
    name: "AUTOLYS_MITOPHAGY",
    genes: &[
        &["PINK1"],
        &["PRKN", "PARK2"],
        &["BNIP3"],
        &["BNIP3L", "NIX"],
    ],
};
