#[derive(Debug, Clone, Copy)]
pub struct GeneSetDef {
    pub name: &'static str,
    pub genes: &'static [&'static str],
}

const AUTOPHAGY_INIT: &[&str] = &["RB1CC1", "BECN1", "ATG13", "ULK1"];
const AUTOPHAGY_ELONG: &[&str] = &["ATG5", "ATG7", "ATG12", "ATG16L1"];
const AUTOPHAGY_LATE: &[&str] = &["LAMP1", "LAMP2", "CTSB", "CTSD"];
const AUTOPHAGY_CARGO: &[&str] = &["SQSTM1", "NBR1", "OPTN"];
const LYSO_VATP: &[&str] = &["ATP6V1A", "ATP6V1B2", "ATP6V0D1", "ATP6V0A1"];
const LYSO_PROT: &[&str] = &["CTSB", "CTSD", "CTSL", "CTSS", "LGMN"];
const LYSO_MEM: &[&str] = &["LAMP1", "LAMP2", "MCOLN1", "NPC1", "NPC2"];
const LYSO_LMP: &[&str] = &["LGALS3", "LGALS8"];
const LYSO_STRESS: &[&str] = &["HSPA1A", "HSPB1", "ATF4", "DDIT3"];
const AUTO_MITOPHAGY: &[&str] = &["BNIP3", "BNIP3L", "FUNDC1", "PINK1", "PRKN"];
const AUTO_AGGREPHAGY: &[&str] = &["SQSTM1", "NBR1", "TAX1BP1", "OPTN", "CALCOCO2"];
const AUTO_ERPHAGY: &[&str] = &["FAM134B", "RTN3", "CCPG1", "TEX264"];
const AUTO_FERRITINOPHAGY: &[&str] = &["NCOA4"];
const AUTO_LIPOPHAGY: &[&str] = &["PNPLA2", "ATGL", "LIPE"];
const PERINUCLEAR_LYSOSOME: &[&str] = &["RILP", "TBC1D15", "TBC1D5", "DYNC1H1"];
const PERIPHERAL_LYSOSOME: &[&str] = &["ARL8B", "KIF5B", "KIF1A", "BORCS5", "BORCS6"];
const FERRITINOPHAGY: &[&str] = &["NCOA4", "FTH1", "FTL"];
const FERROPTOSIS_DEFENSE: &[&str] = &["GPX4", "SLC7A11", "GCLC", "GCLM"];
const IRON_IMPORT_EXPORT: &[&str] = &["TFRC", "SLC40A1"];
const LIPID_PEROXIDATION_CONTEXT: &[&str] = &["ACSL4", "LPCAT3"];
const CHOLESTEROL_EXPORT: &[&str] = &["NPC1", "NPC2", "STARD3"];
const CHOLESTEROL_IMPORT: &[&str] = &["LDLR", "SCARB1"];
const CHOLESTEROL_EFFLUX: &[&str] = &["ABCA1", "ABCG1"];
const ER_LYSO_CONTACT_TETHERING: &[&str] = &["VAPA", "VAPB", "ORP1L", "STARD3"];
const ER_LYSO_CALCIUM_AXIS: &[&str] = &["MCOLN1", "TPCN1", "TPCN2"];
const ER_LYSO_STRESS_ADAPTORS: &[&str] = &["PDZD8", "TMEM55B"];
const SECRETORY_LYSOSOME_TRAFFICKING: &[&str] = &["RAB27A", "RAB27B", "SYTL1", "SYTL2"];
const SECRETORY_FUSION_MACHINERY: &[&str] = &["VAMP7", "STX7", "SNAP23"];
const LYSOSOME_EXOCYTOSIS_REGULATORS: &[&str] = &["MCOLN1", "TFEB"];
const MHC_CLASS_II_CORE: &[&str] = &["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"];
const ANTIGEN_PROCESSING_PROTEASES: &[&str] = &["CTSS", "CTSB", "CTSD"];
const MHC_LOADING_ACCESSORY: &[&str] = &["CD74", "HLA-DM", "HLA-DOA"];
const LIPID_STORAGE_BUFFERING: &[&str] = &["DGAT1", "DGAT2", "PLIN2", "PLIN3"];
const LIPID_MOBILIZATION_UTILIZATION: &[&str] = &["PNPLA2", "LIPE", "CPT1A"];
const LIPID_LIPOPHAGY_CONTEXT: &[&str] = &["RAB7", "LAMP2", "SQSTM1"];
const LYSO_ESCRT_REPAIR: &[&str] = &["CHMP4B", "CHMP2A", "CHMP2B", "VPS4A", "VPS4B", "ALIX"];
const LYSO_CA_RECRUIT: &[&str] = &["MCOLN1", "SYT7"];
const LYSO_LYSOPHAGY_INIT: &[&str] = &["LGALS3", "SQSTM1", "TAX1BP1", "OPTN"];
const LYSO_MEM_STABILIZATION: &[&str] = &["LAMP1", "LAMP2", "TMEM106B"];
const LYSO_CA_RELEASE: &[&str] = &["MCOLN1", "TPCN1", "TPCN2"];
const CA_SIGNALING_ADAPTORS: &[&str] = &["CALM1", "CALM2", "CAMK2A", "CAMK2D"];
const MITO_CA_UPTAKE: &[&str] = &["MCU", "MICU1", "MICU2", "SMDT1"];
const CA_METABOLIC_ACTIVATION: &[&str] = &["PDHA1", "PDP1", "PDP2"];
const LYSO_ROS_SOURCES: &[&str] = &["CTSB", "CTSD", "LGMN"];
const IRON_REDOX_CONTEXT: &[&str] = &["FTH1", "FTL", "HMOX1"];
const ANTIOXIDANT_BUFFERING: &[&str] = &["PRDX1", "PRDX4", "TXN", "TXNRD1", "GSR"];
const ROS_DAMAGE_RESPONDERS: &[&str] = &["LGALS3", "SESN2"];
const LYSO_BIOGENESIS_PROGRAM: &[&str] = &[
    "TFEB", "TFE3", "MITF", "LAMP1", "LAMP2", "ATP6V1A", "ATP6V1B2", "CTSB", "CTSD",
];
const LYSO_ASSEMBLY_MATURATION: &[&str] = &["GNPTAB", "GNPTG", "SORT1", "VPS11", "VPS18"];
const REG_TFEB: &[&str] = &["TFEB", "TFE3", "MITF"];
const REG_MTOR: &[&str] = &["MTOR", "RPTOR", "MLST8", "RHEB", "AKT1"];
const REG_AMPK: &[&str] = &["PRKAA1", "PRKAA2", "PRKAB1", "PRKAG1"];
const REG_RAGULATOR: &[&str] = &["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5"];
const PROLIF: &[&str] = &["MKI67", "TOP2A", "PCNA", "MCM2", "MCM4"];
const APOP_PRO: &[&str] = &["BAX", "BAK1", "CASP3", "CASP7", "CASP8", "PMAIP1", "BBC3"];
const APOP_ANTI: &[&str] = &["BCL2", "BCL2L1"];

const GENE_SETS: [GeneSetDef; 56] = [
    GeneSetDef {
        name: "AUTOPHAGY_INIT",
        genes: AUTOPHAGY_INIT,
    },
    GeneSetDef {
        name: "AUTOPHAGY_ELONG",
        genes: AUTOPHAGY_ELONG,
    },
    GeneSetDef {
        name: "AUTOPHAGY_LATE",
        genes: AUTOPHAGY_LATE,
    },
    GeneSetDef {
        name: "AUTOPHAGY_CARGO",
        genes: AUTOPHAGY_CARGO,
    },
    GeneSetDef {
        name: "LYSO_VATP",
        genes: LYSO_VATP,
    },
    GeneSetDef {
        name: "LYSO_PROT",
        genes: LYSO_PROT,
    },
    GeneSetDef {
        name: "LYSO_MEM",
        genes: LYSO_MEM,
    },
    GeneSetDef {
        name: "LYSO_LMP",
        genes: LYSO_LMP,
    },
    GeneSetDef {
        name: "LYSO_STRESS",
        genes: LYSO_STRESS,
    },
    GeneSetDef {
        name: "AUTO_MITOPHAGY",
        genes: AUTO_MITOPHAGY,
    },
    GeneSetDef {
        name: "AUTO_AGGREPHAGY",
        genes: AUTO_AGGREPHAGY,
    },
    GeneSetDef {
        name: "AUTO_ERPHAGY",
        genes: AUTO_ERPHAGY,
    },
    GeneSetDef {
        name: "AUTO_FERRITINOPHAGY",
        genes: AUTO_FERRITINOPHAGY,
    },
    GeneSetDef {
        name: "AUTO_LIPOPHAGY",
        genes: AUTO_LIPOPHAGY,
    },
    GeneSetDef {
        name: "PERINUCLEAR_LYSOSOME",
        genes: PERINUCLEAR_LYSOSOME,
    },
    GeneSetDef {
        name: "PERIPHERAL_LYSOSOME",
        genes: PERIPHERAL_LYSOSOME,
    },
    GeneSetDef {
        name: "FERRITINOPHAGY",
        genes: FERRITINOPHAGY,
    },
    GeneSetDef {
        name: "FERROPTOSIS_DEFENSE",
        genes: FERROPTOSIS_DEFENSE,
    },
    GeneSetDef {
        name: "IRON_IMPORT_EXPORT",
        genes: IRON_IMPORT_EXPORT,
    },
    GeneSetDef {
        name: "LIPID_PEROXIDATION_CONTEXT",
        genes: LIPID_PEROXIDATION_CONTEXT,
    },
    GeneSetDef {
        name: "CHOLESTEROL_EXPORT",
        genes: CHOLESTEROL_EXPORT,
    },
    GeneSetDef {
        name: "CHOLESTEROL_IMPORT",
        genes: CHOLESTEROL_IMPORT,
    },
    GeneSetDef {
        name: "CHOLESTEROL_EFFLUX",
        genes: CHOLESTEROL_EFFLUX,
    },
    GeneSetDef {
        name: "ER_LYSO_CONTACT_TETHERING",
        genes: ER_LYSO_CONTACT_TETHERING,
    },
    GeneSetDef {
        name: "ER_LYSO_CALCIUM_AXIS",
        genes: ER_LYSO_CALCIUM_AXIS,
    },
    GeneSetDef {
        name: "ER_LYSO_STRESS_ADAPTORS",
        genes: ER_LYSO_STRESS_ADAPTORS,
    },
    GeneSetDef {
        name: "SECRETORY_LYSOSOME_TRAFFICKING",
        genes: SECRETORY_LYSOSOME_TRAFFICKING,
    },
    GeneSetDef {
        name: "SECRETORY_FUSION_MACHINERY",
        genes: SECRETORY_FUSION_MACHINERY,
    },
    GeneSetDef {
        name: "LYSOSOME_EXOCYTOSIS_REGULATORS",
        genes: LYSOSOME_EXOCYTOSIS_REGULATORS,
    },
    GeneSetDef {
        name: "MHC_CLASS_II_CORE",
        genes: MHC_CLASS_II_CORE,
    },
    GeneSetDef {
        name: "ANTIGEN_PROCESSING_PROTEASES",
        genes: ANTIGEN_PROCESSING_PROTEASES,
    },
    GeneSetDef {
        name: "MHC_LOADING_ACCESSORY",
        genes: MHC_LOADING_ACCESSORY,
    },
    GeneSetDef {
        name: "LIPID_STORAGE_BUFFERING",
        genes: LIPID_STORAGE_BUFFERING,
    },
    GeneSetDef {
        name: "LIPID_MOBILIZATION_UTILIZATION",
        genes: LIPID_MOBILIZATION_UTILIZATION,
    },
    GeneSetDef {
        name: "LIPID_LIPOPHAGY_CONTEXT",
        genes: LIPID_LIPOPHAGY_CONTEXT,
    },
    GeneSetDef {
        name: "LYSO_ESCRT_REPAIR",
        genes: LYSO_ESCRT_REPAIR,
    },
    GeneSetDef {
        name: "LYSO_CA_RECRUIT",
        genes: LYSO_CA_RECRUIT,
    },
    GeneSetDef {
        name: "LYSO_LYSOPHAGY_INIT",
        genes: LYSO_LYSOPHAGY_INIT,
    },
    GeneSetDef {
        name: "LYSO_MEM_STABILIZATION",
        genes: LYSO_MEM_STABILIZATION,
    },
    GeneSetDef {
        name: "LYSO_CA_RELEASE",
        genes: LYSO_CA_RELEASE,
    },
    GeneSetDef {
        name: "CA_SIGNALING_ADAPTORS",
        genes: CA_SIGNALING_ADAPTORS,
    },
    GeneSetDef {
        name: "MITO_CA_UPTAKE",
        genes: MITO_CA_UPTAKE,
    },
    GeneSetDef {
        name: "CA_METABOLIC_ACTIVATION",
        genes: CA_METABOLIC_ACTIVATION,
    },
    GeneSetDef {
        name: "LYSO_ROS_SOURCES",
        genes: LYSO_ROS_SOURCES,
    },
    GeneSetDef {
        name: "IRON_REDOX_CONTEXT",
        genes: IRON_REDOX_CONTEXT,
    },
    GeneSetDef {
        name: "ANTIOXIDANT_BUFFERING",
        genes: ANTIOXIDANT_BUFFERING,
    },
    GeneSetDef {
        name: "ROS_DAMAGE_RESPONDERS",
        genes: ROS_DAMAGE_RESPONDERS,
    },
    GeneSetDef {
        name: "LYSO_BIOGENESIS_PROGRAM",
        genes: LYSO_BIOGENESIS_PROGRAM,
    },
    GeneSetDef {
        name: "LYSO_ASSEMBLY_MATURATION",
        genes: LYSO_ASSEMBLY_MATURATION,
    },
    GeneSetDef {
        name: "REG_TFEB",
        genes: REG_TFEB,
    },
    GeneSetDef {
        name: "REG_MTOR",
        genes: REG_MTOR,
    },
    GeneSetDef {
        name: "REG_AMPK",
        genes: REG_AMPK,
    },
    GeneSetDef {
        name: "REG_RAGULATOR",
        genes: REG_RAGULATOR,
    },
    GeneSetDef {
        name: "PROLIF",
        genes: PROLIF,
    },
    GeneSetDef {
        name: "APOP_PRO",
        genes: APOP_PRO,
    },
    GeneSetDef {
        name: "APOP_ANTI",
        genes: APOP_ANTI,
    },
];

pub fn gene_sets() -> &'static [GeneSetDef] {
    &GENE_SETS
}
