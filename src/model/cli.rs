use std::path::PathBuf;

use crate::model::ctx::Mode;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Pipeline,
}

impl Default for RunMode {
    fn default() -> Self {
        Self::Standalone
    }
}

#[derive(Debug, Clone)]
pub struct CliArgs {
    pub input: PathBuf,
    pub cache_path: Option<PathBuf>,
    pub mode: Option<Mode>,
    pub manifest: Option<PathBuf>,
    pub timecourse: bool,
    pub run_mode: RunMode,
}
