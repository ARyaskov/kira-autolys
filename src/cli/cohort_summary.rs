use std::path::Path;

use crate::model::cohort_view::load_cohort_view;
use crate::report::cohort::{classes, metrics, signatures, therapy, vulnerabilities};
use crate::stage_error::StageError;

pub fn run_cohort_summary(input_dir: &Path, by_group: bool) -> Result<(), StageError> {
    let view = load_cohort_view(input_dir)?;
    let out_dir = input_dir.join("cohort");

    classes::write_classes(&view, &out_dir, by_group)?;
    metrics::write_metrics(&view, &out_dir)?;
    signatures::write_signatures(&view, &out_dir)?;
    therapy::write_therapy(&view, &out_dir)?;
    vulnerabilities::write_vulnerabilities(&view, &out_dir)?;

    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/cli/cohort_summary.rs"]
mod tests;
