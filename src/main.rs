use std::path::PathBuf;
use std::time::Instant;

use clap::{Parser, Subcommand, ValueEnum};
use tracing::{error, info, warn};
use tracing_subscriber::EnvFilter;

use kira_autolys::cli::cohort_summary::run_cohort_summary;
use kira_autolys::config::thresholds::Thresholds;
use kira_autolys::model::cli::{CliArgs, RunMode};
use kira_autolys::model::ctx::{Ctx, Mode};
use kira_autolys::{
    StageError, run_stage0, run_stage1, run_stage2, run_stage3, run_stage4, run_stage5, run_stage6,
    run_stage7, run_stage8, run_stage9, run_stage10, run_stage12, run_stage13,
    run_stage14_export_v2, run_stage16, run_stage17, run_stage18, run_stage19, run_stage20,
    run_stage21, run_stage22, run_stage23, run_stage24, run_stage25, run_stage26,
};

#[derive(Parser)]
#[command(name = "kira-autolys")]
#[command(version, about = "kira-autolys cohort utilities")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Run {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        cache: Option<PathBuf>,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum)]
        mode: Option<ModeOpt>,
        #[arg(long)]
        manifest: Option<PathBuf>,
        #[arg(long)]
        timecourse: bool,
        #[arg(long, value_enum)]
        schema: Option<SchemaOpt>,
        #[arg(long)]
        full: bool,
        #[arg(long)]
        no_scores_tsv: bool,
        #[arg(long)]
        fast: bool,
        #[arg(long, value_enum, default_value = "standalone")]
        run_mode: RunModeOpt,
    },
    CohortSummary {
        #[arg(long)]
        input: PathBuf,
        #[arg(long, value_enum)]
        by: Option<GroupBy>,
        #[arg(long, value_enum)]
        schema: Option<SchemaOpt>,
        #[arg(long)]
        full: bool,
    },
}

#[derive(Clone, ValueEnum)]
enum GroupBy {
    SampleGroup,
}

#[derive(Clone, ValueEnum)]
enum ModeOpt {
    Cell,
    Sample,
}

#[derive(Clone, ValueEnum)]
enum SchemaOpt {
    V2,
}

#[derive(Clone, ValueEnum)]
enum RunModeOpt {
    Standalone,
    Pipeline,
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .init();
    if let Err(err) = run() {
        error!(error = %err, "kira-autolys run failed");
        std::process::exit(1);
    }
}

fn run() -> Result<(), StageError> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Run {
            input,
            cache,
            out,
            mode,
            manifest,
            timecourse,
            schema,
            full,
            no_scores_tsv,
            fast,
            run_mode,
        } => {
            let mode = mode.map(|m| match m {
                ModeOpt::Cell => Mode::Cell,
                ModeOpt::Sample => Mode::Sample,
            });
            let run_mode = match run_mode {
                RunModeOpt::Standalone => RunMode::Standalone,
                RunModeOpt::Pipeline => RunMode::Pipeline,
            };
            let cli_args = CliArgs {
                input,
                cache_path: cache.clone(),
                mode,
                manifest,
                timecourse,
                run_mode,
            };
            let stage_out = if matches!(cli_args.run_mode, RunMode::Pipeline) {
                out.join("kira-autolys")
            } else {
                out.clone()
            };
            let mut ctx = Ctx::default();
            ctx.cache_path_override = cache;
            ctx.thresholds = Some(Thresholds::default());
            info!(simd_backend = %simd_backend_label(), "simd backend selected");

            run_stage_with_log("Stage 0 (input)", || run_stage0(&mut ctx, &cli_args))?;
            run_stage_with_log("Stage 1 (normalization)", || run_stage1(&mut ctx))?;
            run_stage_with_log("Stage 2 (autophagy)", || run_stage2(&mut ctx))?;
            run_stage_with_log("Stage 3 (lysosome)", || run_stage3(&mut ctx))?;
            run_stage_with_log("Stage 4 (regulatory)", || run_stage4(&mut ctx))?;
            run_stage_with_log("Stage 5 (survival)", || run_stage5(&mut ctx))?;
            if !fast {
                run_stage_with_log("Stage 7 (damage)", || run_stage7(&mut ctx))?;
                run_stage_with_log("Stage 8 (selectivity)", || run_stage8(&mut ctx))?;
                run_stage_with_log("Stage 9 (coupling)", || run_stage9(&mut ctx))?;
                run_stage_with_log("Stage 10 (cross-organelle)", || run_stage10(&mut ctx))?;
                run_stage_with_log("Stage 12 (classification)", || run_stage12(&mut ctx))?;
                run_stage_with_log("Stage 13 (vulnerabilities)", || run_stage13(&mut ctx))?;
                run_stage_with_log("Stage 16 (positioning)", || run_stage16(&mut ctx))?;
                run_stage_with_log("Stage 17 (ferroptosis)", || run_stage17(&mut ctx))?;
                run_stage_with_log("Stage 18 (cholesterol)", || run_stage18(&mut ctx))?;
                run_stage_with_log("Stage 19 (ER-lysosome)", || run_stage19(&mut ctx))?;
                run_stage_with_log("Stage 20 (secretory)", || run_stage20(&mut ctx))?;
                run_stage_with_log("Stage 21 (antigen)", || run_stage21(&mut ctx))?;
                run_stage_with_log("Stage 22 (lipid buffering)", || run_stage22(&mut ctx))?;
                run_stage_with_log("Stage 23 (lysosomal ROS)", || run_stage23(&mut ctx))?;
                run_stage_with_log("Stage 24 (calcium coupling)", || run_stage24(&mut ctx))?;
                run_stage_with_log("Stage 25 (biogenesis)", || run_stage25(&mut ctx))?;
                run_stage_with_log("Stage 26 (membrane repair)", || run_stage26(&mut ctx))?;
            }

            run_stage_with_log("Stage 6 (export v1)", || {
                run_stage6(&ctx, &stage_out, !no_scores_tsv)
            })?;
            if fast && (full || matches!(schema, Some(SchemaOpt::V2))) {
                return Err(StageError::Validation(
                    "--fast is incompatible with --schema v2/--full".to_string(),
                ));
            }
            if full || matches!(schema, Some(SchemaOpt::V2)) {
                run_stage_with_log("Stage 14 (export v2)", || {
                    run_stage14_export_v2(&ctx, &stage_out)
                })?;
            }
            Ok(())
        }
        Commands::CohortSummary {
            input,
            by,
            schema,
            full,
        } => {
            if let Some(schema_ref) = schema.as_ref() {
                match schema_ref {
                    SchemaOpt::V2 => {}
                }
            }
            if !full && schema.is_none() {
                return Err(StageError::Validation(
                    "--schema v2 or --full is required for cohort-summary".to_string(),
                ));
            }
            let by_group = matches!(by, Some(GroupBy::SampleGroup));
            run_cohort_summary(&input, by_group)
        }
    }
}

fn run_stage_with_log<F>(label: &str, mut f: F) -> Result<(), StageError>
where
    F: FnMut() -> Result<(), StageError>,
{
    info!(stage = label, "stage started");
    let started = Instant::now();
    let result = f();
    let elapsed_ms = started.elapsed().as_millis();
    match &result {
        Ok(()) => info!(
            stage = label,
            elapsed_ms = elapsed_ms as u64,
            "stage finished"
        ),
        Err(err) => warn!(
            stage = label,
            elapsed_ms = elapsed_ms as u64,
            error = %err,
            "stage failed"
        ),
    }
    result
}

fn simd_backend_label() -> &'static str {
    #[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
    {
        return "avx2";
    }
    #[cfg(all(feature = "simd", target_arch = "aarch64", target_feature = "neon"))]
    {
        return "neon";
    }
    #[cfg(not(any(
        all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"),
        all(feature = "simd", target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        "scalar"
    }
}

#[cfg(test)]
#[path = "../tests/src_inline/main_inline.rs"]
mod tests;
