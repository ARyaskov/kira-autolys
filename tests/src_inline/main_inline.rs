use super::*;

#[test]
fn run_mode_defaults_to_standalone() {
    let cli =
        Cli::try_parse_from(["kira-autolys", "run", "--input", "./in", "--out", "./out"]).unwrap();
    match cli.command {
        Commands::Run { run_mode, .. } => {
            assert!(matches!(run_mode, RunModeOpt::Standalone));
        }
        _ => panic!("expected run command"),
    }
}

#[test]
fn run_mode_pipeline_is_accepted() {
    let cli = Cli::try_parse_from([
        "kira-autolys",
        "run",
        "--input",
        "./in",
        "--out",
        "./out",
        "--run-mode",
        "pipeline",
    ])
    .unwrap();
    match cli.command {
        Commands::Run { run_mode, .. } => {
            assert!(matches!(run_mode, RunModeOpt::Pipeline));
        }
        _ => panic!("expected run command"),
    }
}
