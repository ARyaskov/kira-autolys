# kira-autolys

Deterministic, explainable autophagy/lysosome dependency QC for single-cell expression data.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-autolys
```

## Usage examples

Standalone run (single input, per-cell outputs):

```bash
kira-autolys run \
  --input ./data/pbmc3k \
  --out ./out/pbmc3k \
  --mode cell
```

Standalone run (sample mode):

```bash
kira-autolys run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample
```

Pipeline run (shared cache lookup + pipeline artifacts):

```bash
kira-autolys run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample \
  --run-mode pipeline
```

Full/v2 export:

```bash
kira-autolys run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample \
  --schema v2
```

Cohort summary from v2 artifacts:

```bash
kira-autolys cohort-summary \
  --input ./out/inf/kira-autolys \
  --schema v2
```

## Modes

- `--run-mode standalone` (default): standalone behavior and outputs.
- `--run-mode pipeline`: pipeline contract mode for `kira-organelle`.

## Pipeline cache lookup

In pipeline mode, `kira-autolys` searches the input directory for shared cache as specified in [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md):

- no prefix: `kira-organelle.bin`
- prefixed dataset: `<PREFIX>.kira-organelle.bin`

Prefix is detected non-recursively from names like `<PREFIX>_matrix.mtx(.gz)`, `<PREFIX>_features.tsv(.gz)`, `<PREFIX>_barcodes.tsv(.gz)`.

Behavior:

- cache exists and valid: use shared cache path.
- cache missing in pipeline + MTX directory input: hard error (cache required).
- cache missing in other pipeline inputs: warn and fall back to direct input.
- cache exists but invalid: hard error (no silent fallback).

## Pipeline output contract

In pipeline mode, outputs are written to:

- `--out <DIR>` -> `<DIR>/kira-autolys/`

Required artifacts:

- `autolys.tsv` (per-observation contract table)
- `summary.json` (run-level aggregates)
- `panels_report.tsv` (panel audit)
- `pipeline_step.json` (ingestion manifest for `kira-organelle`)

Additional artifacts may be emitted (`scores.tsv`, `report.txt`, `report_summary.json`).
All TSV float values are fixed `%.6f`.

## Shared cache specification

- Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)
- Reader validates header/magic/version/endian/header-size/file-bytes, header CRC64-ECMA, section bounds, string tables, and CSC invariants.
