# Pipeline Overview

kira-autolys is a deterministic, stage‑based pipeline. Each stage reads from `Ctx`, computes defined metrics, and writes back without altering previous outputs.

## Run modes

`kira-autolys run` supports:

* `--run-mode standalone` (default): direct input parsing (`matrix.mtx(.gz)`, `features.tsv/genes.tsv(.gz)`, `barcodes.tsv(.gz)`, or `.h5ad`).
* `--run-mode pipeline`: shared-cache mode for MTX datasets.

Pipeline mode behavior:

* Detect dataset prefix in `--input` directory (non-recursive) using file names like `PREFIX_matrix.mtx(.gz)`, `PREFIX_features.tsv(.gz)`, `PREFIX_barcodes.tsv(.gz)`.
* Resolve cache filename deterministically:
* no prefix: `kira-organelle.bin`
* with prefix: `PREFIX.kira-organelle.bin`
* Open cache read-only with mmap and run downstream normalization/access from cache-backed CSC arrays.
* If cache file is missing, MTX discovery/parsing is validated and then execution hard-fails with a clear pipeline-mode error instructing to provide shared cache.
* Output root for this stage is `<out>/kira-autolys/` with required machine artifacts:
* `autolys.tsv` (per-cell contract table)
* `summary.json` (run-level aggregates)
* `panels_report.tsv` (panel audit)
* `pipeline_step.json` (aggregator manifest)

Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md).

Stages:

1. Stage 0: input detection and gene panel resolution
2. Stage 1: expression loading and normalization
3. Stage 2: autophagy aggregation (AFP components)
4. Stage 3: lysosomal dependency (LDS)
5. Stage 4: regulatory axis (TFEB/mTOR)
6. Stage 5: stress survival mode integration (SSM)
7. Stage 6: export schema v1
8. Stage 7: lysosomal damage (LDI)
9. Stage 8: autophagy selectivity fingerprint
10. Stage 9: AMPK–mTOR–lysosome coupling (LSI)
11. Stage 10: cross‑organelle integration (ERDI and mito coupling)
12. Stage 11: therapy delta analysis
13. Stage 12: dependency class assignment
14. Stage 13: vulnerability hypotheses
15. Stage 14: export schema v2
16. Stage 16–22: specialized lysosomal axes (positioning, ferroptosis, cholesterol, ER contacts, secretion, immune load, lipid buffering)
17. Stage 23: lysosomal ROS handling
18. Stage 24: lysosome–mitochondria calcium coupling
19. Stage 25: lysosome biogenesis pressure
20. Stage 26: lysosomal membrane repair capacity

Stages are designed to be modular, deterministic, and auditable.

## Run modes

`kira-autolys run` supports two deterministic execution modes:

* `--run-mode standalone` (default): reads MTX/H5AD directly.
* `--run-mode pipeline`: resolves dataset prefix in `--input` and reads shared cache file via mmap:
  * `kira-organelle.bin` (standard names)
  * `<PREFIX>.kira-organelle.bin` (prefixed names)

Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md). If the cache file is missing, kira-autolys validates MTX inputs and then returns an explicit error instructing pipeline usage with the shared cache.
