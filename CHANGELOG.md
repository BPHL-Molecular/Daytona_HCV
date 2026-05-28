# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.0] - 2026-05-20

This release restructures the pipeline entry point into two independent SLURM
jobs, adds structured logging throughout, and hardens error handling in
`braken_phy.py`. The scientific pipeline (tools, parameters, references) is
unchanged — outputs are bit-for-bit equivalent to v1.1 for a successful run.

### Added

- **Two-stage SLURM workflow.**
  - `Daytona_HCV_part1_nextflow.sh` — runs the Nextflow stage and generates
    `./output/sum_report.txt`.
  - `Daytona_HCV_part2_phylogeny.sh` — runs `braken_phy.py` and the
    downstream tree-rendering / cleanup steps.
- **Merged pipeline log file.** Both parts write to the same SLURM log files
  (`hepatitis_pipeline.out` / `hepatitis_pipeline.err`); Part 1 uses
  `--open-mode=truncate`, Part 2 uses `--open-mode=append`.
- **Dedicated `braken_phy.py` log files** in `./output/logs/`:
  - `braken_phy_<TS>.log` — shell-captured stdout + stderr.
  - `braken_phy_<TS>.err` — shell-captured stderr only.
  - `braken_phy_py_<TS>.log` — Python-side structured log (`logging` module).
- **`BRAKEN_LOG_TS` environment variable**, exported by Part 2 and consumed
  by `braken_phy.py`, so the three log files above share the same timestamp.
- **Structured logging in `braken_phy.py`** via the standard `logging`
  module: timestamped messages with level (`INFO` / `WARNING` / `ERROR`),
  `INFO` to stdout, `WARNING+` to stderr, all levels to the Python log file.
- **Per-command logging** in `braken_phy.py`: every shell call is preceded
  by `RUN: <cmd>` and followed by `OK exit=0 (X.Xs)` or `FAILED exit=N`.
- **Per-sample input existence checks** in `braken_phy.py` (`kraken.out`,
  `fastq.gz`, reference FASTA) before launching the consensus chain.
- **End-of-run summary** in `braken_phf.py` reporting total samples,
  HCV-processed, no-HCV-skipped, parse errors, failures, and the names of
  processed and failed samples.
- **Sanity check at the start of Part 2** that aborts immediately if
  `./output` is missing (i.e., Part 1 was not run or did not produce
  output).
- `CHANGELOG.md` (this file).

### Changed

- **Pipeline entry point.** The single `Daytona_HCV.sh` script has been
  replaced by the two part scripts above. See the **Migration** section
  below.
- **`braken_phy.py` rewritten** for clarity and robustness:
  - All `os.system()` calls replaced with `subprocess.run(..., shell=True,
    check=False)` and explicit exit-code checks.
  - The three near-duplicate phylogeny branches (`--phy 1`, `2`, `3`) have
    been refactored into a single `run_phy_chain(level)` function.
  - Per-sample work is isolated in `process_sample(sample_id, tax_id, ...)`.
    A failure in one sample logs a full traceback and continues with the
    next sample, instead of silently corrupting the run.
  - Added type hints and module-level docstring.
- **Part 2 memory request reduced** from 300 GB to 100 GB. The phylogeny
  stage (mafft + snp-sites + iqtree + phytreeviz) is far less memory-hungry
  than the Nextflow stage and rarely needs more than a few GB.
- **SLURM banner** added at start and end of each part, including
  `SLURM_JOB_ID` and the parsed command-line argument.
- **`./work` cleanup and `./output → ./output-<dt>` rename** now happen at
  the end of Part 2 only (Part 1 must leave `./output` in place for Part 2
  to consume).

### Fixed

- **`braken_phy.py`: lost samples after the first non-HCV sample.** The
  original code's `else` branch used `break` when a sample had no detected
  HCV reads, which silently halted processing of **every subsequent sample**
  in `sum_report.txt`. Replaced with `continue` and a `WARNING` log entry.
  If you rely on the original behavior, change the marked `continue` near
  the end of the per-sample loop back to `break`.
- **`braken_phy.py`: silent tool failures.** BWA / samtools / ivar / mafft
  / snp-sites / iqtree / phytreeviz exit codes were never checked, so a
  failed tool produced an empty file and the next stage built on garbage.
  The pipeline now fails fast within a sample (and within the phylogeny
  chain) and reports which step failed with its exit code and elapsed time.
- **`braken_phy.py`: `\.` in the sed f-string.** Changed to a raw f-string
  (`rf"..."`) to avoid `SyntaxWarning: invalid escape sequence '\.'` on
  Python 3.12+.
- **`braken_phy.py`: phylogeny stage running on an empty consensus set.**
  Added an `n_hcv == 0` check that aborts before the `cat *.consensus.fa`
  step so the user gets a clear error instead of cryptic downstream failures.

### Removed

- `Daytona_HCV.sh` — superseded by the two part scripts. Keep a local copy
  if you need it as a reference; users should not invoke it any more.

### Migration from v1.1

Replace
```bash
sbatch Daytona_HCV.sh [genotype|subtype|both]
```
with either the manual two-step form
```bash
sbatch Daytona_HCV_part1_nextflow.sh
# wait for it to finish, then:
sbatch Daytona_HCV_part2_phylogeny.sh [genotype|subtype|both]
```
or chain the jobs via SLURM dependencies so Part 2 starts automatically
when Part 1 succeeds:
```bash
P1=$(sbatch --parsable Daytona_HCV_part1_nextflow.sh)
sbatch --dependency=afterok:$P1 Daytona_HCV_part2_phylogeny.sh both
```

If you previously edited `Daytona_HCV.sh` to set
`#SBATCH --mail-user=<your-email>`, apply the same edit to both
`Daytona_HCV_part1_nextflow.sh` and `Daytona_HCV_part2_phylogeny.sh`.

---

## [1.1] - 2025-12-30

- Initial public release with single `Daytona_HCV.sh` entry point.
- Genotype-level and subtype-level SNP-based phylogeny.
- Bootstrap test with 1,000 replicates.
- See repository tags for full release notes.
