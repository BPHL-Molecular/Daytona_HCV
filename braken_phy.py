#!/usr/bin/env python3
"""
braken_phy.py
=============
HCV consensus + phylogeny driver. Reads ./output/sum_report.txt produced by
the Nextflow stage, extracts HCV reads per sample, builds consensus FASTAs,
then runs mafft -> snp-sites -> snp-dists -> iqtree -> phytreeviz at the
genotype level, subtype level, or both.

Logging
-------
This version replaces all bare ``os.system()`` calls and ``print()``
statements with a proper ``logging`` setup:

  * INFO/DEBUG -> stdout  (captured by SLURM stdout / shell tee)
  * WARNING/ERROR -> stderr  (captured by SLURM stderr / shell tee)
  * Everything (DEBUG+) -> ./output/logs/braken_phy_py_<TS>.log

Every shell command is run through :func:`run_cmd`, which:
  * Logs the command before it starts.
  * Times it.
  * Checks the exit code and logs success/failure with the elapsed time.
  * Raises :class:`RuntimeError` on non-zero exit (fail-fast within a sample
    or within the phylogeny chain).

Behaviour changes vs. the original
----------------------------------
1. ``os.system(...)`` -> ``subprocess.run(..., shell=True, check=False)`` with
   explicit exit-code checks. The original swallowed all errors.
2. The original used ``break`` when a sample had no HCV, which silently
   skipped every subsequent sample. This version uses ``continue`` and logs
   a WARNING. To restore the original behaviour, change the marked
   ``continue`` near the end of the per-sample loop back to ``break``.
3. If a single sample's consensus-building chain fails, this version logs
   the failure and moves on to the next sample (instead of aborting). The
   phylogeny chain still runs as long as at least one sample succeeded.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import time
from datetime import datetime


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
def setup_logger(log_dir: str) -> logging.Logger:
    """Create a logger with stdout + stderr + file handlers."""
    os.makedirs(log_dir, exist_ok=True)
    # Allow the wrapping shell script to pin the timestamp so the python
    # log file matches its sibling .log/.err files. Falls back to "now".
    ts = os.environ.get("BRAKEN_LOG_TS") or datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(log_dir, f"braken_phy_py_{ts}.log")

    logger = logging.getLogger("braken_phy")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)-7s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    class _InfoFilter(logging.Filter):
        def filter(self, record):
            return record.levelno < logging.WARNING

    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setLevel(logging.INFO)
    stdout_h.addFilter(_InfoFilter())
    stdout_h.setFormatter(fmt)

    stderr_h = logging.StreamHandler(sys.stderr)
    stderr_h.setLevel(logging.WARNING)
    stderr_h.setFormatter(fmt)

    file_h = logging.FileHandler(log_path, mode="w")
    file_h.setLevel(logging.DEBUG)
    file_h.setFormatter(fmt)

    logger.addHandler(stdout_h)
    logger.addHandler(stderr_h)
    logger.addHandler(file_h)

    logger.info(f"Python-side structured log: {log_path}")
    return logger


# ---------------------------------------------------------------------------
# Command runner
# ---------------------------------------------------------------------------
def run_cmd(cmd: str, logger: logging.Logger, step: str = "", fail_fast: bool = True):
    """Run a shell command, log it, check the exit code.

    On non-zero exit:
      * Logs at ERROR level with command, exit code, elapsed time.
      * Raises :class:`RuntimeError` if ``fail_fast`` is True.
    """
    label = f"[{step}] " if step else ""
    logger.info(f"{label}RUN: {cmd}")
    start = time.time()
    try:
        result = subprocess.run(cmd, shell=True, check=False)
    except Exception:
        logger.exception(f"{label}Exception while launching command")
        if fail_fast:
            raise
        return None
    elapsed = time.time() - start
    if result.returncode != 0:
        logger.error(
            f"{label}FAILED exit={result.returncode} after {elapsed:.1f}s :: {cmd}"
        )
        if fail_fast:
            raise RuntimeError(
                f"Command failed (exit {result.returncode}) at step {step!r}: {cmd}"
            )
    else:
        logger.info(f"{label}OK exit=0 ({elapsed:.1f}s)")
    return result


# ---------------------------------------------------------------------------
# Per-sample consensus chain
# ---------------------------------------------------------------------------
def process_sample(sample_id: str, tax_id: str, current_dir: str,
                   logger: logging.Logger) -> None:
    """Extract HCV reads, align, dedup, and build consensus for one sample."""
    logger.info(f"--- Sample {sample_id} (taxID={tax_id}) ---")
    out = f"{current_dir}/output/extract"
    ref = f"{current_dir}/reference/hcv/KC248195.fasta"
    fq1 = f"{current_dir}/fastqs/hcv/{sample_id}_1.fastq.gz"
    fq2 = f"{current_dir}/fastqs/hcv/{sample_id}_2.fastq.gz"
    kraken = f"{current_dir}/output/{sample_id}/kraken_out/{sample_id}_kraken.out"
    prefix = f"{out}/{sample_id}_{tax_id}"

    # Sanity check that the inputs the original assumes actually exist.
    for required in (kraken, fq1, fq2, ref):
        if not os.path.exists(required):
            logger.error(f"Missing input for {sample_id}: {required}")
            raise FileNotFoundError(required)

    commands = [
        ("extract_kraken_reads",
         f"python ./extract_kraken_reads.py "
         f"-k {kraken} -s {fq1} -s2 {fq2} "
         f"-o {prefix}_1.fq -o2 {prefix}_2.fq -t {tax_id}"),
        ("bwa_mem",
         f"singularity exec docker://staphb/bwa:0.7.17 bwa mem "
         f"{ref} {prefix}_1.fq {prefix}_2.fq > {prefix}_aln.sam"),
        ("samtools_view",
         f"singularity exec docker://staphb/samtools:1.12 samtools view "
         f"-F 4 -u -h -bo {prefix}_aln.bam {prefix}_aln.sam"),
        ("samtools_namesort",
         f"singularity exec docker://staphb/samtools:1.12 samtools sort "
         f"-n -o {prefix}.namesorted.bam {prefix}_aln.bam"),
        ("samtools_fixmate",
         f"singularity exec docker://staphb/samtools:1.12 samtools fixmate "
         f"-m {prefix}.namesorted.bam {prefix}.fixmate.bam"),
        ("samtools_possort",
         f"singularity exec docker://staphb/samtools:1.12 samtools sort "
         f"-o {prefix}.positionsort.bam {prefix}.fixmate.bam"),
        ("samtools_markdup",
         f"singularity exec docker://staphb/samtools:1.12 samtools markdup "
         f"-r {prefix}.positionsort.bam {prefix}.dedup.bam"),
        ("samtools_finalsort",
         f"singularity exec docker://staphb/samtools:1.12 samtools sort "
         f"-o {prefix}.sorted.bam {prefix}.dedup.bam"),
        ("ivar_consensus",
         f"singularity exec docker://staphb/samtools:1.12 samtools mpileup "
         f"-A -B -d 8000 --reference {ref} -Q 0 {prefix}.sorted.bam | "
         f"singularity exec docker://staphb/ivar:latest ivar consensus "
         f"-t 0 -m 10 -n N -p {prefix}.consensus"),
    ]
    for step, cmd in commands:
        run_cmd(cmd, logger, step=f"{sample_id}/{step}")
    logger.info(f"--- Sample {sample_id} done ---")


# ---------------------------------------------------------------------------
# Phylogeny chain
# ---------------------------------------------------------------------------
def run_phy_chain(level: str, current_dir: str, logger: logging.Logger) -> None:
    """Run mafft -> snp-sites -> snp-dists -> iqtree -> phytreeviz at one level."""
    ref_map = {
        "genotype": "hcv_7genotypes_aln",
        "subtype":  "hcv_allsubtypes_aln",
    }
    ref = ref_map[level]
    out = f"{current_dir}/output/extract"
    logger.info(f"===== Phylogeny chain: {level} (ref={ref}) =====")
    commands = [
        (f"mafft_{level}",
         f"singularity exec docker://staphb/mafft:latest mafft "
         f"--addfragments {out}/sum_consensus.fa "
         f"{current_dir}/reference/HCV_align/{ref} > {out}/mafft_{level}_msa"),
        (f"snpsites_{level}",
         f"singularity exec docker://staphb/snp-sites:2.3.3 snp-sites "
         f"-o {out}/SNPs_{level}.fasta -c {out}/mafft_{level}_msa"),
        (f"snpdists_{level}",
         f"singularity exec docker://staphb/snp-dists:0.6.2 snp-dists "
         f"{out}/mafft_{level}_msa > {out}/pairwise_matrix_{level}.tsv"),
        (f"iqtree_{level}",
         f"singularity exec docker://staphb/iqtree:1.6.7 iqtree "
         f"-s {out}/SNPs_{level}.fasta -m MFP+ASC -nt AUTO "
         f"-bb 1000 -alrt 1000 -pre {out}/SNPs_boot_{level} -ntmax 9"),
        (f"phytreeviz_{level}_pdf",
         f"singularity exec docker://staphb/phytreeviz:latest phytreeviz "
         f"-i {out}/SNPs_boot_{level}.contree "
         f"-o {out}/SNPs_boot_{level}.contree.pdf --show_confidence"),
        (f"phytreeviz_{level}_png",
         f"singularity exec docker://staphb/phytreeviz:latest phytreeviz "
         f"-i {out}/SNPs_boot_{level}.contree "
         f"-o {out}/SNPs_boot_{level}.contree.png --show_confidence"),
    ]
    for step, cmd in commands:
        run_cmd(cmd, logger, step=step)
    logger.info(f"===== Phylogeny chain ({level}) finished =====")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description="Parameters for braken_phy.py.")
    parser.add_argument(
        "-p", "--phy", type=int, default=1,
        help="which phylogeny level is generated? Input 1, 2 or 3. "
             "1 means using 7 HCV genotypes as references. "
             "2 means using 61 HCV subtypes as references. "
             "3 means both above.",
    )
    args = parser.parse_args()

    current_dir = os.getcwd()
    log_dir = os.path.join(current_dir, "output", "logs")
    logger = setup_logger(log_dir)

    logger.info("=" * 67)
    logger.info("braken_phy.py started")
    logger.info(f"  working dir : {current_dir}")
    logger.info(f"  --phy       : {args.phy}")
    logger.info(f"  PID         : {os.getpid()}")
    logger.info(f"  SLURM_JOB_ID: {os.environ.get('SLURM_JOB_ID', 'N/A')}")
    logger.info("=" * 67)

    if args.phy not in (1, 2, 3):
        logger.error(
            f"Braken_phy.py Parameter Error! --phy={args.phy} is invalid. "
            "Only int 1 or 2 or 3 is accepted."
        )
        return 2

    sum_report = os.path.join(current_dir, "output", "sum_report.txt")
    if not os.path.isfile(sum_report):
        logger.error(f"Input not found: {sum_report}")
        logger.error("Did you run the Nextflow stage (Part 1) first?")
        return 1

    try:
        with open(sum_report, "r") as kreport:
            lines = kreport.readlines()
    except Exception:
        logger.exception(f"Failed to read {sum_report}")
        return 1

    logger.info(f"Loaded {len(lines)} lines from sum_report.txt (incl. header)")

    n_total = 0
    n_hcv = 0
    n_skip_no_hcv = 0
    n_parse_err = 0
    failed_samples: list[str] = []
    processed_samples: list[str] = []

    for raw in lines[1:]:
        n_total += 1
        try:
            l_parse = raw.strip().split("\t")
            sample_id = l_parse[0].strip()
            species_group = l_parse[1].strip().split(",")
            species_items = species_group[1].strip().split("|")
            tax = species_items[0].strip()
            tax_id = species_items[1].strip()
        except (IndexError, ValueError) as e:
            logger.error(f"Could not parse sum_report line: {raw!r} ({e})")
            n_parse_err += 1
            continue

        substring = "Hepatitis C"
        if substring.lower() in tax.lower():
            try:
                process_sample(sample_id, tax_id, current_dir, logger)
                n_hcv += 1
                processed_samples.append(sample_id)
            except Exception:
                # Log full traceback; do NOT abort the whole pipeline -
                # other samples can still produce useful trees.
                logger.exception(
                    f"Sample {sample_id} failed during consensus building; "
                    "continuing with next sample."
                )
                failed_samples.append(sample_id)
                continue
        else:
            logger.warning(
                f"No Hepatitis C virus genotype found in {sample_id} "
                f"(top species: {tax!r}). Please check this sample."
            )
            n_skip_no_hcv += 1
            # NOTE: original code used `break` here, which would halt the
            # whole loop the moment one sample lacked HCV. That looks like
            # a bug - it silently drops all subsequent samples. We use
            # `continue` instead. To restore original behaviour, replace
            # the next line with `break`.
            continue

    logger.info("-" * 67)
    logger.info(
        f"Per-sample summary: total={n_total}, HCV_processed={n_hcv}, "
        f"no_HCV_skipped={n_skip_no_hcv}, parse_errors={n_parse_err}, "
        f"failures={len(failed_samples)}"
    )
    if processed_samples:
        logger.info(f"Processed samples ({len(processed_samples)}): "
                    f"{', '.join(processed_samples)}")
    if failed_samples:
        logger.error(f"Failed samples ({len(failed_samples)}): "
                     f"{', '.join(failed_samples)}")

    if n_hcv == 0:
        logger.error("No HCV samples were processed successfully. "
                     "Aborting before phylogeny stage.")
        return 1

    # ---- Combine per-sample consensuses --------------------------------
    logger.info("Combining per-sample consensus FASTAs into sum_consensus.fa")
    run_cmd(
        f"cat {current_dir}/output/extract/*.consensus.fa "
        f"> {current_dir}/output/extract/sum_consensus.fa",
        logger, step="cat_consensus",
    )
    # Raw string to avoid SyntaxWarning on \. in newer Python.
    run_cmd(
        rf"sed -i 's/>Consensus_/>/g; s/\.consensus_threshold_.*//g' "
        rf"{current_dir}/output/extract/sum_consensus.fa",
        logger, step="sed_consensus_headers",
    )

    # ---- Phylogeny -----------------------------------------------------
    try:
        if args.phy == 1:
            run_phy_chain("genotype", current_dir, logger)
        elif args.phy == 2:
            run_phy_chain("subtype", current_dir, logger)
        elif args.phy == 3:
            run_phy_chain("genotype", current_dir, logger)
            run_phy_chain("subtype", current_dir, logger)
    except Exception:
        logger.exception("Phylogeny chain failed")
        return 1

    logger.info("=" * 67)
    logger.info("braken_phy.py finished successfully")
    logger.info("=" * 67)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception:
        # Last-resort guard: even if our logger isn't set up, dump traceback.
        logging.basicConfig(level=logging.ERROR)
        logging.getLogger("braken_phy").exception(
            "FATAL: unhandled exception in braken_phy.py"
        )
        sys.exit(1)
