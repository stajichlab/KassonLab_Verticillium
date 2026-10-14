#!/usr/bin/env bash
# 01_nf_funannotate.sh — sbatch launcher for stajichlab/nf_funannotate1 (HPCC).
#
# Drives genome cleaning -> masking -> training -> gene prediction -> functional
# annotation over samples_funannotate.csv (see 00_make_samples_funannotate.sh,
# which regenerates that sheet from samples.csv + the current genome FASTAs).
#
# Nextflow keeps per-launchDir state (work/, .nextflow/, .nextflow.log) that
# is NOT safe to share between concurrently-running pipelines: the project
# root is the launchDir for the still-running nf_aaftf pipeline (its own
# workDir is work/aaftf, isolated by profile_aaftf.config, but .nextflow/ and
# .nextflow.log are launchDir-wide and get rotated/locked by whichever
# `nextflow run` starts most recently). To avoid any interference with the
# live AAFTF run, this launches from its own subdirectory (.nf_launch/annotate)
# with its own work/ and .nextflow/ state, entirely separate from the
# project-root launchDir the AAFTF pipeline uses.
#
# Submit from this project dir. MASK_MODE is REQUIRED -- no default -- so a
# bare invocation fails loudly instead of silently falling back to tantan:
#   MASK_MODE=tantan   sbatch pipeline/annotation/01_nf_funannotate.sh
#   MASK_MODE=earlgrey sbatch pipeline/annotation/01_nf_funannotate.sh [extra nextflow args...]
#
# EarlGrey TE-discovery masking instead of tantan (MASK_MODE=earlgrey, above)
# uses params_earlgrey.json (verified 2026-09-26 via a real stub-run that
# -params-file genuinely switches the workflow graph to
# MASK_GENOME:EARLGREY_BUILD_LIB, not just a no-op param) plus the `earlgrey`
# profile axis for its resource tiers. EarlGrey is a multi-day-per-species
# process (7-day walltime cap on its label) -- much heavier than tantan.
#
# Gene prediction (TRAIN_PREDICT:FUNANNOTATE_PREDICT) runs automatically once
# masking is done -- no extra flag needed. Functional annotation does NOT:
# every ANNOTATE_GENOME step (funannotate annotate, antiSMASH, InterProScan,
# SignalP, DeepTMHMM, funannotate update) is individually gated by its own
# --run_* flag (run_annotate/run_antismash/run_interpro/run_signalp/
# run_deeptmhmm/run_update), all false by default in
# conf/profile_annotate.config -- pass the ones you want, e.g.:
#   MASK_MODE=earlgrey sbatch pipeline/annotation/01_nf_funannotate.sh --run_annotate true
#
# Quick graph validation (stub, runs in seconds on a login node):
#   nextflow run $HOME/projects/nf/nf_funannotate1 -profile test -stub-run
#
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name nf_funannotate
#SBATCH -o logs/slurm/nf_funannotate_%j.out
#SBATCH -e logs/slurm/nf_funannotate_%j.err

set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
PIPELINE="${PIPELINE:-stajichlab/nf_funannotate1}"  # GitHub repo; Nextflow fetches and caches it in ~/.nextflow/assets
REVISION="${REVISION:-}"
LAUNCH_DIR="${LAUNCH_DIR:-$PROJECT/.nf_launch/annotate}"
SAMPLES="${SAMPLES:-$LAUNCH_DIR/samples_funannotate.csv}"
# MASK_MODE has no default: masking silently reverting from earlgrey to tantan
# (or vice versa) between runs produces a DIFFERENT masked genome under
# storeDir-caching, so TRAIN_PREDICT/BUSCO_COMPLETENESS treat it as
# unrelated work rather than resuming -- it must be a deliberate choice on
# every invocation, not an accidental fallback.
#   MASK_MODE=tantan    -> PROFILES=annotate,slurm,ucr_hpcc,singularity
#   MASK_MODE=earlgrey  -> PROFILES=annotate,earlgrey,slurm,ucr_hpcc,singularity
#                          PARAMS_FILE=pipeline/annotation/params_earlgrey.json
# Either can still be overridden explicitly via PROFILES=/PARAMS_FILE=.
MASK_MODE="${MASK_MODE:?Set MASK_MODE=tantan or MASK_MODE=earlgrey (no default). Example: MASK_MODE=earlgrey sbatch pipeline/annotation/01_nf_funannotate.sh}"
case "$MASK_MODE" in
    tantan)
        PROFILES="${PROFILES:-annotate,slurm,ucr_hpcc,singularity}"
        PARAMS_FILE="${PARAMS_FILE:-}"
        ;;
    earlgrey)
        PROFILES="${PROFILES:-annotate,earlgrey,slurm,ucr_hpcc,singularity}"
        PARAMS_FILE="${PARAMS_FILE:-pipeline/annotation/params_earlgrey.json}"
        ;;
    *)
        echo "MASK_MODE must be 'tantan' or 'earlgrey', got: '$MASK_MODE'" >&2
        exit 1
        ;;
esac
# Resolve relative to PROJECT (the caller's cwd), not LAUNCH_DIR -- this
# script cd's into LAUNCH_DIR below before invoking nextflow.
case "$PARAMS_FILE" in
    /*|"") ;;
    *) PARAMS_FILE="$PROJECT/$PARAMS_FILE" ;;
esac

# ucr_hpcc alone (native Lmod modules) is NOT sufficient for the 'genome_clean'
# label: its `AAFTF` module only sets AAFTF_DB + pulls in singularity as a
# dependency -- it never actually puts the AAFTF/taxonkit CLIs on PATH
# (confirmed via `module show AAFTF`, 2026-09-26). The real AAFTF+taxonkit+
# FCS-GX toolchain only exists bundled in conf/provision_singularity.config's
# container_genome_clean (AAFTF.sif). Layer `singularity` on top of
# `ucr_hpcc` (not instead of it) so that container overrides just the labels
# that need it (genome_clean, ...) while everything else keeps ucr_hpcc's
# native-module beforeScripts/SLURM safety net.

mkdir -p "$PROJECT/logs/slurm" "$LAUNCH_DIR/logs/nextflow"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow apptainer 2>/dev/null || module load nextflow

# Regenerate the funannotate samplesheet from samples.csv every run, so it
# always reflects the latest AAFTF genome output (see 00_make_samples_funannotate.sh
# for where GENOME currently resolves to and why that will need updating).
# GENOME paths in the sheet are absolute, so it's safe to launch from
# LAUNCH_DIR instead of PROJECT.
"$PROJECT/pipeline/annotation/00_make_samples_funannotate.sh" "$PROJECT/samples.csv" "$SAMPLES"

cd "$LAUNCH_DIR"

nextflow run "${PIPELINE}" ${REVISION:+-r "${REVISION}"} \
    -profile "${PROFILES}" \
    -resume \
    ${PARAMS_FILE:+-params-file "$PARAMS_FILE"} \
    --samples "$SAMPLES" \
    --target  "$PROJECT/results/annotation" \
    --run_sra_fetch true \
    --proteins /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/swissprot_fungi.faa \
    --run_annotate true --run_antismash true --run_interpro true \
    --sbt_template /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/template.sbt \
    --sbt_dir "$PROJECT/lib/sbt" \
    "$@"
