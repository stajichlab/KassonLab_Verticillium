#!/usr/bin/env bash
# 01_matpredict.sh — predict the mating type (MAT1-1 / MAT1-2) of each Verticillium
# genome with MATPredict (https://github.com/stajichlab/MATPredict).
#
# Input : pipeline/mating_type/matpredict_samples.tsv (sample, taxid, genome path).
#         The genomes are the funannotate scaffolds in
#         results/annotation/*/predict_results/*.scaffolds.fa.
#         Taxids match pipeline/annotation/funannotate_meta.tsv
#         (27337 V. dahliae, 1051616 V. nonalfalfae).
# Output: results/mating_type/<sample>/ (detection_report.yaml, detected_loci.gff3, report.html)
#         results/mating_type/mating_type_summary.tsv (written by the summary step, below)
#
# Run from the nf/ directory (the paths above are relative to it).
# One array task per genome (4 genomes in the sheet):
#   sbatch --array=1-4 pipeline/mating_type/01_matpredict.sh
# Then the summary (login node is fine, it only reads YAML):
#   bash pipeline/mating_type/02_summarize.sh
#
# MATPredict runs the tblastn search on the genome, not on the funannotate
# proteins, so an unannotated or mis-modelled MAT gene is still found.
#
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 0-02:30:00
#SBATCH --job-name matpredict
#SBATCH -o logs/slurm/matpredict_%A_%a.out
#SBATCH -e logs/slurm/matpredict_%A_%a.err

set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?submit with sbatch --array=1-N}"

# Absolute paths. Do not use BASH_SOURCE: it fails in SLURM work dirs.
PROJ=/bigdata/stajichlab/shared/projects/KassonLab/Verticillium/nf
MATPREDICT=$HOME/projects/MATPredict
SAMPLESHEET=$PROJ/pipeline/mating_type/matpredict_samples.tsv
OUTROOT=$PROJ/results/mating_type

# MATPredict reads db/ and its NCBI cache from its own checkout.
export MATPREDICT_DB_ROOT=$MATPREDICT/db
export MATPREDICT_CACHE_DIR=$MATPREDICT/.matpredict_cache
export MATPREDICT_NCBI_EMAIL=${MATPREDICT_NCBI_EMAIL:-jason.stajich@ucr.edu}
export MATPREDICT_HTML=1

# Row N of the sheet, skipping the header.
read -r SAMPLE TAXID GENOME < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLESHEET")
: "${SAMPLE:?no row ${SLURM_ARRAY_TASK_ID} in $SAMPLESHEET}"

OUTDIR=$OUTROOT/$SAMPLE
mkdir -p "$OUTDIR"
echo "sample=$SAMPLE taxid=$TAXID genome=$GENOME out=$OUTDIR"

cd "$MATPREDICT"
pixi run matpredict detect \
    --genome "$PROJ/$GENOME" \
    --taxid "$TAXID" \
    --sample "$SAMPLE" \
    --out-dir "$OUTDIR"
