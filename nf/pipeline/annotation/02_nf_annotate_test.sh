#!/usr/bin/env bash
# 02_nf_annotate_test.sh — test the annotate step of nf_funannotate1 on ONE genome
# (MTK05, V. dahliae MD4) using the funannotate-1.9.0.sif image.
# Reads the existing prediction through a symlink at results/annotation_test/<id>/predict_results,
# so results/annotation is not touched. Own launch dir/work; no -resume carry-over.
#   PIPELINE=/path/to/checkout sbatch pipeline/annotation/02_nf_annotate_test.sh [extra args]
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 32G
#SBATCH -t 3-00:00:00
#SBATCH --job-name nf_annot_test
#SBATCH -o logs/slurm/nf_annot_test_%j.out
#SBATCH -e logs/slurm/nf_annot_test_%j.err
set -euo pipefail
PROJECT="${PROJECT:-$PWD}"
PIPELINE="${PIPELINE:-/bigdata/stajichlab/jstajich/projects/nf/nf_funannotate1}"
LAUNCH_DIR="$PROJECT/.nf_launch/annotate_test"
source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow apptainer 2>/dev/null || module load nextflow
mkdir -p "$LAUNCH_DIR/logs/nextflow"
cd "$LAUNCH_DIR"
nextflow run "$PIPELINE" \
    -profile annotate,slurm,ucr_hpcc,singularity \
    -c "$PIPELINE/conf/site_ucr_hpcc.config" \
    -c "$PIPELINE/conf/site_ucr_hpcc_singularity.config" \
    -resume \
    --samples "$LAUNCH_DIR/samples_test.csv" \
    --target  "$PROJECT/results/annotation_test" \
    --run_sra_fetch false \
    --run_annotate true --run_antismash true --run_interpro true \
    --proteins /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/swissprot_fungi.faa \
    --sbt_template /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/template.sbt \
    --sbt_dir "$PROJECT/lib/sbt" \
    "$@"
