#!/usr/bin/env bash
# 04_redo_sqn.sh — re-run funannotate's tbl2asn step with the current lib/sbt/<BioSample>.sbt
# and write a new .sqn, without re-running funannotate annotate.
# One array task per row of samples_funannotate.csv (rows 2..N):
#   sbatch --array=1-4 pipeline/annotation/04_redo_sqn.sh
# Reads  results/annotation/<Species_Strain>/annotate_misc/tbl2asn/genome.{tbl,fsa}
# Writes results/annotation/<Species_Strain>/annotate_results/<Species_Strain>.sqn
#   (old file kept as .sqn.bak the first time). The "Annotated using 1.9.0" comment is
#   rewritten to "Annotated using funannotate-1.9.0" in the .sqn and the .gbk.
# Then regenerates, from the same genome.tbl and the current <Species_Strain>.scaffolds.fa
#   (see regen_funannotate_outputs.py): .gbk (from tbl2asn genome.gbf), .tbl, .gff3,
#   .proteins.fa, .mrna-transcripts.fa, .cds-transcripts.fa, .agp, .contigs.fsa, .stats.json.
#   Old versions are kept as <file>.bak the first time. .scaffolds.fa is not rewritten.
#   Run 05_remove_organelle_contigs.sh first, so genome.tbl and .scaffolds.fa agree.
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -t 04:00:00
#SBATCH --job-name redo_sqn
#SBATCH -o logs/slurm/redo_sqn_%A_%a.out
#SBATCH -e logs/slurm/redo_sqn_%A_%a.err
set -euo pipefail
PROJECT="${PROJECT:-$PWD}"
SIF="${SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/funannotate-1.9.0.sif}"
OLD_TEXT="Annotated using 1.9.0"
NEW_TEXT="Annotated using funannotate-1.9.0"
SCRATCH="${SCRATCH:?SCRATCH is not set}"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || true

# Row N+1 of samples_funannotate.csv: SPECIES,STRAIN,ASMID,...
row=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PROJECT/samples_funannotate.csv")
species=$(cut -d, -f1 <<<"$row")
strain=$(cut -d, -f2 <<<"$row")
asmid=$(cut -d, -f3 <<<"$row")
transl=$(cut -d, -f6 <<<"$row")
id="${species// /_}_${strain}"
biosample=$(awk -F, -v a="$asmid" '$1==a {print $10}' "$PROJECT/samples.csv")
sbt="$PROJECT/lib/sbt/${biosample}.sbt"
base="$PROJECT/results/annotation/$id"
res="$base/annotate_results"
echo "task=$SLURM_ARRAY_TASK_ID id=$id asm=$asmid biosample=$biosample sbt=$sbt"
[ -n "$biosample" ] && [ -f "$sbt" ] || { echo "ERROR: SBT not found: $sbt" >&2; exit 1; }
for f in genome.tbl genome.fsa; do
    [ -f "$base/annotate_misc/tbl2asn/$f" ] || { echo "ERROR: missing $base/annotate_misc/tbl2asn/$f" >&2; exit 1; }
done

# Work on a copy in node-local scratch, so the funannotate output is untouched on failure.
work="$SCRATCH/redo_sqn_$id"
rm -rf "$work"; mkdir -p "$work"
cp "$base/annotate_misc/tbl2asn/genome.tbl" "$base/annotate_misc/tbl2asn/genome.fsa" "$work/"
cp "$sbt" "$work/"
cd "$work"

# Remove CDS models with an in-frame stop (tbl2asn SEQ_FEAT.InternalStop) before tbl2asn.
# A splice site is moved only where PASA/transcript introns support the new site;
# otherwise the gene is dropped. Report: internal_stops.tsv (installed at the end).
FIXSTOPS="$PROJECT/pipeline/annotation/06_fix_internal_stops.py"
[ -f "$FIXSTOPS" ] || { echo "ERROR: missing $FIXSTOPS" >&2; exit 1; }
evidence=()
for f in "$base/training/funannotate_train.pasa.gff3" "$base/training/funannotate_train.transcripts.gff3"; do
    [ -f "$f" ] && evidence+=("$f")
done
/usr/bin/python3.12 "$FIXSTOPS" --tbl "$work/genome.tbl" --fasta "$work/genome.fsa" \
    --table "${transl:-1}" --introns "${evidence[@]}" \
    --out-tbl "$work/genome.fixed.tbl" --report "$work/internal_stops.tsv"
mv "$work/genome.fixed.tbl" "$work/genome.tbl"

# Same arguments funannotate annotate logged (logfiles/funannotate-annotate.log).
apptainer exec -B /bigdata --env PATH=/pixi/.pixi/envs/base/bin:/usr/bin:/bin "$SIF" \
    /pixi/.pixi/envs/base/bin/python3.8 \
    /pixi/.pixi/envs/base/lib/python3.8/site-packages/funannotate/aux_scripts/tbl2asn_parallel.py \
    -i "$work/genome.tbl" -f "$work/genome.fsa" -o "$work" \
    --sbt "$work/$(basename "$sbt")" -d "$work/discrepency.report.txt" \
    -s "$species" -t "-l paired-ends" -v 1 -c "${SLURM_CPUS_PER_TASK:-8}" --strain "$strain"

[ -s genome.sqn ] || { echo "ERROR: tbl2asn produced no genome.sqn" >&2; exit 1; }

# Gate: no internal stops left in the tbl, and tbl2asn (genome.val) flags none we missed.
/usr/bin/python3.12 "$FIXSTOPS" --tbl "$work/genome.tbl" --fasta "$work/genome.fsa" \
    --table "${transl:-1}" --val "$work/genome.val" --check \
    || { echo "ERROR: internal stop codons remain, see $work/genome.val" >&2; exit 1; }
sed -i "s/$OLD_TEXT/$NEW_TEXT/" genome.sqn

[ -s genome.gbf ] || { echo "ERROR: tbl2asn produced no genome.gbf" >&2; exit 1; }
sed -i "s/$OLD_TEXT/$NEW_TEXT/" genome.gbf

# Regenerate gff3/proteins/transcripts/agp/contigs/stats in the work dir first, so a failure
# leaves annotate_results untouched.
mkdir -p "$work/regen"
apptainer exec -B /bigdata --env PATH=/pixi/.pixi/envs/base/bin:/usr/bin:/bin "$SIF" \
    /pixi/.pixi/envs/base/bin/python3.8 "$PROJECT/pipeline/annotation/regen_funannotate_outputs.py" \
    "$work/genome.tbl" "$res/$id.scaffolds.fa" "$work/regen" "$id" "$res/$id.stats.json" "${transl:-1}"
for f in gff3 proteins.fa mrna-transcripts.fa cds-transcripts.fa agp contigs.fsa stats.json; do
    [ -s "$work/regen/$id.$f" ] || { echo "ERROR: regeneration produced no $id.$f" >&2; exit 1; }
done

# Install: keep the first copy of each replaced file as .bak
install_file() { # <src> <dest>
    [ -f "$2.bak" ] || { [ ! -f "$2" ] || cp -p "$2" "$2.bak"; }
    cp "$1" "$2"
}
install_file genome.sqn "$res/$id.sqn"
echo "file copied to $res/$id.sqn"
cp discrepency.report.txt "$res/$id.discrepency.report.txt" 2>/dev/null || true
cp internal_stops.tsv "$res/$id.internal_stops.tsv"
install_file genome.gbf "$res/$id.gbk"
install_file genome.tbl "$res/$id.tbl"
for f in gff3 proteins.fa mrna-transcripts.fa cds-transcripts.fa agp contigs.fsa stats.json; do
    install_file "$work/regen/$id.$f" "$res/$id.$f"
done
echo "regenerated: $(grep -c $'\tgene\t' "$res/$id.gff3") genes, $(grep -c '>' "$res/$id.proteins.fa") proteins"
echo "done: $res/$id.sqn"
grep -m1 "Annotated using" "$res/$id.sqn" || true
grep -m1 "Annotated using" "$res/$id.gbk" || true
grep -c "contact" "$res/$id.sqn" || true
