#!/usr/bin/env bash
# 00_make_samples_funannotate.sh — build the nf_funannotate1 sample sheet.
#
# Inputs:
#   samples.csv                          Sample,IlluminaSample,Species,SeqStrain,Strain,
#                                        Project,Description,AssemblyFocus,BioProject,
#                                        BioSample,LOCUSTAG (NCBI locus tag prefix)
#   pipeline/annotation/funannotate_meta.tsv   per-Sample values samples.csv lacks:
#                                        BUSCO_LINEAGE, NCBI_TAXONID
# Output: SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME,SBT,BIOSAMPLE
#
# Fixed values:
#   ASMID        = Sample (unique per row; strain names are not unique, e.g. Chile4)
#   TRANSL_TABLE = 1 (standard genetic code, used by fungi)
#   GENOME       = nf_aaftf sorted assembly, <outdir>/sort/<Sample>.sorted.fasta.gz
#   SBT          = <BioSample>.sbt when samples.csv has a BioSample for the row, else blank.
#                  The pipeline looks the name up under --sbt_dir (lib/sbt) and falls back to
#                  --sbt_template if the file is missing.
#   BIOSAMPLE    = BioSample from samples.csv; the pipeline checks the SBT names the same one.
#
# Exits non-zero if any LOCUSTAG, BUSCO_LINEAGE, or NCBI_TAXONID is blank.
#
# Usage: 00_make_samples_funannotate.sh <samples.csv> <out.csv> [meta.tsv] [aaftf_results_dir]

set -euo pipefail

IN="${1:?usage: $0 <samples.csv> <out.csv> [meta.tsv] [aaftf_results_dir]}"
OUT="${2:?usage: $0 <samples.csv> <out.csv> [meta.tsv] [aaftf_results_dir]}"
HERE="$(cd "$(dirname "$IN")" && pwd)"
META="${3:-$HERE/pipeline/annotation/funannotate_meta.tsv}"
RESULTS="${4:-$HERE/results}"

[ -f "$IN" ]   || { echo "missing samples file: $IN" >&2; exit 1; }
[ -f "$META" ] || { echo "missing metadata file: $META" >&2; exit 1; }

# Map Sample -> "BUSCO_LINEAGE NCBI_TAXONID" from the metadata TSV.
declare -A BUSCO TAXID
while IFS=$'\t' read -r S B T; do
    [ "$S" = "Sample" ] && continue
    [ -n "$S" ] || continue
    BUSCO[$S]="$B"; TAXID[$S]="$T"
done < <(tr -d '\r' < "$META")

tmp="$OUT.tmp"
rm -f "$tmp.bad"
{
    echo "SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME,SBT,BIOSAMPLE"
    tail -n +2 "$IN" | tr -d '\r' | while IFS=, read -r SAMPLE _ILLUMINA SPECIES STRAIN _OFFICIAL _PROJECT _DESC _FOCUS _BIOPROJECT BIOSAMPLE LOCUSTAG; do
        [ -n "$SAMPLE" ] || continue
        L="$LOCUSTAG"; B="${BUSCO[$SAMPLE]:-}"; T="${TAXID[$SAMPLE]:-}"
        if [ -z "$L" ] || [ -z "$B" ] || [ -z "$T" ]; then
            echo "incomplete metadata for $SAMPLE in $META (LOCUSTAG='$L' BUSCO_LINEAGE='$B' NCBI_TAXONID='$T')" >&2
            # The loop runs in a subshell, so record failures in a file.
            echo "$SAMPLE" >> "$tmp.bad"
        fi
        echo "${SPECIES},${STRAIN},${SAMPLE},${L},${B},1,${T},${RESULTS}/sort/${SAMPLE}.sorted.fasta.gz,${BIOSAMPLE:+${BIOSAMPLE}.sbt},${BIOSAMPLE}"
    done
} > "$tmp"

if [ -s "$tmp.bad" ]; then
    echo "aborting: fill LOCUSTAG, BUSCO_LINEAGE, NCBI_TAXONID for: $(paste -sd' ' "$tmp.bad")" >&2
    rm -f "$tmp" "$tmp.bad"
    exit 1
fi
mv "$tmp" "$OUT"
echo "wrote $OUT ($(tail -n +2 "$OUT" | wc -l) samples)"
