#!/usr/bin/env bash
# 05_remove_organelle_contigs.sh — remove depth-OUTLIER contigs (likely organelle/contaminant)
# from the funannotate outputs before tbl2asn (see 04_redo_sqn.sh).
#
# Per sample (rows 2..N of samples_funannotate.csv):
#   1. Read flagged contigs from results/depth/<ASMID>.depth.txt ("** OUTLIER" lines).
#   2. blastn each against lib/MT_Genome_ref.fa to label it mito / non-mito (report only).
#   3. Sequester the sequences and their .tbl feature blocks in results/organelle_contigs/<ASMID>/.
#   4. Remove them from these files (original kept once as <file>.pre_organelle):
#        annotate_results/<id>.scaffolds.fa
#        annotate_results/<id>.tbl
#        annotate_misc/tbl2asn/genome.fsa   (input of 04_redo_sqn.sh)
#        annotate_misc/tbl2asn/genome.tbl   (input of 04_redo_sqn.sh)
#
# Not updated: .gff3, .gbk, .proteins.fa, .cds/.mrna-transcripts.fa, .sqn, .agp, .contigs.fsa.
# Re-running is safe: filtering starts from the .pre_organelle copy.
#
# Usage (from the project dir, nf/):
#   DRY_RUN=1 bash pipeline/annotation/05_remove_organelle_contigs.sh   # report only
#   bash pipeline/annotation/05_remove_organelle_contigs.sh             # apply
#   KEEP="MTK08_406 MTK08_425" bash ...   # flagged contigs to keep (space separated)
set -euo pipefail
PROJECT="${PROJECT:-$PWD}"
DRY_RUN="${DRY_RUN:-0}"
KEEP="${KEEP:-}"
MTREF="$PROJECT/lib/MT_Genome_ref.fa"
OUT="$PROJECT/results/organelle_contigs"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load ncbi-blast/2.16.0+ 2>/dev/null || true
command -v blastn >/dev/null || { echo "ERROR: blastn not found" >&2; exit 1; }

# filter_fasta <in> <out> <ids-file> : write records NOT in ids-file
filter_fasta() {
    awk 'NR==FNR{k[$1]=1;next} /^>/{id=substr($1,2); p=!(id in k)} p' "$3" "$1" > "$2"
}
# filter_tbl <in> <out> <ids-file> : drop "Feature <id>" blocks in ids-file
filter_tbl() {
    awk 'NR==FNR{k[$1]=1;next} /^>Feature/{p=!($2 in k)} p' "$3" "$1" > "$2"
}
# extract_fasta / extract_tbl : the complement
extract_fasta() { awk 'NR==FNR{k[$1]=1;next} /^>/{id=substr($1,2); p=(id in k)} p' "$3" "$1" > "$2"; }
extract_tbl()   { awk 'NR==FNR{k[$1]=1;next} /^>Feature/{p=($2 in k)} p' "$3" "$1" > "$2"; }

# apply <file> <kind> <ids-file> : back up once, then rewrite from the backup
apply() {
    local f="$1" kind="$2" ids="$3" bak="$1.pre_organelle"
    [ -f "$f" ] || { echo "  skip (missing): $f"; return; }
    [ -f "$bak" ] || cp -p "$f" "$bak"
    "filter_$kind" "$bak" "$f.tmp" "$ids"
    mv "$f.tmp" "$f"
}

tail -n +2 "$PROJECT/samples_funannotate.csv" | while IFS=, read -r species strain asmid _; do
    id="${species// /_}_${strain}"
    depth="$PROJECT/results/depth/$asmid.depth.txt"
    base="$PROJECT/results/annotation/$id"
    res="$base/annotate_results"
    [ -f "$depth" ] || { echo "$asmid: no depth report, skipping"; continue; }
    [ -f "$res/$id.scaffolds.fa" ] || { echo "$asmid: no scaffolds in $res, skipping"; continue; }

    sample_out="$OUT/$asmid"
    mkdir -p "$sample_out"
    ids="$sample_out/remove.ids"
    # depth report rows: "<contig> <len with commas> <depth> ** OUTLIER (...)"
    grep -F '** OUTLIER' "$depth" | awk '{print $1}' > "$sample_out/flagged.ids"
    : > "$ids"
    for c in $(cat "$sample_out/flagged.ids"); do
        case " $KEEP " in *" $c "*) echo "$asmid: keeping $c (KEEP)";; *) echo "$c" >> "$ids";; esac
    done

    # Source for the mito check and sequestering: the pre-filter scaffolds.
    src="$res/$id.scaffolds.fa"; [ -f "$src.pre_organelle" ] && src="$src.pre_organelle"
    extract_fasta "$src" "$sample_out/$asmid.removed.fa" "$ids"
    blastn -task blastn -query "$sample_out/$asmid.removed.fa" -subject "$MTREF" -evalue 1e-5 \
        -outfmt '6 qseqid qlen length pident' > "$sample_out/$asmid.mito_blast.tsv"

    # report: contig, length, depth, mito-aligned bp (union not computed; sum of HSPs), max %id, genes in .tbl
    tblsrc="$res/$id.tbl"; [ -f "$tblsrc.pre_organelle" ] && tblsrc="$tblsrc.pre_organelle"
    {
        printf 'sample\tcontig\tlength\tdepth\tmito_hsp_bp\tmito_max_pident\tn_genes_in_tbl\n'
        while read -r c; do
            len=$(grep -F "$c " "$depth" | grep -F '** OUTLIER' | awk -v c="$c" '$1==c{gsub(",","",$2); print $2}')
            dep=$(grep -F "$c " "$depth" | grep -F '** OUTLIER' | awk -v c="$c" '$1==c{print $3}')
            hsp=$(awk -v c="$c" '$1==c{s+=$3; if($4>m)m=$4} END{printf "%d\t%.1f", s, m}' "$sample_out/$asmid.mito_blast.tsv")
            ng=$(awk -v c="$c" '/^>Feature/{p=($2==c);next} p && $3=="gene"' "$tblsrc" | wc -l)
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$asmid" "$c" "$len" "$dep" "$hsp" "$ng"
        done < "$ids"
    } > "$sample_out/$asmid.report.tsv"
    echo "== $asmid ($id)"; column -t -s$'\t' "$sample_out/$asmid.report.tsv"

    [ "$DRY_RUN" = 1 ] && continue
    extract_tbl "$tblsrc" "$sample_out/$asmid.removed.tbl" "$ids"
    apply "$res/$id.scaffolds.fa" fasta "$ids"
    apply "$res/$id.tbl" tbl "$ids"
    apply "$base/annotate_misc/tbl2asn/genome.fsa" fasta "$ids"
    apply "$base/annotate_misc/tbl2asn/genome.tbl" tbl "$ids"
    echo "  records now: scaffolds=$(grep -c '>' "$res/$id.scaffolds.fa") tbl=$(grep -c '^>Feature' "$res/$id.tbl") (removed $(wc -l < "$ids"))"
done
