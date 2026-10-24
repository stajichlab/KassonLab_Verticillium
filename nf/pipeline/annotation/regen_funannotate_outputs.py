#!/usr/bin/env python3
"""Regenerate funannotate annotate_results files from a (filtered) .tbl and scaffolds FASTA.

Runs inside the funannotate-1.9.0 container (python3.8 with the funannotate package).
Mirrors the end of `funannotate annotate`:
  lib.tbl2allout -> .gff3 .proteins.fa .mrna-transcripts.fa .cds-transcripts.fa
  parse_scaffolds_makeagp -> .agp .contigs.fsa
  lib.annotation_summary -> .stats.json (command/organism/software kept from the old file)

usage: regen_funannotate_outputs.py TBL SCAFFOLDS OUTDIR PREFIX OLD_STATS_JSON TRANSL_TABLE
Writes OUTDIR/PREFIX.{gff3,proteins.fa,mrna-transcripts.fa,cds-transcripts.fa,agp,contigs.fsa,stats.json}
and OUTDIR/PREFIX.scaffolds.fa.regen (a copy of SCAFFOLDS, not used).
"""
import json
import os
import sys

from funannotate import library as lib
from funannotate.aux_scripts.fasta2agp import parse_scaffolds_makeagp

tbl, fasta, outdir, prefix, old_stats, transl = sys.argv[1:7]
out = os.path.join(outdir, prefix)

# every contig in the tbl must be in the FASTA, or CDS extraction is wrong
fa_ids = set()
with open(fasta) as fh:
    for line in fh:
        if line.startswith(">"):
            fa_ids.add(line[1:].split()[0])
with open(tbl) as fh:
    missing = {l.split()[1] for l in fh if l.startswith(">Feature")} - fa_ids
if missing:
    sys.exit("ERROR: %d tbl contigs missing from %s, e.g. %s" % (len(missing), fasta, sorted(missing)[:3]))

lib.tbl2allout(
    tbl, fasta,
    out + ".gff3", out + ".proteins.fa", out + ".mrna-transcripts.fa",
    out + ".cds-transcripts.fa", out + ".scaffolds.fa.regen",
    transl_table=int(transl),
)
parse_scaffolds_makeagp(fasta, out + ".agp", out + ".contigs.fsa")

new_stats = out + ".stats.json"
lib.annotation_summary(fasta, new_stats, tbl=tbl)
with open(new_stats) as fh:
    new = json.load(fh)
if os.path.isfile(old_stats):
    with open(old_stats) as fh:
        old = json.load(fh)
    for key in ("command", "organism", "software"):
        if key in old:
            new[key] = old[key]
with open(new_stats, "w") as fh:
    json.dump(new, fh, indent=4)
print("regenerated outputs for", prefix)
