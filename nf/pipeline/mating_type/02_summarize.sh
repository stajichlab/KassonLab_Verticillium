#!/usr/bin/env bash
# 02_summarize.sh — collect the MATPredict calls into results/mating_type/mating_type_summary.tsv.
# One row per called locus. A genome with no called locus gets one row with call "none".
# Run from any directory after 01_matpredict.sh has finished.
set -euo pipefail

PROJ=/bigdata/stajichlab/shared/projects/KassonLab/Verticillium/nf
OUTROOT=$PROJ/results/mating_type

/usr/bin/python3.12 - "$OUTROOT" <<'PY'
import sys, glob, os, yaml
root = sys.argv[1]
cols = ["sample", "routing_mode", "family", "idiomorph", "confidence", "locus_class",
        "contig", "start", "end", "two_idiomorphs", "withheld_loci"]
rows = []
for rep in sorted(glob.glob(os.path.join(root, "*", "detection_report.yaml"))):
    sample = os.path.basename(os.path.dirname(rep))
    d = yaml.safe_load(open(rep))
    base = {"sample": sample, "routing_mode": d.get("routing_mode"),
            "two_idiomorphs": ",".join(map(str, d.get("two_idiomorphs") or [])) or "no",
            "withheld_loci": len(d.get("suppressed_loci") or [])}
    det = d.get("detected") or []
    if not det:
        rows.append({**base, "family": "NA", "idiomorph": "none"})
    for l in det:
        rows.append({**base, **{k: l.get(k) for k in
                     ("family", "idiomorph", "confidence", "locus_class", "contig", "start", "end")}})
out = os.path.join(root, "mating_type_summary.tsv")
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join("" if r.get(c) is None else str(r[c]) for c in cols) + "\n")
print(f"wrote {out} ({len(rows)} rows)")
PY
