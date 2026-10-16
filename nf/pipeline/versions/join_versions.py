#!/usr/bin/env python3
"""Join tool-version records from nf_aaftf and nf_funannotate1 into one TSV per sample.

nf_aaftf records are keyed by ASMID (e.g. MTK05):
    <aaftf-root>/<ASMID>/versions/*.versions.yml
nf_funannotate1 records are keyed by SPECIES_STRAIN (e.g. Verticillium_dahliae_MD4):
    <annot-root>/<SPECIES_STRAIN>/versions/*.versions.yml
samples_funannotate.csv (SPECIES, STRAIN, ASMID) links the two.

Writes OUT_DIR/<ASMID>.tool_versions.tsv and OUT_DIR/all_samples.tool_versions.tsv.
Columns: asmid, annotation_id, workflow, step, tool, version, container,
pipeline_rev, host, recorded, retroactive, note.

Usage (from the project root):
    pipeline/versions/join_versions.py \
        --samples samples_funannotate.csv \
        --aaftf-root results/sample_versions \
        --annot-root results/annotation \
        --out-dir results/tool_versions
"""
import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merge_versions import parse  # noqa: E402

COLS = ["asmid", "annotation_id", "workflow", "step", "tool", "version", "container",
        "pipeline_rev", "host", "recorded", "retroactive", "note"]


def read_records(vdir, workflow):
    rows = []
    for f in sorted(glob.glob(os.path.join(vdir, "*.versions.yml"))):
        proc, meta, tools = parse(f)
        step = proc.split(":")[-1]
        for tool, ver in tools:
            rows.append([workflow, step, tool, ver, meta.get("container", ""),
                         meta.get("pipeline", ""), meta.get("host", ""),
                         meta.get("recorded", ""), meta.get("retroactive", "false"),
                         meta.get("note", "")])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--aaftf-root", required=True)
    ap.add_argument("--annot-root", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    with open(a.samples, newline="") as fh:
        samples = [(r["ASMID"], r["SPECIES"].replace(" ", "_") + "_" + r["STRAIN"])
                   for r in csv.DictReader(fh)]

    allrows = []
    for asmid, ann_id in samples:
        rows = read_records(os.path.join(a.aaftf_root, asmid, "versions"), "nf_aaftf")
        n_aaftf = len(rows)
        rows += read_records(os.path.join(a.annot_root, ann_id, "versions"), "nf_funannotate1")
        n_ann = len(rows) - n_aaftf
        out = [[asmid, ann_id] + r for r in rows]
        with open(os.path.join(a.out_dir, f"{asmid}.tool_versions.tsv"), "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            fh.writelines("\t".join(r) + "\n" for r in out)
        allrows += out
        print(f"{asmid}\t{ann_id}\tnf_aaftf={n_aaftf}\tnf_funannotate1={n_ann}")
        if not n_aaftf or not n_ann:
            print(f"  WARN: no {'nf_aaftf' if not n_aaftf else 'nf_funannotate1'} records for {asmid}",
                  file=sys.stderr)
    with open(os.path.join(a.out_dir, "all_samples.tool_versions.tsv"), "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(r) + "\n" for r in allrows)
    print(f"wrote {len(allrows)} records for {len(samples)} samples to {a.out_dir}")


if __name__ == "__main__":
    main()
