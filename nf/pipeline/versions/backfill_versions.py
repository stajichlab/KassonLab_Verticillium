#!/usr/bin/env python3
"""Backfill tool-version records for samples that finished before version recording existed.

For each sample and step it finds output files as evidence that the step ran,
then runs the step's version commands inside the container that the pipeline
config names TODAY, and writes a record with `retroactive: "true"`.

Limit: this records the image that exists now, not the one that ran. Each
record's note gives the image mtime and the output mtime, and says whether the
image is newer than the output. If it is newer, the versions may not match the run.
Run this BEFORE replacing an image (for example the AAFTF release update).

Real records (retroactive: "false") are never overwritten.

Tool lists below mirror the record_versions.sh calls in the nf_aaftf and
nf_funannotate1 modules. Keep them in sync.

Usage (from the project root):
    pipeline/versions/backfill_versions.py --samples samples_funannotate.csv [--dry-run]
"""
import argparse
import csv
import datetime
import glob
import os
import subprocess
import sys

SIF_AAFTF = "/bigdata/stajichlab/shared/lib/singularity_cache/AAFTF.sif"
SIF_DIR = "/bigdata/stajichlab/shared/lib/singularity_cache"
SIF_FUN = f"{SIF_DIR}/funannotate-1.9.0.sif"
SIF_ANTISMASH = f"{SIF_DIR}/antismash-standalone-8.0.4-procps.sif"
FCS_ADAPTOR_SIF = "/bigdata/stajichlab/shared/lib/AAFTF_DB/fcs-adaptor.0.5.5.sif"
DEFAULT_HELPER = os.path.expanduser("~/projects/nf/nf_funannotate1/bin/record_versions.sh")

AAFTF = "AAFTF=AAFTF --version"
BBDUK = "bbduk=bbduk.sh --version 2>&1 | grep -i 'BBTools version'"
BWA = "bwa=bwa 2>&1 | grep -i '^Version'"
MM2 = "minimap2=minimap2 --version"
SAM = "samtools=samtools --version | head -1"
GX = ("fcs_gx=command -v gx >/dev/null && ls $(dirname $(command -v gx))/../conda-meta "
      "2>/dev/null | grep ncbi-fcs-gx | head -1")
TAXK = "taxonkit=taxonkit version"
FUN = ["funannotate=funannotate --version", "augustus=augustus --version 2>&1 | head -1",
       "minimap2=minimap2 --version", "samtools=samtools --version | head -1",
       "tantan=tantan --version 2>&1 | head -1", "diamond=diamond --version",
       "bedtools=bedtools --version"]

# (workflow, step, evidence glob, image, tools). {a}=ASMID, {i}=annotation id.
# Paths are relative to --results (AAFTF) or --annot-root / --launch-dir.
AAFTF_STEPS = [
    ("TRIM", "trim/{a}_1P.fastq.gz", SIF_AAFTF, [AAFTF, "fastp=fastp --version", BBDUK, "pigz=pigz --version"]),
    ("FILTER", "filter/{a}_filtered_1.fastq.gz", SIF_AAFTF,
     [AAFTF, BBDUK, BWA, MM2, SAM, "filter_aligner=echo bbduk"]),
    ("ASSEMBLE", "asm/{a}.spades.fasta.gz", SIF_AAFTF,
     [AAFTF, "spades=spades.py --version", "megahit=megahit --version", "assembler=echo spades"]),
    ("RMDUP", "rmdup/{a}.rmdup.fasta.gz", SIF_AAFTF, [AAFTF, MM2, "rmdup_tools=echo minimap2"]),
    ("VECSCREEN", "vecscreen/{a}.vecscreen.fasta.gz", SIF_AAFTF, [AAFTF, "blastn=blastn -version | head -1"]),
    ("CONTAM_CLEAN", "contam_clean/{a}.contam_clean.fasta.gz", SIF_AAFTF,
     [AAFTF, GX, TAXK, "HOST:fcs_gx_db=/srv/projects/db/ncbi-fcs/0.5.4/gxdb"]),
    ("FCS_SCREEN", "fcs_screen/{a}.fcs_screen.fasta.gz", None,
     [f"fcs_adaptor_sif=basename $(readlink -f {FCS_ADAPTOR_SIF})"]),
    ("SORT", "sort/{a}.sorted.fasta.gz", SIF_AAFTF, [AAFTF]),
    ("COMPRESS", "sort/{a}.sorted.fasta.gz", SIF_AAFTF, ["bgzip=bgzip --version | head -1"]),
    ("POLISH", "polish/{a}.polished.fasta.gz", SIF_AAFTF,
     [AAFTF, "pilon=pilon --version 2>&1 | head -1", BWA, SAM, "polisher=echo polca"]),
    ("SOURPURGE", "sourpurge/{a}.sourpurge.fasta.gz", SIF_AAFTF,
     [AAFTF, "sourmash=sourmash --version", TAXK, "HOST:aaftf_db=/bigdata/stajichlab/shared/lib/AAFTF_DB"]),
    ("DEPTH", "depth/{a}.depth.txt", SIF_AAFTF, [AAFTF, MM2, SAM]),
    ("ASSESS", "assess/{a}.stats.txt", SIF_AAFTF, [AAFTF]),
]
ANNOT_STEPS = [
    ("GENOME_CLEAN", "launch:input_clean_genomes/{a}.fa.gz", SIF_AAFTF,
     [AAFTF, TAXK, GX, "pigz=pigz --version"]),
    ("MASKREPEAT_TANTAN_RUN", "launch:input_clean_genomes/{a}.masked.fasta.gz", SIF_FUN, [FUN[0], FUN[4]]),
    ("FUNANNOTATE_TRAIN", "annot:{i}/training/funannotate_train.transcripts.gff3", SIF_FUN,
     FUN + ["hisat2=hisat2 --version | head -1", "trinity=Trinity --version 2>&1 | head -1"]),
    ("FUNANNOTATE_PREDICT", "annot:{i}/predict_results/{i}.gbk", SIF_FUN, FUN),
    ("BUSCO_COMPLETENESS", "annot:{i}/busco_completeness/{i}/short_summary.specific.*.txt",
     f"{SIF_DIR}/depot.galaxyproject.org-singularity-busco-6.1.0--pyhdfd78af_1.img", ["busco=busco --version"]),
    ("ANTISMASH_RUN", "annot:{i}/antismash_local/{i}.gbk", SIF_ANTISMASH,
     ["antismash=antismash --version 2>&1 | head -1"]),
    ("FUNANNOTATE_ANNOTATE", "annot:{i}/annotate_results/{i}.gbk", SIF_FUN, FUN[:2] + [FUN[5]]),
]


def newest_mtime(pattern):
    hits = glob.glob(pattern)
    return max((os.stat(h).st_mtime for h in hits), default=None)


def iso(t):
    return datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime("%Y-%m-%d")


def is_real_record(path):
    try:
        return 'retroactive: "false"' in open(path).read()
    except OSError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="samples_funannotate.csv (SPECIES, STRAIN, ASMID)")
    ap.add_argument("--results", default="results", help="AAFTF results dir")
    ap.add_argument("--annot-root", default="results/annotation")
    ap.add_argument("--launch-dir", default=".nf_launch/annotate")
    ap.add_argument("--aaftf-versions", default="results/sample_versions")
    ap.add_argument("--helper", default=DEFAULT_HELPER)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    helper = os.path.realpath(a.helper)
    if not os.path.exists(helper):
        sys.exit(f"helper not found: {helper}")
    binds = ",".join(sorted({os.path.dirname(helper), "/bigdata", "/scratch"}))
    with open(a.samples, newline="") as fh:
        samples = [(r["ASMID"], r["SPECIES"].replace(" ", "_") + "_" + r["STRAIN"])
                   for r in csv.DictReader(fh)]

    n_written = n_skipped = 0
    for asmid, ann_id in samples:
        for wf, steps in (("NF_AAFTF", AAFTF_STEPS), ("NF_FUNANNOTATE1", ANNOT_STEPS)):
            for step, ev, image, tools in steps:
                root = {"launch:": a.launch_dir, "annot:": a.annot_root}.get(
                    next((p for p in ("launch:", "annot:") if ev.startswith(p)), ""), a.results)
                pat = os.path.join(root, ev.split(":", 1)[-1].format(a=asmid, i=ann_id))
                out_mtime = newest_mtime(pat)
                if out_mtime is None:
                    continue
                dest = (os.path.join(a.aaftf_versions, asmid, "versions") if wf == "NF_AAFTF"
                        else os.path.join(a.annot_root, ann_id, "versions"))
                rec = os.path.join(dest, f"{step}.versions.yml")
                if is_real_record(rec):
                    print(f"SKIP  {asmid} {step}: real record exists")
                    n_skipped += 1
                    continue
                if image and not os.path.exists(image):
                    print(f"SKIP  {asmid} {step}: image not found: {image}")
                    n_skipped += 1
                    continue
                if image:
                    img_real = os.path.realpath(image)
                    img_mtime = os.stat(img_real).st_mtime
                    newer = img_mtime > out_mtime
                    note = (f"retroactive; image {os.path.basename(img_real)} mtime {iso(img_mtime)}; "
                            f"step output mtime {iso(out_mtime)}; image_newer_than_output={str(newer).lower()}")
                    container = img_real
                else:
                    note = f"retroactive; native (no container); step output mtime {iso(out_mtime)}"
                    container = "native"
                # HOST:name=path -> resolve on the host (the path may not be visible in the container)
                tools = [f"{t[5:].split('=')[0]}=echo {os.path.realpath(t.split('=', 1)[1])}"
                         if t.startswith("HOST:") else t for t in tools]
                cmd = [helper, os.path.abspath(dest), f"{wf}:{step}", asmid if wf == "NF_AAFTF" else ann_id,
                       container, "unknown (retroactive)"] + tools
                if image:
                    cmd = ["apptainer", "exec", "--env", "RECORD_RETROACTIVE=1", "--env", f"RECORD_NOTE={note}",
                           "-B", binds, img_real] + cmd
                env = dict(os.environ, RECORD_RETROACTIVE="1", RECORD_NOTE=note)
                print(f"{'PLAN ' if a.dry_run else 'WRITE'} {asmid} {step} -> {rec}  [{note}]")
                if not a.dry_run:
                    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
                    if r.returncode != 0 or not os.path.exists(rec):
                        print(f"  FAILED rc={r.returncode}: {r.stderr.strip()[-300:]}", file=sys.stderr)
                        continue
                n_written += 1
    print(f"{'planned' if a.dry_run else 'written'}={n_written} skipped={n_skipped}")


if __name__ == "__main__":
    main()
