# Tool versions

Per-task version records come from `record_versions.sh` in nf_aaftf and nf_funannotate1
(branch `feat/tool-versions`). Records:

- nf_aaftf: `results/sample_versions/<ASMID>/versions/<STEP>.versions.yml`
- nf_funannotate1: `results/annotation/<Species_Strain>/versions/<STEP>.versions.yml`

Scripts here:

1. `backfill_versions.py` — for steps that finished before recording existed. Runs the
   version commands in today's container and marks the record `retroactive: "true"`.
   The note field gives the image and output mtimes. Never overwrites a real record.
   Run it BEFORE replacing an image (for example the AAFTF release update).
2. `join_versions.py` — joins both pipelines' records into `results/tool_versions/<ASMID>.tool_versions.tsv`
   and `all_samples.tool_versions.tsv`, using `samples_funannotate.csv` to link ASMID to Species_Strain.

```
pipeline/versions/backfill_versions.py --samples samples_funannotate.csv
pipeline/versions/join_versions.py --samples samples_funannotate.csv \
    --aaftf-root results/sample_versions --annot-root results/annotation --out-dir results/tool_versions
```

`merge_versions.py` is a copy of the one in the pipeline repos' `bin/`.
