# Mating-type prediction with MATPredict

Date set up: 2026-10-10. Tool: MATPredict v0.6.1 (`~/projects/MATPredict`, pixi environment).

## Goal

Predict the mating type (MAT1-1 or MAT1-2) of the four annotated *Verticillium* genomes.

## Input

Funannotate scaffold FASTA files: `results/annotation/*/predict_results/*.scaffolds.fa`.

| Sample | Species | NCBI taxid |
|---|---|---|
| Verticillium_dahliae_MD4 | *V. dahliae* | 27337 |
| Verticillium_dahliae_WVU | *V. dahliae* | 27337 |
| Verticillium_nonalfalfae_S13 | *V. nonalfalfae* | 1051616 |
| Verticillium_nonalfalfae_S32 | *V. nonalfalfae* | 1051616 |

The taxids come from `pipeline/annotation/funannotate_meta.tsv`. The sheet is `matpredict_samples.tsv`.

## Method

`matpredict detect` runs per genome:

1. It routes by taxid to the curated MAT families. Both species fall under the Ascomycota `MAT` family (scope: leotiomyceta).
2. It searches the genome with tblastn, using the curated MAT and flanking proteins.
3. It clusters hits, builds gene models (miniprot/exonerate) and applies the evidence gates.
4. It types the locus from the gene content (MAT1-1-1 alpha-box, MAT1-2-1 HMG).

The search uses the genome, not the funannotate proteins. A MAT gene missed by the annotation can still be found.

## How to run

Run from `nf/`:

```bash
sbatch --array=1-4 pipeline/mating_type/01_matpredict.sh   # one task per genome
bash pipeline/mating_type/02_summarize.sh                  # after all tasks finish
```

Each task uses 4 CPUs, 16 GB and a 2.5 h limit. These are the limits of `run_detection_batch.slurm` in MATPredict. They have not been tuned for these genomes.

The script sets `MATPREDICT_DB_ROOT` and `MATPREDICT_CACHE_DIR` to the MATPredict checkout. It sets `MATPREDICT_NCBI_EMAIL` to the user's address if the variable is unset.

## Output

- `results/mating_type/<sample>/detection_report.yaml`: calls, evidence, withheld loci, `two_idiomorphs` flag.
- `results/mating_type/<sample>/detected_loci.gff3`: gene features at each locus.
- `results/mating_type/<sample>/report.html`: readable report.
- `results/mating_type/mating_type_summary.tsv`: one row per called locus.
- Logs: `logs/slurm/matpredict_<jobid>_<task>.{out,err}`.

## How to read the calls

- `idiomorph` is MAT1-1, MAT1-2 or `undetermined`. Check `confidence` and `locus_class`.
- A genome with no call gets `none` in the summary. Read `routing_mode` and `withheld_loci` before you treat `none` as a real absence.
- A non-empty `two_idiomorphs` is a flag for review. It can mean a homothallic strain, or a collapsed or mixed assembly.
- *Verticillium* (Glomerellales) has no curated record of its own. Records of related Hypocreales and Sordariales genomes carry the search. Lineage distance can reduce sensitivity. This has not been measured for *Verticillium* here.

## Status

Job 29819424 submitted 2026-10-10. It was pending (`AssocGrpCpuLimit`) when this was written. No results yet.

## Results

To be filled in after the run (do not enter values that are not in `mating_type_summary.tsv`).
