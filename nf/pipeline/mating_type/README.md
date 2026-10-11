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

The result depends on the MATPredict reference database at the time of the run. See [Reference database](#reference-database).

## How to run

Run from `nf/`:

```bash
sbatch --array=1-4 pipeline/mating_type/01_matpredict.sh   # one task per genome
bash pipeline/mating_type/02_summarize.sh                  # after all tasks finish
bash pipeline/mating_type/02_summarize.sh results/mating_type_final   # optional: another results directory
```

Each task uses 4 CPUs, 16 GB and a 2.5 h limit. These are the limits of `run_detection_batch.slurm` in MATPredict. They have not been tuned for these genomes.

The script sets `MATPREDICT_DB_ROOT` and `MATPREDICT_CACHE_DIR` to the MATPredict checkout (`~/projects/MATPredict`). The database in that checkout is what the run uses. It sets `MATPREDICT_NCBI_EMAIL` to the user's address if the variable is unset.

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
- `withheld_loci` counts candidate loci that the evidence gates rejected. It is not a count of MAT loci.
- The identity of the core MAT hit (in `gene_evidence`) shows how close the genome is to the reference that matched. Low identity means a distant reference, not necessarily a wrong call.

## Reference database

At the first run (2026-10-10), the database had no *Verticillium* or Glomerellales record. The core MAT hits matched references from other orders at 34 to 50% identity.

Three *V. dahliae* records were then curated and accepted into `~/projects/MATPredict/db/Ascomycota/Glomerellales/` (2026-10-11, signed off by J. Stajich):

| Record | Source | Content |
|---|---|---|
| `27337_tv103_MAT_MAT1-2` | GenBank `AB505214.1`, strain TV103 (Usami et al. 2009, doi:10.1007/s10327-009-0197-6) | MAT1-2-1 only; the deposit does not annotate flanking genes |
| `27337_nbrc6119_MAT_MAT1-1` | GenBank `AB505215.1`, strain NBRC 6119 (same paper; the strain is named in the Fig. 2 legend, not on the deposit) | MAT1-1-1 only |
| `27337_vdls-17_MAT_MAT1-2` | RefSeq `GCF_000150675.1`, scaffold `NW_009276925.1`, strain VdLs.17 (genome-derived) | Full locus: SLA2, MAT1-2-1, APN2, COX13 |

Notes on the records:

- The VdLs.17 MAT1-2-1 protein (`XP_009651392.1`, 232 aa) contains the TV103 protein (`BAH66364.1`, 156 aa) at 100% identity over all 156 aa. The RefSeq model is 76 aa longer at the N terminus.
- *V. dahliae* JR2 (`GCA_000400815.2`, `CP009079.1:3830361-3865909`) has the same arrangement and a MAT1-2 call. It is not a record.
- Both Usami deposits come from one paper. NBRC 6119 is the isolate that lacks MAT1-2-1 (PCR and Southern blot), and the MAT1-1 idiomorph was cloned from its genomic library.
- Short et al. 2014 (PMID 25383550) report that 12 of 1,120 *V. dahliae* isolates (1%) were MAT1-1. The MD4 MAT1-1 call is therefore the rare type.

Eight more records from other species are proposed in `db/candidates/Ascomycota/` but not accepted, so they do not affect these runs: *Fulvia fulva* (2), *Penicillium chrysogenum* (2), *Eupenicillium crustaceum* (2, homothallic) and *Fusarium sacchari* (2).

## Status

| Run | Database | Directory | Job |
|---|---|---|---|
| Baseline | no *Verticillium* record | `results/mating_type/` | 29819424 |
| With 2 Usami records | frozen worktree `run-e9060daf32-vdahliae` | `results/mating_type_vref/` | 29858243 |
| Final | live database, 3 *V. dahliae* records | `results/mating_type_final/` | 29858995 |

All tasks finished and exited 0 (3 to 7 min each). The summary table of the final run is `results/mating_type_final/mating_type_summary.tsv`.

## Results

The call is the same in all three runs. The evidence changes.

| Sample | Call | Locus in final run (contig:start-end) | Core hit identity: baseline → final | Margin: baseline → final | Withheld loci: baseline → final |
|---|---|---|---|---|---|
| Verticillium_dahliae_MD4 | MAT1-1 | MTK05_279:10347-28213 | MAT1-1-1: 38.4% → 100% | 46.6 → 607.9 | 12 → 16 |
| Verticillium_dahliae_WVU | MAT1-2 | MTK08_229:12812-30522 | MAT1-2-1: 43.2% → 100% | 18.1 → 63.8 | 13 → 16 |
| Verticillium_nonalfalfae_S13 | MAT1-2 | MTK06_1:270937-314050 | MAT1-2-1: 44.1% → 94.8% | 64.0 → 255.4 | 13 → 18 |
| Verticillium_nonalfalfae_S32 | MAT1-2 | MTK07_16:4337-32180 | MAT1-2-1: 44.1% → 94.8% | 64.0 → 255.4 | 15 → 20 |

All four: confidence `high`, routing `lineage`, `two_idiomorphs` empty. In the final run, the core hits have status `polished_agree`.

### Independent check (BLAST, no MATPredict)

The two Usami deposits were searched against the four genomes with `blastn` and `tblastn` (2026-10-11).

| Sample | Best nucleotide match | Best protein match |
|---|---|---|
| MD4 | `AB505215` (MAT1-1), 99.6% over 5,849 bp (the whole deposit) | MAT1-1-1, 99.7% and 100% (two hits, 311 and 110 aa) |
| WVU | `AB505214` (MAT1-2), 99.4% over 8,504 bp (the whole deposit) | MAT1-2-1, 90.1% |
| S13, S32 | `AB505214` (MAT1-2), 88.5% over 7,474 bp plus 89.7% over 1,088 bp | MAT1-2-1, 84.3% |

WVU, S13 and S32 have no MAT1-1-1 hit above 29% identity. MD4 has no strong MAT1-2-1 hit in the top four hits. These results agree with all four MATPredict calls. The identity values come from two tools with different alignments (`tblastn` here, `exonerate` in MATPredict), so they are not directly comparable. In both, the *V. nonalfalfae* genomes (S13, S32) have lower identity than the *V. dahliae* genomes, as expected for a reference from another species.

### Remaining weak hits

These are the same in the baseline and final runs.

- MD4 has a short MAT1-2-1 hit (126 bp, 50% identity). It is probably background. It is not flagged as a second idiomorph. In the BLAST check, no MAT1-2-1 hit appears in the top four protein hits for MD4.
- S13 and S32 have a small MAT1-1-1 hit (132 bp, 34%, 12% coverage, status `unpolished`).
- The added references raised the number of withheld loci in all four genomes (see the table). The calls did not change.

## MATPredict development notes

The added references raised the number of withheld loci (see the table above). No call changed in these four genomes.
The cause analysis, the MATPredict regression panels (two runs) and the open curation items are MATPredict work. They are in
`~/projects/MATPredict/docs/notes/2026-10-11_verticillium-reference-records.md`.

## Open items

1. The MD4, WVU, S13 and S32 calls are final for this project. They do not depend on the open MATPredict items.
2. The MATPredict curation work (eight proposed records, regression-panel review) is tracked in the note above.
