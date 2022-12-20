#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 2 --mem 4gb --out logs/stats.log

module load AAFTF

SAMPLEFILE=samples.csv
INDIR=asm
OUTDIR=genomes

mkdir -p $OUTDIR
IFS=, # set the delimiter to be ,
tail -n +2 $SAMPLEFILE | while read ID BASE SPECIES STRAIN LOCUSTAG TYPESTRAIN
do
    for type in AAFTF
    do
	if [ ! -f $INDIR/$type/$ID.sorted.fasta ]; then
		continue
	fi
	rsync -a $INDIR/$type/$ID.sorted.fasta $OUTDIR/$ID.$type.fasta
	if [[ ! -f $OUTDIR/$ID.$type.stats.txt || $OUTDIR/$ID.$type.fasta -nt $OUTDIR/$ID.$type.stats.txt ]]; then
    	    AAFTF assess -i $OUTDIR/$ID.$type.fasta -r $OUTDIR/$ID.$type.stats.txt
	fi
    done
done
