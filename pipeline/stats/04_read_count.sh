#!/usr/bin/bash -l

#SBATCH --nodes 1 --ntasks 24 --mem 24G -p batch -J readcount --out logs/bbcount.%a.log --time 48:00:00
module load BBMap
hostname
MEM=24
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
SAMPLEFILE=samples.csv

ASM=genomes
OUTDIR=$(realpath mapping_report)
SAMPLES=samples.csv
mkdir -p $OUTDIR

IFS=, # set the delimiter to be ,
IFS=, # set the delimiter to be ,
tail -n +2 $SAMPLES | sed -n ${N}p | while read ID BASE SPECIES STRAIN LOCUSTAG TYPESTRAIN
do
    
    LEFT=$(realpath $INDIR/${BASE}_R1.fastq.gz)
    RIGHT=$(realpath $INDIR/${BASE}_R2.fastq.gz)
    
    echo "$LEFT $RIGHT"
    for type in AAFTF
    do
	SORTED=$(realpath $ASM/${ID}.AAFTF.fasta)
	REPORTOUT=${ID}.${type}
	if [ -s $SORTED ]; then
	    pushd $SCRATCH
	    if [ ! -s $OUTDIR/${REPORTOUT}.bbmap_covstats.txt ]; then
		bbmap.sh -Xmx${MEM}g ref=$SORTED in=$LEFT in2=$RIGHT covstats=$OUTDIR/${REPORTOUT}.bbmap_covstats.txt  statsfile=$OUTDIR/${REPORTOUT}.bbmap_summary.txt
	    fi
	    popd
	fi
    done
done
