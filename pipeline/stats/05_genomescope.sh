#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 24 --mem 64gb --out logs/genomescope.%a.log -a 1-90

module load workspace/scratch
module load samtools
module load jellyfish
module load R

GENOMESCOPE=genomescope

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi
FASTQFOLDER=input
SAMPLEFILE=samples.csv
MAX=$(wc -l $SAMPLEFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in samplefile=$SAMPLEFILE"
  exit
fi
mkdir -p $GENOMESCOPE
JELLYFISHSIZE=1000000000
IFS=,
KMER=21
READLEN=150 # note this assumes all projects are 150bp reads which they may not be
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE ILLUMINASAMPLE SPECIES STRAIN PROJECT DESCRIPTION ASMFOCUS
do
    # for this project has only a single file base  but this needs fixing otherse
    jellyfish count -C -m $KMER -s $JELLYFISHSIZE -t $CPU -o $SCRATCH/$STRAIN.jf <(pigz -dc $FASTQFOLDER/${BASE}_${ILLUMINASAMPLE}_R[12]_001.fastq.gz)
    jellyfish histo -t $CPU $SCRATCH/$STRAIN.jf > $GENOMESCOPE/$STRAIN.histo
    Rscript scripts/genomescope.R $GENOMESCOPE/$STRAIN.histo $KMER $READLEN $GENOMESCOPE/$STRAIN/
done
