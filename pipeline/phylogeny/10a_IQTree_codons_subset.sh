#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 6 --mem 96gb --time 7-0:00:00 -p intel --out logs/iqtree_CDS.%A.log
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load iqtree/2.2.0
NUM=$(wc -l Phylogeny/prefix.tab | awk '{print $1}')
source Phylogeny/config.txt

ALN=Phylogeny/$PREFIX.subset.${NUM}_taxa.$HMM.cds.fasaln
PART=Phylogeny/$PREFIX.subset.${NUM}_taxa.$HMM.cds.partitions.txt

perl -ip -e 's/^DNA/CODON/' $PART
iqtree2 -s $ALN --prefix $PREFIX.subset.${NUM}_taxa.$HMM.cds -T AUTO --threads-max $CPU -p $PART -m TESTMERGE -rcluster 10
