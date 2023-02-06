#!/bin/bash -l
#SBATCH -p short -N 1 -n 128 --mem 96gb  -N 1 --out logs/make_CDS_trees.%A.log
module load fasttree
module load iqtree/2.2.0

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
pushd Phylogeny
make -f ../PHYling_Unified/util/makefiles/Makefile.trees HMM=fungi_odb10 -j $CPU CDS
