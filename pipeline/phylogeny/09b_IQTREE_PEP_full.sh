#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 16 --mem 64gb -p intel --out logs/fastIQTREE_run.%A.log

module load miniconda3
conda activate /bigdata/stajichlab/shared/condaenv/phyling
module load iqtree/2.2.0
NUM=$(wc -l Phylogeny/prefix.tab | awk '{print $1}')
source Phylogeny/config_1.txt

ALN=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.fasaln
TREE1=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.iqft.contree
TREE2=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.iqft_long.tre

if [ ! -s $TREE1 ]; then
	iqtree2 -fast -nt AUTO -m MFP -mset WAG,LG -pre Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.iqft -s $ALN
fi
if [ -s $TREE1 ]; then
	perl PHYling_Unified/util/rename_tree_nodes.pl $TREE1 Phylogeny/prefix.tab > $TREE2
fi
