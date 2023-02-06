#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 32 --mem 64gb -p intel --out logs/fasttree_run.%A.log

conda activate /bigdata/stajichlab/shared/condaenv/phyling
module load fasttree
NUM=$(wc -l Phylogeny/prefix.tab | awk '{print $1}')
source Phylogeny/config.txt

ALN=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.cds.fasaln
TREE1=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.cds.ft_gtr.tre
TREE2=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.cds.ft_gtr_long.tre
if [ ! -s $TREE1 ]; then
	which FastTreeMP
	FastTreeMP -nt -gtr -gamma < $ALN > $TREE1
	echo "ALN is $ALN"
fi
	if [ -s $TREE1 ]; then
		perl PHYling_Unified/util/rename_tree_nodes.pl $TREE1 Phylogeny/prefix.tab > $TREE2
	fi
ALN=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.fasaln
TREE1=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.ft_lg.tre
TREE2=Phylogeny/$PREFIX.${NUM}_taxa.$HMM.aa.ft_lg_long.tre

if [ ! -s $TREE1 ]; then
	FastTreeMP -lg -gamma < $ALN > $TREE1
fi
if [ -s $TREE1 ]; then
		 perl PHYling_Unified/util/rename_tree_nodes.pl $TREE1 Phylogeny/prefix.tab > $TREE2
fi
