#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 128 --mem 96gb --time 2:00:00 -p short -C ryzen --out logs/fasttree_run.%A.log

module load fasttree
conda activate /bigdata/stajichlab/shared/condaenv/phyling
NUM=$(wc -l prefix.tab | awk '{print $1}')
source config.txt

ALN=$PREFIX.subset.${NUM}_taxa.$HMM.cds.fasaln
TREE1=$PREFIX.subset.${NUM}_taxa.$HMM.cds.ft_gtr.tre
TREE2=$PREFIX.subset.${NUM}_taxa.$HMM.cds.ft_gtr_long.tre
if [ ! -s $TREE1 ]; then
	FastTreeMP -nt -gtr -gamma < $ALN > $TREE1
	echo "ALN is $ALN"
fi
	if [ -s $TREE1 ]; then
		perl PHYling_unified/util/rename_tree_nodes.pl $TREE1 prefix.tab > $TREE2
	fi
ALN=$PREFIX.subset.${NUM}_taxa.$HMM.aa.fasaln
TREE1=$PREFIX.subset.${NUM}_taxa.$HMM.aa.ft_lg.tre
TREE2=$PREFIX.subset.${NUM}_taxa.$HMM.aa.ft_lg_long.tre

if [ ! -s $TREE1 ]; then
	FastTreeMP -lg -gamma < $ALN > $TREE1
fi
if [ -s $TREE1 ]; then
		 perl PHYling_unified/util/rename_tree_nodes.pl $TREE1 prefix.tab > $TREE2
fi
