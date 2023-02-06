#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 48 --mem 64gb -p intel --out logs/raxml_run.%A.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load miniconda3
conda activate /bigdata/stajichlab/shared/condaenv/phyling
module load raxml
OUTGROUP=CatherCBS279.70

NUM=$(wc -l Phylogeny/prefix.tab | awk '{print $1}')
source Phylogeny/config_1.txt

ALN=$PREFIX.${NUM}_taxa.$HMM.aa.fasaln
FULLPREF=$PREFIX.${NUM}_taxa.$HMM.aa.raxml
TREE1=$PREFIX.${NUM}_taxa.$HMM.aa.raxml
TREE2=$PREFIX.${NUM}_taxa.$HMM.aa.iqft_long.tre
pushd Phylogeny
raxmlHPC-PTHREADS-AVX2 -s $ALN -m PROTGAMMAAUTO -n $FULLPREF  -T $CPU -b 121 -# autoFC -p 227 -o $OUTGROUP

#raxml-ng --all --msa testAA.fa --model LG+G8+F --tree pars{10} --bs-trees 200

# if [ ! -s $TREE2 ]; then
#	perl PHYling_Unified/util/rename_tree_nodes.pl $TREE1 Phylogeny/prefix.tab > $TREE2
#fi
