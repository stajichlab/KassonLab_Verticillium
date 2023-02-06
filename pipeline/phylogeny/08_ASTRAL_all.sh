#!/bin/bash -l
#SBATCH -p gpu --gres=gpu:1 -N 1 -n 8 --mem 64gb --out logs/ASTRAL.log

source Phylogeny/config.txt

NUM=$(wc -l Phylogeny/prefix.tab | awk '{print $1}')
MEM=64g
module load java
module load ASTRAL

# usethe subset folder instead

CDSGENETREES=Phylogeny/gene_trees/${PREFIX}.ALL.${NUM}_taxa.$HMM.cds.gene_trees.tre
PEPGENETREES=Phylogeny/gene_trees/${PREFIX}.ALL.${NUM}_taxa.$HMM.aa.gene_trees.tre

CDSCONSTREE=Phylogeny/gene_trees/$(basename $CDSGENETREES .gene_trees.tre)".astral.tre"
PEPCONSTREE=Phylogeny/gene_trees/$(basename $PEPGENETREES .gene_trees.tre)".astral.tre"

if [ ! -f $CDSGENETREES ]; then
    cat Phylogeny/gene_trees/$HMM/*.cds.clipkit.FT.tre > $CDSGENETREES
fi

if [ ! -f $PEPGENETREES ]; then
    cat Phylogeny/gene_trees/$HMM/*.aa.clipkit.FT.tre > $PEPGENETREES
fi

echo "$CDSGENETREES ==> $CDSCONSTREE"

if [ ! -s $CDSCONSTREE ]; then
    java -Xmx${MEM} -jar $ASTRALJAR -i $CDSGENETREES -o $CDSCONSTREE
fi

echo "$PEPGENETREES ==> $PEPCONSTREE"
if [ ! -s $PEPCONSTREE ]; then
    java -Xmx${MEM} -jar $ASTRALJAR -i $PEPGENETREES -o $PEPCONSTREE
fi
conda activate /bigdata/stajichlab/shared/condaenv/phyling

PHYling_Unified/util/rename_tree_nodes.pl $PEPCONSTREE Phylogeny/prefix.tab > Phylogeny/gene_trees/$(basename $PEPCONSTREE .tre).long.tre
PHYling_Unified/util/rename_tree_nodes.pl $CDSCONSTREE Phylogeny/prefix.tab > Phylogeny/gene_trees/$(basename $CDSCONSTREE .tre).long.tre

