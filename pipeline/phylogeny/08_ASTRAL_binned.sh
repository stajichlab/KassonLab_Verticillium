#!/bin/bash -l
#SBATCH -p gpu --gres=gpu:1 -N 1 -n 16 --mem 64gb --out logs/ASTRAL_binned.log

source Phylogeny/config.txt

NUM=$(wc -l Phylogeny/prefix.tab | awk '{print $1}')
MEM=64g
module load java
module load ASTRAL

# could make another loop with cds/aa but I think I sometimes am inconsistent with AA vs PEP 
# usethe subset folder instead
for PART in 75p 90p
do
    CDSGENETREES=Phylogeny/gene_trees/${PREFIX}.binned_$PART.${NUM}_taxa.$HMM.cds.gene_trees.tre
    PEPGENETREES=Phylogeny/gene_trees/${PREFIX}.binned_$PART.${NUM}_taxa.$HMM.aa.gene_trees.tre
    
    CDSCONSTREE=Phylogeny/gene_trees/$(basename $CDSGENETREES .gene_trees.tre)".astral.tre"
    PEPCONSTREE=Phylogeny/gene_trees/$(basename $PEPGENETREES .gene_trees.tre)".astral.tre"
    PEPTREESCORE=Phylogeny/gene_trees/$(basename $PEPGENETREES .gene_trees.tre)".astral_scores.tre"
    CDSTREESCORE=Phylogeny/gene_trees/$(basename $CDSGENETREES .gene_trees.tre)".astral_scores.tre"
    
    if [ ! -f $CDSGENETREES ]; then
	cat Phylogeny/gene_trees/binned_trees/$HMM/$PART/*.cds.clipkit.FT.tre > $CDSGENETREES
    fi
    if [ ! -f $PEPGENETREES ]; then
	cat Phylogeny/gene_trees/binned_trees/$HMM/$PART/*.aa.clipkit.FT.tre > $PEPGENETREES
    fi
    
    echo "$CDSGENETREES ==> $CDSCONSTREE"
    
    if [ ! -s $CDSCONSTREE ]; then
	java -Xmx${MEM} -jar $ASTRALJAR -i $CDSGENETREES -o $CDSCONSTREE
	java -Xmx${MEM} -jar $ASTRALJAR -q $CDSCONSTREE -i $CDSGENETREES -o $CDSTREESCORE
    fi
    
    echo "$PEPGENETREES ==> $PEPCONSTREE"
    if [ ! -s $PEPCONSTREE ]; then
	java -Xmx${MEM} -jar $ASTRALJAR -i $PEPGENETREES -o $PEPCONSTREE
	java -Xmx${MEM} -jar $ASTRALJAR -q $PEPCONSTREE -i $PEPGENETREES -o $PEPTREESCORE
    fi

    LONGPEPTREE=Phylogeny/gene_trees/$(basename $PEPCONSTREE .tre).long.tre
    LONGCDSTREE=Phylogeny/gene_trees/$(basename $CDSCONSTREE .tre).long.tre
    if [ ! -s $LONGPEPTREE ]; then
	# need bioperl modules to run this script
	conda activate /bigdata/stajichlab/shared/condaenv/phyling	
	PHYling_Unified/util/rename_tree_nodes.pl $PEPCONSTREE Phylogeny/prefix.tab > $LONGPEPTREE	
	PHYling_Unified/util/rename_tree_nodes.pl $CDSCONSTREE Phylogeny/prefix.tab > $LONGCDSTREE
	conda deactivate    
    fi
done
