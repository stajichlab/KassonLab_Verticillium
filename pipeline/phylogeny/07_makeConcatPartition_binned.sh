#!/bin/bash -l
#SBATCH -p short --out logs/make_subsetAln.log

TOPDIR=Phylogeny
#Rscript scripts/explore_trees.R 
conda activate /bigdata/stajichlab/shared/condaenv/phyling

source Phylogeny/config.txt

mkdir -p Phylogeny/gene_trees/binned_trees/$HMM Phylogeny/gene_trees/binned_aln/$HMM Phylogeny/gene_trees/combined

subsetdirTREE=$(realpath Phylogeny/gene_trees/binned_trees/$HMM)
subsetdirALN=$(realpath Phylogeny/gene_trees/binned_aln/$HMM)
subsetdirBIN=$(realpath Phylogeny/gene_trees/combined)
mkdir -p $subsetdirTREE $subsetdirALN

ALN_OUTDIR=$(realpath $TOPDIR/$ALN_OUTDIR)
#for n in $(seq 10 10 100)
for n in 75p 90p
do
    # PEPIN=$BINFILES/filtered_q${n}LenEvoRate_AA.tsv
    # CDSIN=$BINFILES/filtered_q${n}LenEvoRate_CDS.tsv

    mkdir -p $subsetdirTREE/${n} $subsetdirALN/${n}
    PEPIN=$TOPDIR/gene_trees/gene_trees.$HMM.filtered_tax${n}_PEP.tsv
    CDSIN=$TOPDIR/gene_trees/gene_trees.$HMM.filtered_tax${n}_CDS.tsv
    
    cut -f1,3 $PEPIN | tail -n +2 | while read TREE ALN
    do
	if [ ! -f $subsetdirTREE/${n}/$TREE ]; then ln -s ${ALN_OUTDIR}/$HMM/$TREE $subsetdirTREE/${n}/; fi
	if [ ! -f $subsetdirALN/${n}/$ALN ]; then ln -s ${ALN_OUTDIR}/$HMM/$ALN $subsetdirALN/${n}/; fi
    done
    
    cut -f1,3 $CDSIN | tail -n +2 | while read TREE ALN
    do
	if [ ! -f $subsetdirTREE/${n}/$TREE ]; then ln -s ${ALN_OUTDIR}/$HMM/$TREE $subsetdirTREE/${n}/; fi
	if [ ! -f $subsetdirALN/${n}/$ALN ]; then ln -s ${ALN_OUTDIR}/$HMM/$ALN $subsetdirALN/${n}/; fi
    done

    pushd Phylogeny
    taxa=$(wc -l prefix.tab | awk '{print $1}')
    if [ ! -s $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.cds.fasaln ]; then
	../PHYling_Unified/util/combine_multiseq_aln.py -d $subsetdirALN/${n} --moltype DNA \
							--expected expected_prefixes.lst --ext cds.clipkit \
							-p $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.cds.partitions.txt \
							-o $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.cds.fasaln
	perl -i -p -e 's/DNA/GTR/' $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.cds.partitions.txt
    fi
    if [ ! -s $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.aa.fasaln ]; then
	../PHYling_Unified/util/combine_multiseq_aln.py -d $subsetdirALN/${n} --expected expected_prefixes.lst \
							--ext aa.clipkit \
						    -p $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.aa.partitions.txt \
						    -o $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.aa.fasaln
	perl -i -p -e 's/PROT/WAG/' $subsetdirBIN/$PREFIX.subset_${n}.${taxa}_taxa.$HMM.aa.partitions.txt
    fi
    popd
done
