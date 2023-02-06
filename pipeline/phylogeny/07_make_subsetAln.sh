#!/bin/bash -l
#SBATCH -p short --out logs/make_subsetAln.log

BINFILES=Phylogeny/aln_bin/bin_lists
mkdir -p $BINFILES
Rscript scripts/explore_trees.R $BINFILES
conda activate /bigdata/stajichlab/shared/condaenv/phyling

source config.txt

subsetdirCDS=Phylogeny/aln_bin/CDS
subsetdirPEP=Phylogeny/aln_bin/PEP
subsetdirBIN=Phylogeny/aln_bin/bins

mkdir -p $subsetdirCDS $subsetdirPEP

#for n in $(seq 10 10 100)
for n in 380
do
 PEPIN=$BINFILES/filtered_q${n}LenEvoRate_AA.tsv
 CDSIN=$BINFILES/filtered_q${n}LenEvoRate_CDS.tsv

 cut -f1,2 $PEPIN | tail -n +2 | while read TREE ALN
 do
     rsync -a ${ALN_OUTDIR}/$HMM/$TREE $subsetdirPEP/${n}
     rsync -a ${ALN_OUTDIR}/$HMM/$ALN $subsetdirPEP/${n}
 done

 cut -f1,2 $CDSIN | tail -n +2 | while read TREE ALN
 do
     rsync -a ${ALN_OUTDIR}/$HMM/$TREE $subsetdirCDS/${n}
     rsync -a ${ALN_OUTDIR}/$HMM/$ALN $subsetdirCDS/${n}
 done

 taxa=$(wc -l prefix.tab | awk '{print $1}')
 ./PHYling_Unified/util/combine_multiseq_aln.py -d $subsetdirCDS/${n} --moltype DNA --expected expected_prefixes.lst --ext cds.clipkit \
						-p $subsetdirBIN/$PREFIX.subset.${taxa}_taxa.$HMM.cds.partitions.txt \
						-o $PREFIX.subset.${taxa}_taxa.$HMM.cds.fasaln
 perl -i -p -e 's/DNA/GTR/' $subsetdirBIN/$PREFIX.subset.${taxa}_taxa.$HMM.cds.partitions.txt
 ./PHYling_Unified/util/combine_multiseq_aln.py -d $subsetdirPEP/${n} --expected expected_prefixes.lst --ext aa.clipkit \
						-p $subsetdirBIN/$PREFIX.subset.${taxa}_taxa.$HMM.aa.partitions.txt \
						-o $PREFIX.subset.${taxa}_taxa.$HMM.aa.fasaln
 perl -i -p -e 's/PROT/WAG/' $subsetdirBIN/$PREFIX.subset.${taxa}_taxa.$HMM.aa.partitions.txt

done
