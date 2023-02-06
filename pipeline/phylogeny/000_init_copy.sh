#!/bin/bash -l
#SBATCH -p short -N 1 -n 1 --mem 1gb --out logs/phylo_init.log

TARGET=Phylogeny/genomes
ZYGOLIFEASM=zygolife
ZYGOLIFELIST=lib/zygolife_strains_phylogeny.tsv
HKIASM=genomes
HKILIST=lib/hki_strains_phylogeny.tsv

mkdir -p $TARGET

IFS=,

cat $HKILIST | while read ID STRAIN SPECIES
do
    ASM=$HKIASM/$ID.AAFTF.fasta
    NAME=$(echo -n "${SPECIES}_${STRAIN}" | perl -p -e 's/\s+/_/g')
    if [ ! -f $TARGET/$NAME.dna.fasta ]; then
	rsync -a $ASM $TARGET/$NAME.dna.fasta
    fi
done

cat $ZYGOLIFELIST | while read SPECIES
do
    ASM=$ZYGOLIFEASM/$SPECIES.sorted.fasta
    if [ -s $ASM ]; then
	if [ ! -f $TARGET/$SPECIES.asm.fasta ]; then
	    rsync -aL $ASM $TARGET/$SPECIES.dna.fasta
	fi
    else
	echo "did not find $ASM"
    fi
done
