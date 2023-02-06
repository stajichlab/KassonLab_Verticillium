#!/bin/bash -l
#SBATCH -p short -N 1 -n 96 --mem 96gb --out logs/phykit_summarize_cds.%A.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

source Phylogeny/config.txt
TREEEXT=cds.clipkit.FT.tre
INDIR=$(realpath Phylogeny/$ALN_OUTDIR/$HMM)
mkdir -p Phylogeny/gene_trees/$HMM
pushd Phylogeny/gene_trees/$HMM
ln -s $INDIR/*.${TREEEXT} .
ln -s $INDIR/*.$(basename $TREEEXT .FT.tre) .
find -L . -size 0 | xargs rm -f

module load parallel
summarize() {
module load phykit
treefile=$(basename $1)
aln=$(basename $treefile .FT.tre)
taxct=$(grep -c "^>" $aln)
len=$(phykit aln_len $aln)
bss=$(phykit bss $treefile)
evorate=$(phykit evolutionary_rate $treefile)
meanBSS=$(echo $bss | grep mean: | awk '{print $2}')
medianBSS=$(echo $bss | grep median: | awk '{print $2}')
treeness=$(phykit treeness $treefile)
echo -e "$treefile\t$taxct\t$aln\t$len\t$meanBSS\t$medianBSS\t$evorate\t$treeness"
}
export -f summarize
echo -e "TREE\tTAXCOUNT\tALN\tALNLEN\tmean_BSS\tmedian_BSS\tevoRate\ttreeness" > ../gene_trees.$HMM.summarize_CDS.tsv
parallel -j $CPU summarize ::: $(find . -name "*.${TREEEXT}" ) >> ../gene_trees.$HMM.summarize_CDS.tsv
