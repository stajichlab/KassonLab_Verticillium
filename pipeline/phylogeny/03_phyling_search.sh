#!/bin/bash -l
#SBATCH --ntasks 2 --mem 24G --time 2:00:00 -p short -N 1 --out logs/phyling_init_search.log
module load miniconda3
conda activate /bigdata/stajichlab/shared/condaenv/phyling

pushd Phylogeny
if [ ! -f config.txt ]; then
    echo "Need config.txt for PHYling"
    exit
fi

source config.txt
if [ ! -z $PREFIX ]; then
    rm -rf aln/$PREFIX
fi
# probably should check to see if allseq is newer than newest file in the folder?
echo " I will remove prefix.tab to make sure it is regenerated"
pwd
rm prefix.tab
../PHYling_Unified/PHYling init
../PHYling_Unified/PHYling search -q slurm
#../PHYling_Unified/PHYling search -q parallel
