#!/bin/bash -l
#SBATCH --ntasks 96 --mem 96G --time 2:00:00 -p short -N 1 --out logs/align_parallel.log

module load miniconda3
conda activate /bigdata/stajichlab/shared/condaenv/phyling
#module load parallel -- installed in the conda env now
#module load biopython -- installed in the conda env now
#module load hmmer/3 --installed in the conda env
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
../PHYling_Unified/PHYling aln -c -q parallel
