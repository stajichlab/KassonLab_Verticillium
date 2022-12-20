#!/usr/bin/bash -l
mkdir -p input
for SAMPLE in $(tail -n +2 samples.csv  | cut -d, -f1)
do
	ln -s /bigdata/stajichlab/shared/projects/SeqData/Kasson/2022_12_Verticillium_Ento/${SAMPLE}_*/${SAMPLE}*.gz input/
done

