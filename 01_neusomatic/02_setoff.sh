#!/bin/bash

path="sbatch/tmp/"

#echo in brackets is due to read failing when input not terminated with newline, echo is used so the last line of the file can be read
(head -2 /home/path/sample_pairs.csv | tail -1 -; echo;) | while IFS=, read -r normal tumour rest
do
	sample="${normal}_${tumour}"
	echo ${sample}

	pren=$(sbatch sbatch/tmp/preneu_${sample}.sh)
	echo $pren
	pre=$(echo ${pren} | tail -c 9)
	echo $pre


	neun=$(sbatch --dependency=afterok:${pre} sbatch/tmp/neu_${sample}.sh)
	echo $neun
	neu=$(echo ${neun} | tail -c 9)
	echo $neu

	post=$(sbatch --dependency=afterok:${neu} sbatch/tmp/postneu_${sample}.sh)
	echo $post

done
 
<<< $(tail -44 /home/path/sample_pairs.csv | head -3 -)
