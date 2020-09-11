#!/bin/bash

cd /path/to/project/neu-out

for i in $(ls);
do
	#filter
	cat ${i}/NeuSomatic.vcf | awk '$7!="FAILED"' > ${i}/filered.vcf
	cat ${i}/NeuSomatic.vcf | awk '$7=="PASS"' > ${i}/PASS_filtered.vcf
	cat ${i}/NeuSomatic.vcf | awk '$6 > 10' > ${i}/10_filtered.vcf
	cat ${i}/NeuSomatic.vcf | awk '$6 > 20' > ${i}/20_filtered.vcf
	cat ${i}/NeuSomatic.vcf | awk '$6 > 40' > ${i}/40_filtered.vcf

	#show number of samples in filtered
	echo "${i} no fails, only pass, >10, >20, >40"
	cut -f10 ${i}/filered.vcf | sort | uniq | wc -l
	cut -f10 ${i}/PASS_filtered.vcf | sort | uniq | wc -l
	cut -f10 ${i}/10_filtered.vcf | sort | uniq | wc -l
	cut -f10 ${i}/20_filtered.vcf | sort | uniq | wc -l
	cut -f10 ${i}/40_filtered.vcf | sort | uniq | wc -l;

done
