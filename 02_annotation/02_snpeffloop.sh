#!/bin/bash

for i in $(ls neu-out);
do
	echo $i
	mkdir -p snpeff-out/${i}
	cd snpeff-out/${i}
	java -Xmx4g -jar ~/Applications_Downloaded/snpEff/snpEff.jar -v -cancer GRCh38.86 ~/Desktop/uni/diss/neu-out/${i}/PASS_filtered.vcf > ${i}.PASS_filtered.ann.vcf
	cat ${i}.PASS_filtered.ann.vcf | awk '$6 > 10' > ${i}.10_filtered.ann.vcf
	cat ${i}.PASS_filtered.ann.vcf | awk '$6 > 20' > ${i}.20_filtered.ann.vcf
	cat ${i}.PASS_filtered.ann.vcf | awk '$6 > 40' > ${i}.40_filtered.ann.vcf
	cd  ~/Desktop/uni/diss/

done
