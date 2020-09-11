#!/bin/bash

for i in $(ls snpeff-out);
	do cd snpeff-out/${i}
	sort -k6 -n  ${i}.PASS_filtered.ann.vcf | cut -f8 | awk -F "|" '{ print $4 }'  > gene_list_extract_PASS
	sort -k6 -n  ${i}.10_filtered.ann.vcf | cut -f8 | awk -F "|" '{ print $4 }'  > gene_list_extract_10
	sort -k6 -n  ${i}.20_filtered.ann.vcf | cut -f8 | awk -F "|" '{ print $4 }'  > gene_list_extract_20
	sort -k6 -n  ${i}.40_filtered.ann.vcf | cut -f8 | awk -F "|" '{ print $4 }'  > gene_list_extract_40
	cd ../../
done
