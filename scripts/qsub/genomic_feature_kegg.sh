cd $PBS_O_WORKDIR 

awk -F "\t" -v OFS="\t" '$5 == "contrapositive" && $3 >= 8 {print}' data/tsv/network_edge_stats_kegg.tsv | \
	cut -f 1,2 | \
	perl -ne 's/\t/\n/;s/ko://g;print' | \
	sort | uniq | \
	#head -n2 | \
	xargs -I {} -P 12 -n 1 sh -c "bash scripts/bash/genomic_feature_cog_kegg.sh {}"
