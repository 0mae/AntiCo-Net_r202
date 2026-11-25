#PBS -S /bin/bash
#PBS -l ncpus=40

cd $HOME/crest/cpnet

cat <(tail -n +2 data/tsv/network_edge_stats_cog.tsv) <(tail -n +2 data/tsv/network_edge_stats_kegg.tsv) | \
	awk -F "\t" -v OFS="\t" '$5 == "contrapositive" {print}' | \
	cut -f 1,2 | \
	perl -ne 's/\t/\n/;s/ko://g;print' | \
	sort | uniq | \
	xargs -I {} -P 40 -n 1 sh -c "bash scripts/bash/gmgc10_cog_kegg.sh {}"
