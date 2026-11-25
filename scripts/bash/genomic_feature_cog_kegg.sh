#!/usr/bin/bash

if [[ $1 =~ COG[0-9][0-9][0-9][0-9] ]]; then 

        grep -w $1 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
                #head | \ 
                cut -f 2,10 | \
                sed '1s/^/00__accession\tOG\n/' | \
		join -t "$(printf '\011')" -1 1 -2 1 \
			- \
			<(cat data/tsv/genomes_analyzed_{ar,bac}.tsv | sed '1s/^/00__/;/^accession/d' | sort -k 1,1 | cut -f 1,3-5,8-14,26,27,30,31,39,40,46,89-92,107-119) | \
		sed '1s/^00__//' | \
                gzip > \
                data/genomic_feature_cog/genomic_feature_$1.tsv.gz

elif [[ $1 =~ K[0-9][0-9][0-9][0-9][0-9] ]]; then

        grep -w $1 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
                #head | \
                cut -f 2,10 | \
                sed '1s/^/00__accession\tOG\n/' | \
		join -t "$(printf '\011')" -1 1 -2 1 \
			- \
			<(cat data/tsv/genomes_analyzed_{ar,bac}.tsv | sed '1s/^/00__/;/^accession/d' | sort -k 1,1 | cut -f 1,3-5,8-14,26,27,30,31,39,40,46,89-92,107-119) | \
		sed '1s/^00__//' | \
                gzip > \
                data/genomic_feature_kegg/genomic_feature_$1.tsv.gz

fi
