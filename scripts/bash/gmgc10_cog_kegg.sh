#!/usr/bin/bash

if [[ $1 =~ COG[0-9][0-9][0-9][0-9] ]]; then 

        zgrep -w $1 raw/gmgc10/GMGC10_cog_kegg_env.tsv.gz | \
                #head | \ 
                cut --complement -f 2 | \
                cat <(zcat raw/gmgc10/GMGC10_cog_kegg_env.tsv.gz | head -n1 | cut --complement -f 2) - | \
                sed '1s/eggNOG_OGs/eggnog/; 1s/$/\tOGs/' | \
                sed "2,\$s/$/\t$1/" | \
                gzip > \
                data/gmgc10_cog/gmgc10_$1.tsv.gz

elif [[ $1 =~ K[0-9][0-9][0-9][0-9][0-9] ]]; then

        zgrep -w $1 raw/gmgc10/GMGC10_cog_kegg_env.tsv.gz | \
                #head | \
                cut --complement -f 3 | \
                cat <(zcat raw/gmgc10/GMGC10_cog_kegg_env.tsv.gz | head -n1 | cut --complement -f 3) - | \
                sed '1s/KEGG_ko/eggnog/; 1s/$/\tOGs/' | \
                sed "2,\$s/$/\t$1/" | \
                gzip > \
                data/gmgc10_kegg/gmgc10_$1.tsv.gz

fi
