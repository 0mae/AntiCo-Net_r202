if [[ $1 == "cog" ]]; then

    if [[ $3 == "Domain" ]]; then
        cut -f 1,2,10,3 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tcog_id\n") - \
            > data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_$2.tsv
    elif [[ $3 == "Phylum" ]]; then
        cut -f 1,2,10,4 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tcog_id\n") - \
            > data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_$2.tsv
    elif [[ $3 == "Class" ]]; then
        cut -f 1,2,10,5 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tcog_id\n") - \
            > data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_$2.tsv
    elif [[ $3 == "Order" ]]; then
        cut -f 1,2,10,6 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tcog_id\n") - \
            > data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_$2.tsv
    elif [[ $3 == "Family" ]]; then
        cut -f 1,2,10,7 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tcog_id\n") - \
            > data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_$2.tsv
    elif [[ $3 == "Genus" ]]; then
        cut -f 1,2,10,8 data/tsv/proteins_raw_gtdbtax_cog2020_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tcog_id\n") - \
            > data/proteins_raw_gtdbtax_cog2020_tidy/proteins_raw_gtdbtax_cog2020_tidy_$2.tsv
    fi

elif [[ $1 == "kegg" ]]; then

    if [[ $3 == "Domain" ]]; then
        cut -f 1,2,10,3 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tko_id\n") - \
            > data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_$2.tsv
    elif [[ $3 == "Phylum" ]]; then
        cut -f 1,2,10,4 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tko_id\n") - \
            > data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_$2.tsv
    elif [[ $3 == "Class" ]]; then
        cut -f 1,2,10,5 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tko_id\n") - \
            > data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_$2.tsv
    elif [[ $3 == "Order" ]]; then
        cut -f 1,2,10,6 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tko_id\n") - \
            > data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_$2.tsv
    elif [[ $3 == "Family" ]]; then
        cut -f 1,2,10,7 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tko_id\n") - \
            > data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_$2.tsv
    elif [[ $3 == "Genus" ]]; then
        cut -f 1,2,10,8 data/tsv/proteins_raw_gtdbtax_ko_tidy.tsv | \
            grep -w $2 | \
            cat <(echo -ne "count\taccession\ttax\tko_id\n") - \
            > data/proteins_raw_gtdbtax_ko_tidy/proteins_raw_gtdbtax_ko_tidy_$2.tsv
    fi

fi
