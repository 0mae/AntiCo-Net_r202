cat data/tsv/genomes_analyzed_ar.tsv <(tail -n +2 data/tsv/genomes_analyzed_bac.tsv) | \
    cut -f 1,111-117 | \
    sort -k 1,1 | \
    join -t "$(printf '\011')" -1 1 -2 1 -a 1 <(
        sort -k 1,1 data/tsv/proteins_raw.tsv | perl -ne 's/\t$//; print'
    ) - \
    > data/tsv/proteins_raw_gtdbtax.tsv
