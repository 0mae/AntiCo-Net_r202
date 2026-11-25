seq -w 0 249 | \
  xargs -I {} sh -c "
    cat data/diamond_blast_kegg_filter_joined_assigned/proteins_raw.filtered.cog2020.assigned.tsv.{} | sed -e '1d'
  " | \
  cat <(
    head -n1 data/diamond_blast_kegg_filter_joined_assigned/proteins_raw.filtered.cog2020.assigned.tsv.000
  ) - | \
  sed '1s/^query/00_prodigal_acc/' | \
  sort -k 1,1 | \
  join -t "$(printf '\011')" -1 2 -2 1 <(cut --complement -f 6 data/tsv/proteins_raw_gtdbtax.tsv | sed '1s/prodigal_acc/00_prodigal_acc/' | sort -k 2,2) - | \
  sed '1s/00_prodigal_acc/gene_id/' | \
  awk -F "\t" -v OFS="\t" '{print $6,$0}' | \
  sed 's/_/\t/' | \
  sed '1s/prodigal/00000/;1s/ID/00000/' | \
  sort -k 4,4 -k 1,1n -k 2,2n | \
  cut --complement -f 1,2 \
  > data/tsv/proteins_raw_gtdbtax_ko.tsv
