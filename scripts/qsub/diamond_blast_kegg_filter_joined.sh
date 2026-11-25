source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

seq -w 0 249 | xargs -I {} -P 12 -n 1 sh -c "
    bash scripts/bash/join_filtered_cog2020.sh \
      db/kegg_fixed/kegg_geneid_koid.tsv \
      data/diamond_blast_kegg_filter/proteins_raw.filtered.tsv.{} \
      kegg \
      > data/diamond_blast_kegg_filter_joined/proteins_raw.filtered.cog2020.tsv.{} \
      2>qsub_out/diamond_blast_kegg_filter_joined/{}.err
  "
