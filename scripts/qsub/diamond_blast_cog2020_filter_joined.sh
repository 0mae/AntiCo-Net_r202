source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

seq -w 0 249 | xargs -I {} -P 12 -n 1 sh -c "bash scripts/bash/join_filtered_cog2020.sh db/cog2020_fixed/cog-20.seqname_cogid.tsv data/diamond_blast_cog2020_filter/proteins_raw.filtered.tsv.{} cog > data/diamond_blast_cog2020_filter_joined/proteins_raw.filtered.cog2020.tsv.{} 2>qsub_out/diamond_blast_cog2020_filter_joined/{}.err"
