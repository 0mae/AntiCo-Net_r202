source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

seq -w 0 249 | xargs -I {} -P 12 -n 1 sh -c "Rscript --vanilla --slave scripts/r/assign_og.R data/diamond_blast_cog2020_filter_joined/proteins_raw.filtered.cog2020.tsv.{} cog data/diamond_blast_cog2020_filter_joined_assigned/ 2>qsub_out/diamond_blast_cog2020_filter_joined_assigned/{}.err"
