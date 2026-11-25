source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

seq -w 0 249 | xargs -I {} -P 36 -n 1 sh -c "Rscript --vanilla --slave scripts/r/search_result_filter.R data/diamond_blast_cog2020/proteins_raw.tsv.{} diamond data/diamond_blast_cog2020_filter/ 2>qsub_out/diamond_blast_cog2020_filter/{}.err"
