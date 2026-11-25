source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

cut -f 1 data/tsv/taxa_summarized_50.tsv | \
  tail -n +2 | \
  xargs -I {} -P 12 -n 1 sh -c "
    Rscript --vanilla --slave scripts/r/generate_contingency_tsv.R \
      cog \
      {} \
      data/contingency_tsv_cog2020/contingency_cog2020_{}.tsv \
      2>qsub_out/contingency_tsv_cog2020/{}.err
  "
