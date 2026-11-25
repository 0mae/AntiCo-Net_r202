source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

cat data/tsv/taxa_summarized_50.tsv | \
  tail -n +2 | \
  cut -f 1,2 | \
  sed 's/\t/ /' | \
  xargs -I {} -P 12 -n 1 sh -c "bash scripts/bash/split_genomes_ogs_tidy_tsv.sh kegg {}"
